#!/usr/bin/env python3
"""
compare_pep_fastas.py - validate a generated protein FASTA against a reference.

Compares record sets, headers, sequences, record order and decompressed byte
content of two pep FASTA files. Either input may be plain text or gzip; the
comparison is pure Python, so no external diff or shell features are needed.

Exit codes:
    0  exact parity, including byte-identical decompressed content
    1  records match but the byte content differs (formatting only)
    2  biological, header or ordering mismatch
    3  file read or parse error
"""

import argparse
import gzip
import json
import os
import sys
from datetime import datetime, timezone

EXIT_EXACT = 0
EXIT_FORMATTING = 1
EXIT_MISMATCH = 2
EXIT_ERROR = 3

CONTEXT_FLANK = 10
BYTE_CHUNK = 1 << 20


class FastaError(Exception):
    """Unreadable or malformed FASTA input."""


def open_fasta(path, mode='rt'):
    """Open a FASTA file transparently for gzip and plain text."""
    try:
        if path.endswith('.gz'):
            return gzip.open(path, mode)
        return open(path, mode)
    except OSError as exc:
        raise FastaError("cannot open {0}: {1}".format(path, exc))


def parse_fasta(path):
    """Yield (record_id, header, sequence) for each record in a FASTA file."""
    handle = open_fasta(path, 'rt')
    try:
        header = None
        chunks = []
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header.split(None, 1)[0], header, ''.join(chunks)
                header = line[1:]
                chunks = []
            elif header is None:
                if line.strip():
                    raise FastaError(
                        "{0}: sequence data before the first header".format(path))
            else:
                chunks.append(line)
        if header is not None:
            yield header.split(None, 1)[0], header, ''.join(chunks)
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise FastaError("cannot read {0}: {1}".format(path, exc))
    finally:
        handle.close()


def load_fasta(path):
    """
    Read a FASTA file into ({record_id: (header, sequence)}, [record_id, ...]).

    Duplicate record IDs are reported rather than silently collapsed.
    """
    records = {}
    order = []
    duplicates = []
    for record_id, header, sequence in parse_fasta(path):
        if record_id in records:
            duplicates.append(record_id)
        records[record_id] = (header, sequence)
        order.append(record_id)
    return records, order, duplicates


def bytes_identical(path_a, path_b):
    """Compare the decompressed byte content of two FASTA files."""
    handle_a = open_fasta(path_a, 'rb')
    handle_b = open_fasta(path_b, 'rb')
    try:
        while True:
            chunk_a = handle_a.read(BYTE_CHUNK)
            chunk_b = handle_b.read(BYTE_CHUNK)
            if chunk_a != chunk_b:
                return False
            if not chunk_a:
                return True
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise FastaError("cannot read inputs for byte comparison: {0}".format(exc))
    finally:
        handle_a.close()
        handle_b.close()


def first_difference(expected, observed):
    """Index of the first differing character, or None if one is a prefix."""
    for i, (a, b) in enumerate(zip(expected, observed)):
        if a != b:
            return i
    if len(expected) != len(observed):
        return min(len(expected), len(observed))
    return None


def context(sequence, position, flank=CONTEXT_FLANK):
    """A short window of sequence around a position, for the report."""
    if position is None:
        return ''
    start = max(0, position - flank)
    return sequence[start:position + flank + 1]


def ordering_mismatches(expected_order, observed_order, shared):
    """
    Count positions where the two files disagree on record order.

    Only records present in both files are considered, so a missing record
    does not inflate the ordering count on its own.
    """
    expected_shared = [r for r in expected_order if r in shared]
    observed_shared = [r for r in observed_order if r in shared]
    mismatches = sum(
        1 for a, b in zip(expected_shared, observed_shared) if a != b
    )
    mismatches += abs(len(expected_shared) - len(observed_shared))
    return mismatches, expected_shared, observed_shared


def compare(expected_path, observed_path):
    """Compare two FASTA files and return (summary dict, report rows)."""
    expected, expected_order, expected_dupes = load_fasta(expected_path)
    observed, observed_order, observed_dupes = load_fasta(observed_path)

    expected_ids = set(expected)
    observed_ids = set(observed)
    missing = sorted(expected_ids - observed_ids)
    extra = sorted(observed_ids - expected_ids)
    shared = expected_ids & observed_ids

    order_count, expected_shared, observed_shared = ordering_mismatches(
        expected_order, observed_order, shared)

    rows = []
    for record_id in missing:
        rows.append({
            'category': 'missing',
            'record_id': record_id,
            'detail': 'present in expected, absent from observed',
            'expected': expected[record_id][0],
            'observed': '',
        })
    for record_id in extra:
        rows.append({
            'category': 'extra',
            'record_id': record_id,
            'detail': 'present in observed, absent from expected',
            'expected': '',
            'observed': observed[record_id][0],
        })

    header_count = 0
    sequence_count = 0
    for record_id in expected_order:
        if record_id not in shared:
            continue
        expected_header, expected_seq = expected[record_id]
        observed_header, observed_seq = observed[record_id]

        if expected_header != observed_header:
            header_count += 1
            position = first_difference(expected_header, observed_header)
            rows.append({
                'category': 'header_mismatch',
                'record_id': record_id,
                'detail': 'first_diff_char={0}'.format(position),
                'expected': expected_header,
                'observed': observed_header,
            })

        if expected_seq != observed_seq:
            sequence_count += 1
            position = first_difference(expected_seq, observed_seq)
            rows.append({
                'category': 'sequence_mismatch',
                'record_id': record_id,
                'detail': 'first_diff_aa_pos={0} expected_len={1} observed_len={2}'.format(
                    'NA' if position is None else position + 1,
                    len(expected_seq), len(observed_seq)),
                'expected': context(expected_seq, position),
                'observed': context(observed_seq, position),
            })

    for position, (want, got) in enumerate(zip(expected_shared, observed_shared)):
        if want != got:
            rows.append({
                'category': 'order_mismatch',
                'record_id': want,
                'detail': 'shared_record_index={0}'.format(position),
                'expected': want,
                'observed': got,
            })
            break

    for record_id in expected_dupes:
        rows.append({
            'category': 'duplicate_in_expected', 'record_id': record_id,
            'detail': 'record id appears more than once', 'expected': '', 'observed': '',
        })
    for record_id in observed_dupes:
        rows.append({
            'category': 'duplicate_in_observed', 'record_id': record_id,
            'detail': 'record id appears more than once', 'expected': '', 'observed': '',
        })

    records_match = not (missing or extra or header_count or sequence_count
                         or order_count or expected_dupes or observed_dupes)
    identical_bytes = bytes_identical(expected_path, observed_path) if records_match else False

    if records_match and identical_bytes:
        status, exit_code = 'exact_parity', EXIT_EXACT
    elif records_match:
        status, exit_code = 'formatting_difference', EXIT_FORMATTING
    else:
        status, exit_code = 'mismatch', EXIT_MISMATCH

    summary = {
        'expected_file': os.path.abspath(expected_path),
        'observed_file': os.path.abspath(observed_path),
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'expected_record_count': len(expected_order),
        'observed_record_count': len(observed_order),
        'missing_record_count': len(missing),
        'extra_record_count': len(extra),
        'header_mismatch_count': header_count,
        'sequence_mismatch_count': sequence_count,
        'ordering_mismatch_count': order_count,
        'duplicate_expected_count': len(expected_dupes),
        'duplicate_observed_count': len(observed_dupes),
        'missing_record_ids': missing[:100],
        'extra_record_ids': extra[:100],
        'decompressed_bytes_identical': identical_bytes,
        'status': status,
        'passed': status == 'exact_parity',
        'exit_code': exit_code,
    }
    return summary, rows


def write_report(path, rows):
    """Write the mismatch-only TSV report."""
    with open(path, 'w') as fh:
        fh.write("category\trecord_id\tdetail\texpected\tobserved\n")
        for row in rows:
            fh.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                row['category'], row['record_id'], row['detail'],
                row['expected'].replace('\t', ' '),
                row['observed'].replace('\t', ' '),
            ))


def print_summary(summary, out=sys.stderr):
    print("expected file            {0}".format(summary['expected_file']), file=out)
    print("observed file            {0}".format(summary['observed_file']), file=out)
    print("expected records         {0:,}".format(summary['expected_record_count']), file=out)
    print("observed records         {0:,}".format(summary['observed_record_count']), file=out)
    print("missing records          {0:,}".format(summary['missing_record_count']), file=out)
    print("extra records            {0:,}".format(summary['extra_record_count']), file=out)
    print("header mismatches        {0:,}".format(summary['header_mismatch_count']), file=out)
    print("sequence mismatches      {0:,}".format(summary['sequence_mismatch_count']), file=out)
    print("ordering mismatches      {0:,}".format(summary['ordering_mismatch_count']), file=out)
    print("duplicate ids (expected) {0:,}".format(summary['duplicate_expected_count']), file=out)
    print("duplicate ids (observed) {0:,}".format(summary['duplicate_observed_count']), file=out)
    print("bytes identical          {0}".format(summary['decompressed_bytes_identical']), file=out)
    print("status                   {0}".format(summary['status']), file=out)


class _ArgumentParser(argparse.ArgumentParser):
    """argparse exits 2 on a usage error; 2 means mismatch here."""

    def error(self, message):
        self.print_usage(sys.stderr)
        self.exit(EXIT_ERROR, "{0}: error: {1}\n".format(self.prog, message))


def main(argv=None):
    parser = _ArgumentParser(
        prog='compare_pep_fastas.py',
        description='Compare a generated pep FASTA against a reference dump.',
    )
    parser.add_argument('--expected', required=True,
                        help='reference FASTA (.fa or .fa.gz)')
    parser.add_argument('--observed', required=True,
                        help='generated FASTA (.fa or .fa.gz)')
    parser.add_argument('--report', help='TSV path for mismatch rows')
    parser.add_argument('--summary', help='JSON path for the summary')
    args = parser.parse_args(argv)

    try:
        summary, rows = compare(args.expected, args.observed)
        if args.report:
            write_report(args.report, rows)
        if args.summary:
            with open(args.summary, 'w') as fh:
                json.dump(summary, fh, indent=2, sort_keys=True)
                fh.write('\n')
    except FastaError as exc:
        print("error: {0}".format(exc), file=sys.stderr)
        return EXIT_ERROR
    except OSError as exc:
        print("error: {0}".format(exc), file=sys.stderr)
        return EXIT_ERROR

    print_summary(summary)
    return summary['exit_code']


if __name__ == '__main__':
    sys.exit(main())
