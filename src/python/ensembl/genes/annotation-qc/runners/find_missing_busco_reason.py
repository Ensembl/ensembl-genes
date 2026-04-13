"""
BUSCO Missing Annotation Scanner (advanced)

For BUSCOs complete in reference but missing in test:
	• map BUSCO → reference transcript + coordinates
	• scan multiple GTF/GFF files
	• compute:
		- presence
		- length coverage
		- intron count
		- multi-hit reporting

Output:
	missing_busco_annotation_presence.csv
"""

import argparse
from collections import defaultdict
import pandas as pd


# ---------------------------------------------------
# BUSCO
# ---------------------------------------------------

def parse_busco_table(path):
	buscos = {}
	with open(path) as f:
		for line in f:
			if line.startswith("#") or not line.strip():
				continue
			parts = line.strip().split("\t")
			buscos[parts[0]] = {
				"status": parts[1],
				"sequence": parts[2] if len(parts) > 2 else None
			}
	return buscos


# ---------------------------------------------------
# Attributes
# ---------------------------------------------------

def parse_attributes(attr_string):
	attrs = {}

	if "=" in attr_string:  # GFF3
		for a in attr_string.split(";"):
			if "=" in a:
				k, v = a.split("=", 1)
				attrs[k.strip()] = v.strip()
	else:  # GTF
		for a in attr_string.split(";"):
			a = a.strip()
			if not " " in a:
				continue
			k, v = a.split(" ", 1)
			attrs[k.strip()] = v.strip('"')

	return attrs


# ---------------------------------------------------
# Annotation parsing
# ---------------------------------------------------

def parse_annotation(path, cds_only=False):
	"""
	cds_only=True  →  transcript start/end are set to CDS bounds only;
					  transcripts with no CDS lines are dropped entirely.
					  Use this for the *reference* annotation when the GFF
					  contains UTR features that would otherwise inflate
					  the reference span.
	cds_only=False →  original behaviour (start/end = full genomic extent).
	"""
	transcripts = defaultdict(lambda: {
		"chrom": None,
		"strand": None,
		"start": 10**12,
		"end": 0,
		"exons": [],
		"cds_exons": [],
	})

	id_map = {}

	with open(path) as f:
		for line in f:
			if line.startswith("#"):
				continue

			parts = line.strip().split("\t")
			if len(parts) < 9:
				continue

			chrom, source, feature, start, end, score, strand, phase, attrs = parts
			start, end = int(start), int(end)

			attr = parse_attributes(attrs)

			# GTF always has transcript_id on every line → use it directly.
			# GFF3 sub-features (exon, CDS, UTR …) carry their own ID but
			# link back to the transcript via Parent=.  Parent can be a
			# comma-separated list (shared exons) — handle all of them.
			# Fall back to ID only when there is no Parent (i.e. the line
			# IS the transcript/mRNA/gene itself).
			if attr.get("transcript_id"):
				parent_ids = [attr["transcript_id"]]          # GTF
			elif attr.get("Parent"):
				parent_ids = [p.strip() for p in attr["Parent"].split(",")]  # GFF3 sub-feature(s)
			elif attr.get("ID"):
				parent_ids = [attr["ID"]]                     # GFF3 transcript/mRNA/gene
			elif attr.get("protein_id"):
				parent_ids = [attr["protein_id"]]
			else:
				parent_ids = []

			pid = attr.get("protein_id")

			if not parent_ids:
				continue

			for tid in parent_ids:

				t = transcripts[tid]

				t["chrom"] = chrom
				t["strand"] = strand
				t["start"] = min(t["start"], start)
				t["end"] = max(t["end"], end)

				if feature.lower() == "exon":
					t["exons"].append((start, end))

				# track CDS separately so we can drop UTRs without collapsing structure
				if feature.lower() == "cds":
					t["cds_exons"].append((start, end))

				id_map[tid] = tid
				if pid:
					id_map[pid] = tid

	# --- replace exons with CDS exons if requested (preserve structure, just drop UTRs) ---
	if cds_only:
		for tid in list(transcripts.keys()):
			t = transcripts[tid]

			# drop transcripts that have no CDS
			if not t["cds_exons"]:
				del transcripts[tid]
				continue

			# keep CDS exon structure (do NOT collapse)
			t["exons"] = sorted(t["cds_exons"])
			t["start"] = min(s for s, _ in t["exons"])
			t["end"] = max(e for _, e in t["exons"])

	return transcripts, id_map


# ---------------------------------------------------
# Map BUSCO → reference transcript
# ---------------------------------------------------

def map_buscos(buscos, transcripts, id_map):

	for bid, info in buscos.items():
		seq = info["sequence"]

		if seq in id_map:
			tid = id_map[seq]
			t = transcripts[tid]

			info["ref_transcript"] = tid
			info["chrom"] = t["chrom"]
			info["start"] = t["start"]
			info["end"] = t["end"]
			info["ref_exons"] = sorted(t["exons"])   # preserve exon structure
		else:
			info["ref_transcript"] = None
			info["chrom"] = info["start"] = info["end"] = None
			info["ref_exons"] = []

	return buscos


# ---------------------------------------------------
# Overlap helpers
# ---------------------------------------------------

def overlap(a_start, a_end, b_start, b_end):
	return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def clip_exons(exons, region_start, region_end):
	"""Clip a list of exons to a region window and drop any that fall outside entirely."""
	clipped = []
	for es, ee in exons:
		cs = max(es, region_start)
		ce = min(ee, region_end)
		if cs <= ce:
			clipped.append((cs, ce))
	return sorted(clipped)


def classify_region(chrom, start, end, ref_exons, transcripts):
	"""
	Compare each candidate's exon structure against the reference.

	Behaviour:
	  • exact same clipped exon structure  → PRESENT_FULL
	  • otherwise choose closest transcript by exon coverage + intron similarity
		→ PRESENT_SHORT
	  • no overlap → ABSENT

	Coverage is computed on *exonic bases only* (not genomic span).
	Intron suffix is only used for SHORT results.
	"""

	if start is None or end is None:
		return "ABSENT", 0, 0, 0, ""

	# ---- reference exon length (exonic only) ----
	ref_exon_set = tuple(sorted(ref_exons))
	ref_len = sum(e - s + 1 for s, e in ref_exon_set)

	if ref_len == 0:
		return "ABSENT", 0, 0, 0, ""

	ref_introns = max(0, len(ref_exon_set) - 1)

	full_matches = []
	candidates = []  # (tid, cov, introns, exon_overlap)

	for tid, t in transcripts.items():

		if t["chrom"] != chrom:
			continue

		# clip to BUSCO region
		clipped = clip_exons(t["exons"], start, end)
		if not clipped:
			continue

		introns = max(0, len(clipped) - 1)

		# ---- exon-based overlap coverage ----
		exon_overlap = 0
		for rs, re in ref_exon_set:
			for ts, te in clipped:
				exon_overlap += overlap(rs, re, ts, te)

		cov = round(exon_overlap / ref_len * 100, 1)

		candidates.append((tid, cov, introns, exon_overlap))

		# exact structure match
		# exact structure match inside region
		if tuple(clipped) == ref_exon_set:

			# check if transcript extends outside BUSCO region
			raw_start = min(s for s, e in t["exons"])
			raw_end = max(e for s, e in t["exons"])

			if raw_start < start or raw_end > end:
				full_matches.append((tid, cov, introns, exon_overlap, "LONGER"))
			else:
				full_matches.append((tid, cov, introns, exon_overlap, "FULL"))

	if not candidates:
		return "ABSENT", 0, 0, 0, ""

	hits = len(candidates)
	names = ",".join([c[0] for c in candidates])

	# ---- exact structure match: always FULL (no intron suffix) ----

	if full_matches:
		# choose closest in case multiple (coverage first, intron diff second)
		best = max(full_matches, key=lambda x: (x[1], -abs(x[2] - ref_introns)))
		tid, cov, introns, _, match_type = best
		names = tid

		if match_type == "LONGER":
			return "PRESENT_LONGER", cov, introns, len(candidates), names
		else:
			return "PRESENT_FULL", cov, introns, len(candidates), names

	# ---- otherwise choose closest candidate ----
	best = max(candidates, key=lambda x: (x[1], -abs(x[2] - ref_introns)))
	tid, cov, introns, _ = best
	names = tid  # only the transcript that determined the status

	status = "PRESENT_SHORT"
	if introns > 0:
		status += "_INTRONS"

	return status, cov, introns, len(candidates), names


# ---------------------------------------------------
# Main
# ---------------------------------------------------

def analyze(buscos_ref, buscos_test, ref_anno, annos, names, cds_only=False):

	complete_ref = {b for b, i in buscos_ref.items() if i["status"] in ["Complete", "Duplicated"]}
	complete_test = {b for b, i in buscos_test.items() if i["status"] in ["Complete", "Duplicated"]}

	missing = sorted(complete_ref - complete_test)

	print(f"Missing BUSCOs: {len(missing)}")

	# cds_only=True strips UTRs from the reference so that the reference
	# span used for coverage / exact-match comparison is CDS-only.
	ref_tx, ref_map = parse_annotation(ref_anno, cds_only=cds_only)
	buscos_ref = map_buscos(buscos_ref, ref_tx, ref_map)

	# Test annotation (last in files)
	test_tx, _ = parse_annotation(args.test_annotation, cds_only=args.cds_only)

	# Optional layers (all except test)
	layer_files = files[:-1]
	layer_names = names[:-1]
	parsed_layers = [parse_annotation(f, cds_only=False)[0] for f in layer_files]
	parsed_layers.append(test_tx)  # add test annotation at the end

	rows = []

	for bid in missing:

		info = buscos_ref[bid]
		chrom = info["chrom"]
		start = info["start"]
		end = info["end"]
		# clip stored ref exons to [start, end] once (handles cds_only trimming)
		ref_exons = clip_exons(info["ref_exons"], start, end) if chrom else []

		row = {
			"BUSCO_ID": bid,
			"Reference_Transcript": info["ref_transcript"],
			"Chrom": chrom,
			"Start": start,
			"End": end
		}

		for name, transcripts in zip(names, parsed_layers):

			if chrom:
				s, cov, intr, hits, tnames = classify_region(chrom, start, end, ref_exons, transcripts)
			else:
				s, cov, intr, hits, tnames = "ABSENT", 0, 0, 0, ""

			row[f"{name}_status"] = s
			row[f"{name}_coverage_pct"] = cov
			row[f"{name}_introns"] = intr
			row[f"{name}_hits"] = hits
			row[f"{name}_transcripts"] = tnames

		# --- Add hierarchical error column ---
		def get(row, key):
			return row.get(key, "")

		if get(row, "test_status") == "PRESENT_FULL":
			row["status_check"] = "non_can_hit"

		elif (get(row, "miniprot_ortho_sel_status") == "PRESENT_FULL" and
		      get(row, "protein_sel_status") != "PRESENT_FULL"):
			row["status_check"] = "miniprot_hit"

		elif (get(row, "protein_sel_status") == "PRESENT_FULL" or
		      get(row, "uniprot_sel_status") == "PRESENT_FULL"):
			row["status_check"] = "genebuild_error"

		elif ((get(row, "protein_status") == "PRESENT_FULL" and
		       get(row, "protein_sel_status") != "PRESENT_FULL") or
		      (get(row, "uniprot_status") == "PRESENT_FULL" and
		       get(row, "uniprot_sel_status") != "PRESENT_FULL") or
		      (get(row, "miniprot_ortho_status") == "PRESENT_FULL" and
		       get(row, "miniprot_ortho_sel_status") != "PRESENT_FULL")):
			row["status_check"] = "select_error"

		elif ((get(row, "protein_status") == "PRESENT_SHORT_INTRONS" and
		       get(row, "protein_sel_status") == "PRESENT_SHORT_INTRONS") or
		      (get(row, "uniprot_status") == "PRESENT_SHORT_INTRONS" and
		       get(row, "uniprot_sel_status") == "PRESENT_SHORT_INTRONS") or
		      (get(row, "miniprot_ortho_status") == "PRESENT_SHORT_INTRONS" and
		       get(row, "miniprot_ortho_sel_status") == "PRESENT_SHORT_INTRONS")):
			row["status_check"] = "short_intron_error"

		else:
			row["status_check"] = ""

		rows.append(row)

	df = pd.DataFrame(rows)
	df.to_csv("missing_busco_annotation_presence.csv", index=False)
	# --- Summary of test_status ---
	if "test_status" in df.columns:
		summary = (
			df["test_status"]
			.value_counts(dropna=False)
			.rename_axis("test_status")
			.reset_index(name="count")
		)

		print("\n=== TEST STATUS SUMMARY ===")
		print(summary.to_string(index=False))
	else:
		print("\nNo test_status column found.")


# ---------------------------------------------------
# CLI
# ---------------------------------------------------

if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument("--busco_ref", required=True)
	parser.add_argument("--busco_test", required=True)
	parser.add_argument("--ref_annotation", required=True)

	parser.add_argument("--transcriptomic", required=False)
	parser.add_argument("--protein", required=False)
	parser.add_argument("--protein_sel", required=False)
	parser.add_argument("--uniprot", required=False)
	parser.add_argument("--uniprot_sel", required=False)
	parser.add_argument("--miniprot_ortho", required=False)
	parser.add_argument("--miniprot_ortho_sel", required=False)
	parser.add_argument("--test_annotation", required=True)

	parser.add_argument(
		"--cds_only", action="store_true", default=False,
		help="Apply CDS-only mode to ALL annotations (reference + test). "
		     "Removes UTRs and keeps only CDS exon structure."
	)

	args = parser.parse_args()

	buscos_ref = parse_busco_table(args.busco_ref)
	buscos_test = parse_busco_table(args.busco_test)

	optional_inputs = [
		("transcriptomic", args.transcriptomic),
		("protein", args.protein),
		("protein_sel", args.protein_sel),
		("uniprot", args.uniprot),
		("uniprot_sel", args.uniprot_sel),
		("miniprot_ortho", args.miniprot_ortho),
		("miniprot_ortho_sel", args.miniprot_ortho_sel),
	]

	# keep only provided ones
	filtered = [(n, f) for n, f in optional_inputs if f]

	# always include test annotation at the end
	filtered.append(("test", args.test_annotation))

	names = [n for n, _ in filtered]
	files = [f for _, f in filtered]

	print("Using annotation layers:", ", ".join(names))

	analyze(buscos_ref, buscos_test, args.ref_annotation, files, names, cds_only=args.cds_only)