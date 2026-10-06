"""
Local dashboard server (standard library only).

Binds to 127.0.0.1 by default and serves the static single-page app plus a small
JSON/CSV API over one or more dashboard datasets (read-only SQLite). No hosted
services, no comparison work at request time.

Endpoints (all GET):
    /api/experiments                         experiments served
    /api/meta?exp=                           configuration, policies, runs, validation, metric dictionary
    /api/overview?exp=&reference=[&format=csv]
    /api/genes?exp=&reference=&policy=&...   reference-gene table (+format=csv)
    /api/predictions?exp=&reference=&policy=&run=&...   prediction table (+format=csv)
    /api/gene?exp=&reference=&policy=&gene_id=
    /api/locus?exp=&reference=&chrom=&start=&end=&policy=&runs=&gene_id=&isoforms=
    /api/chromosomes?exp=&reference=
    /api/document?exp=&index=
"""

from __future__ import annotations

import csv
import io
import json
import mimetypes
import sqlite3
import threading
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from ensembl.genes.annotation_qc.dashboard import queries

STATIC = Path(__file__).with_name("static")


class Datasets:
    """Experiment id -> SQLite path; one read-only connection per thread and dataset."""

    def __init__(self, dbs: dict[str, Path]):
        self.dbs = dbs
        self.local = threading.local()
        self._meta: dict[str, dict] = {}

    def con(self, exp: str) -> sqlite3.Connection:
        if exp not in self.dbs:
            raise queries.QueryError(f"unknown experiment {exp!r}")
        cache = getattr(self.local, "cons", None)
        if cache is None:
            cache = self.local.cons = {}
        if exp not in cache:
            cache[exp] = queries.connect(self.dbs[exp])
        return cache[exp]

    def meta(self, exp: str) -> dict:
        if exp not in self._meta:
            self._meta[exp] = queries.meta(self.con(exp))
        return self._meta[exp]


def _csv(
    rows: list[dict], flatten_outcomes: bool = False, runs: list[str] | None = None
) -> str:
    out = io.StringIO()
    if not rows:
        return ""
    if flatten_outcomes:
        flat = []
        for r in rows:
            base = {k: v for k, v in r.items() if k != "outcomes"}
            for run in runs or sorted(r["outcomes"]):
                o = r["outcomes"].get(run) or {}
                for key in (
                    "cds_status",
                    "cds_exact",
                    "cds_overlap",
                    "cds_matched_id",
                    "counterpart_count",
                    "any_pair_chain",
                ):
                    base[f"{run}.{key}"] = o.get(key, "NA" if not o else "")
            flat.append(base)
        rows = flat
    writer = csv.DictWriter(
        out, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue()


def make_handler(datasets: Datasets):
    class Handler(BaseHTTPRequestHandler):
        server_version = "annotation-qc-dashboard"

        def log_message(self, fmt, *args):  # quieter default logging
            if getattr(self.server, "verbose", False):
                super().log_message(fmt, *args)

        def _send(
            self,
            status: int,
            body: bytes,
            content_type: str,
            filename: str | None = None,
        ):
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            if filename:
                self.send_header(
                    "Content-Disposition", f'attachment; filename="{filename}"'
                )
            self.end_headers()
            self.wfile.write(body)

        def _json(self, data, status=200):
            self._send(
                status, json.dumps(data, default=str).encode(), "application/json"
            )

        def do_GET(self):  # noqa: N802
            url = urlparse(self.path)
            p = {k: v[-1] for k, v in parse_qs(url.query).items()}
            try:
                if url.path in ("/", "/index.html"):
                    return self._static("index.html")
                if url.path.startswith("/static/"):
                    return self._static(url.path[len("/static/") :])
                if url.path == "/api/experiments":
                    return self._json(
                        [
                            {
                                "id": e,
                                **{
                                    k: datasets.meta(e)["experiment"].get(k)
                                    for k in ("title", "description", "synthetic")
                                },
                            }
                            for e in datasets.dbs
                        ]
                    )
                exp = p.get("exp") or next(iter(datasets.dbs))
                con = datasets.con(exp)
                if url.path == "/api/meta":
                    meta = dict(datasets.meta(exp))
                    meta["context_documents"] = [
                        {
                            "index": i,
                            "name": d["name"],
                            "path": d["path"],
                            "available": d["text"] is not None,
                        }
                        for i, d in enumerate(meta.get("context_documents") or [])
                    ]
                    return self._json(meta)
                if url.path == "/api/document":
                    docs = datasets.meta(exp).get("context_documents") or []
                    doc = docs[int(p.get("index", -1))]
                    return self._json(doc)
                if url.path == "/api/overview":
                    data = queries.overview(con, p["reference"])
                    if p.get("format") == "csv":
                        return self._send(
                            200,
                            _csv(data["metrics"]).encode(),
                            "text/tab-separated-values",
                            f"{exp}_{p['reference']}_metrics.tsv",
                        )
                    return self._json(data)
                if url.path == "/api/chromosomes":
                    return self._json(queries.chromosomes(con, p["reference"]))
                if url.path == "/api/genes":
                    as_csv = p.get("format") == "csv"
                    data = queries.genes(con, p, csv=as_csv)
                    if as_csv:
                        runs = [
                            r for r in (p.get("runs") or "").split(",") if r
                        ] or None
                        return self._send(
                            200,
                            _csv(data["rows"], True, runs).encode(),
                            "text/tab-separated-values",
                            f"{exp}_{p['reference']}_{p['policy']}_reference_genes.tsv",
                        )
                    return self._json(data)
                if url.path == "/api/predictions":
                    as_csv = p.get("format") == "csv"
                    data = queries.query_genes(con, p, csv=as_csv)
                    if as_csv:
                        return self._send(
                            200,
                            _csv(data["rows"]).encode(),
                            "text/tab-separated-values",
                            f"{exp}_{p['run']}_{p['policy']}_predictions.tsv",
                        )
                    return self._json(data)
                if url.path == "/api/gene":
                    return self._json(
                        queries.gene_detail(
                            con, p["reference"], p["policy"], p["gene_id"]
                        )
                    )
                if url.path == "/api/locus":
                    return self._json(queries.locus(con, p))
                return self._json({"error": f"unknown path {url.path}"}, 404)
            except (queries.QueryError, KeyError, IndexError) as error:
                return self._json({"error": f"{type(error).__name__}: {error}"}, 400)
            except Exception as error:  # noqa: BLE001
                traceback.print_exc()
                return self._json({"error": f"internal error: {error}"}, 500)

        def _static(self, rel: str):
            path = (STATIC / rel).resolve()
            if STATIC.resolve() not in path.parents or not path.is_file():
                return self._json({"error": "not found"}, 404)
            ctype = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
            if path.suffix == ".js":
                ctype = "text/javascript"
            return self._send(
                200,
                path.read_bytes(),
                ctype + ("; charset=utf-8" if ctype.startswith("text") else ""),
            )

    return Handler


def serve(
    dbs: dict[str, Path], host: str = "127.0.0.1", port: int = 8765, verbose=False
) -> None:
    missing = {e: str(p) for e, p in dbs.items() if not Path(p).exists()}
    if missing:
        raise SystemExit(
            f"Dashboard dataset not built for {missing}; run `annotation-qc benchmark build --config …` first"
        )
    httpd = ThreadingHTTPServer((host, port), make_handler(Datasets(dbs)))
    httpd.verbose = verbose
    print(
        f"Annotation QC dashboard: http://{host}:{port}/  (Ctrl+C to stop)", flush=True
    )
    if host not in ("127.0.0.1", "localhost", "::1"):
        print(
            "WARNING: bound to a non-loopback address; the dashboard has no authentication.",
            flush=True,
        )
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()
