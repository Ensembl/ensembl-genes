"""``annotation-qc dashboard --config EXPERIMENT.json [--config …] [--port 8765]``"""

from __future__ import annotations

import sys


def cmd_dashboard(args):
    from ensembl.genes.annotation_qc.benchmark.config import (
        ConfigError,
        load_experiment,
    )
    from ensembl.genes.annotation_qc.benchmark.workspace import Workspace
    from ensembl.genes.annotation_qc.dashboard.server import serve

    dbs = {}
    for path in args.config:
        try:
            experiment = load_experiment(path)
        except ConfigError as error:
            sys.exit(str(error))
        if experiment.id in dbs:
            sys.exit(f"experiment id {experiment.id!r} is configured twice")
        dbs[experiment.id] = Workspace(experiment).dashboard_db()
    serve(dbs, args.host, args.port, args.verbose)


def register(subparsers):
    parser = subparsers.add_parser(
        "dashboard", help="Serve the local annotation-QC dashboard (localhost)."
    )
    parser.add_argument(
        "--config",
        action="append",
        required=True,
        help="Experiment configuration; repeat for several.",
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Bind address (default 127.0.0.1, local only).",
    )
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--verbose", action="store_true", help="Log every request.")
    parser.set_defaults(func=cmd_dashboard)
