#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
Filter the gene symbol assignments from the GSC network to create the subset of them
to be loaded to the genome assembly core database.
"""


# standard library imports
import argparse
import pathlib
import sys

# third party imports
import pandas as pd

from loguru import logger

# project imports


logging_format = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{message}</level>"


def filter_gene_symbols(symbol_assignments, threshold):
    gene_symbols_csv_path = pathlib.Path(symbol_assignments)

    all_symbols = pd.read_csv(gene_symbols_csv_path, sep="\t")

    filtered_symbols = all_symbols.loc[all_symbols["probability"] >= threshold]

    filtered_symbols_csv_path = pathlib.Path(
        f"{gene_symbols_csv_path.parent}/{gene_symbols_csv_path.stem}_filtered.csv"
    )
    filtered_symbols.to_csv(filtered_symbols_csv_path, sep="\t", index=False)
    logger.info(f"filtered assignments saved at {filtered_symbols_csv_path}")


def main():
    """
    main function
    """
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument(
        "--symbol_assignments",
        help="gene symbol assignments CSV file path",
    )
    argument_parser.add_argument(
        "--threshold",
        default=0.9,
        type=float,
        help="gene symbol assignment probability threshold for including to the core db",
    )

    args = argument_parser.parse_args()

    # set up logger
    logger.remove()
    logger.add(sys.stderr, format=logging_format)

    if args.symbol_assignments:
        filter_gene_symbols(args.symbol_assignments, args.threshold)
    else:
        print("Error: missing argument.")
        print(__doc__)
        argument_parser.print_help()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted with CTRL-C, exiting...")
