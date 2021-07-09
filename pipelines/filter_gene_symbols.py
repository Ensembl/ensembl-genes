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
Filter the gene symbols assigned by the classifier to create the subset of them
to be loaded to the genome assembly core database.
"""


# standard library imports
import argparse
import pathlib

# third party imports
import pandas as pd

# project imports


def filter_gene_symbols(gene_symbols_tsv, threshold):
    """
    """
    gene_symbols_tsv_path = pathlib.Path(gene_symbols_tsv)

    all_symbols = pd.read_csv(gene_symbols_tsv_path, sep="\t")

    filtered_symbols = all_symbols.loc[all_symbols["probability"] >= threshold]

    filtered_symbols_tsv_path = pathlib.Path(
        f"{gene_symbols_tsv_path.parent}/{gene_symbols_tsv_path.stem}_filtered.csv"
    )
    filtered_symbols.to_csv(filtered_symbols_tsv_path, sep="\t", index=False)


def main():
    """
    main function
    """
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument(
        "--gene_symbols_tsv",
        help="gene symbols assignments TSV file path",
    )
    argument_parser.add_argument(
        "--threshold",
        default=0.9,
        type=float,
        help="gene symbol assignment probability threshold for including to the core db",
    )

    args = argument_parser.parse_args()

    if args.gene_symbols_tsv:
        filter_gene_symbols(args.gene_symbols_tsv, args.threshold)
    else:
        print("Error: missing argument.")
        print(__doc__)
        argument_parser.print_help()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted with CTRL-C, exiting...")
