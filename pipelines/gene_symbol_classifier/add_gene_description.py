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
Add gene description to the assignments and save to a new CSV.
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


def add_gene_description(assignments_csv):
    assignments_path = pathlib.Path(assignments_csv)

    assignments = pd.read_csv(assignments_path, sep="\t")

    # generate gene description values
    assignments["gene_description"] = [
        "{} [Ensembl NN prediction with score {}%]".format(
            row[4],
            100 * round(row[3], 4),
        )
        for row in assignments.itertuples()
    ]

    # reorder dataframe columns
    columns = [
        "stable_id",
        "symbol",
        "gene_description",
        "probability",
        "description",
        "source",
    ]
    assignments = assignments[columns]

    assignments_with_description_path = pathlib.Path(
        f"{assignments_path.parent}/{assignments_path.stem}_description.csv"
    )
    assignments.to_csv(assignments_with_description_path, sep="\t", index=False)
    logger.info(f"assignments CSV with gene description saved at {assignments_with_description_path}")


def main():
    """
    main function
    """
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument(
        "--assignments",
        help="gene symbol assignments CSV file path",
    )

    args = argument_parser.parse_args()

    # set up logger
    logger.remove()
    logger.add(sys.stderr, format=logging_format)

    if args.assignments:
        add_gene_description(args.assignments)
    else:
        print("Error: missing argument.")
        print(__doc__)
        argument_parser.print_help()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted with CTRL-C, exiting...")
