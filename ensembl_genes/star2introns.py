#!/usr/bin/env python3

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" This script will process all junctions files created by STAR
    and store the information in an Ensembl core database
    using a tabulated file to group the data by sample names

    It expects the first column of the tabulated file to be the sample
    name and the second column of the file to be the accession/id of the
    sample. The second column value should match the following regex: (^[^_]+)
    on the junction file

Examples:
    star2introns --tsv_file /hps/nobackup/genebuild/salmon/rnaseq/salmo_salar.tsv
        --intron_db mysql+pymysql://rw_user:password@genebuild6:3306/salmo_salar_gca905237065v2_104
        --junctions_dir /hps/nobackup/genebuild/salmon/rnaseq/output
        --species salmo_salar --verbose

"""

from typing import Dict, Optional
import sys
import logging
import argparse
from re import search
from pathlib import Path
from math import ceil
import csv
import sqlalchemy as db
from sqlalchemy.engine import Engine  # Needed for typing


def get_engine(url: str) -> Engine:
    """Create SQLAlchemy engine.

    Args:
        url: URI to the database.

    Returns:
        SQLAlchemy engine.
    """
    return db.create_engine(url)


def get_analyses(csv_file: str, species: str) -> Dict[str, str]:
    """Parse the csv_file to generate the logic names.

    Args:
        csv_file: Tabulated separated file where row[0] is the sample id and row[1] is the file id.
        species: scientific name of the species, joined by underscores.

    Returns:
        A dictionary where key is the file id and value is the logic_name.
    """

    with open(csv_file, newline="", encoding="utf-8") as csvfile:
        samplereader = csv.reader(csvfile, delimiter="\t")
        analyses = {}
        for row in samplereader:
            logic_name = "_".join([species, row[0], "rnaseq_daf"])
            analyses[row[1]] = logic_name
    return analyses


def fetch_slice_ids(engine: Engine) -> Dict[str, int]:
    """Retrieve the seq_region names and dbIDs.

    Args:
        engine: SQLAlchemy Engine object to connect to the database.

    Returns
        A dictionary where key is the sequence name and value is its dbID.
    """

    slice_ids = {}
    with engine.connect() as connection:
        rows = connection.execute(
            db.text(
                "SELECT sr.seq_region_id, sr.name"
                + " FROM seq_region sr, seq_region_attrib sra, attrib_type at"
                + " WHERE sr.seq_region_id = sra.seq_region_id"
                + " AND sra.attrib_type_id = at.attrib_type_id"
                + " AND at.code = 'toplevel'"
            )
        )
        for row in rows:
            slice_ids[row[1]] = row[0]
    return slice_ids


def process_file(
    filename: Path,
    analyses: Dict[str, str],
    daf_table: Optional[Dict[str, Dict[str, Dict[str, int]]]] = None,
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """Parses a STAR junctions file

    It stores the information in a complex dictionary structure:
    seq_region->position->logic_name

    Args:
        filename: Path of file to process.
        analyses: Dictionary of logic_name to assign depending on the filename.
        daf_table: Dictionary for storing the depth of the intron depending on the sample.

    Returns:
        A structure of dictionaries storing the intron information:
            dict[seq_region,
                dict[intron_id,
                    dict[analysis, depth]
                    ]
                ]
    """

    logger = logging.getLogger("star2introns")
    logger.debug("Processing %s", filename.name)
    file_id = search("^([^_]+)", filename.name)
    if daf_table is None:
        daf_table = {}
    if file_id:
        try:
            logic_name = analyses[file_id.group(1)]
            with open(filename, newline="", encoding="utf-8") as csvfile:
                intronreader = csv.reader(csvfile, delimiter="\t")
                for row in intronreader:
                    seq_region = row[0]
                    if row[3] == "2":
                        intron_strand = -1
                    else:
                        intron_strand = 1
                    intron_id = f"{row[1]}:{row[2]}:{intron_strand}:{row[4]}"
                    depth = int(row[6]) + ceil(int(row[7]) / 2)
                    if seq_region not in daf_table:
                        daf_table[seq_region] = {}
                    if intron_id not in daf_table[seq_region]:
                        daf_table[seq_region][intron_id] = {}
                    if logic_name in daf_table[seq_region][intron_id]:
                        daf_table[seq_region][intron_id][logic_name] += depth
                    else:
                        daf_table[seq_region][intron_id][logic_name] = depth
        except KeyError:
            logger.error("Could not find analysis for file %s", filename)

    return daf_table


def write_output(  # pylint: disable=too-many-locals
    engine: Engine,
    analyses: Dict[str, str],
    slices: Dict[str, int],
    daf_table: Dict[str, Dict[str, Dict[str, int]]],
    batch_size: int,
) -> None:
    """Stores the analyses and the junction information into an Ensembl database.

    Before inserting the data into the dna_align_feature table we disable the indexes
    to speed up the load. The indexes are enabled after all the rows have been inserted.

    Args:
        engine: SQLAlchemy Engine object to connect to the database.
        analyses: Dictionary where key is the file id and value is the logic_name.
        slices: Dictionary where key is the sequence name and value is the dbID.
        daf_table: Dictionary with the data to store in the dna_align_feature table.
        batch_size: The number of rows of data for each insert.
    """

    metadata = db.MetaData()

    analysis_table = db.Table("analysis", metadata, autoload=True, autoload_with=engine)
    dna_align_feature_table = db.Table(
        "dna_align_feature", metadata, autoload=True, autoload_with=engine
    )
    logger = logging.getLogger("star2introns")
    with engine.connect() as connection:
        logger.info("Inserting analyses")
        analyses_id = {}
        analysis_insert = (
            analysis_table.insert(bind=db.bindparam("logic_name"))
            .prefix_with("IGNORE")
            .values({"created": db.sql.func.now()})
        )
        analysis_query = db.select([analysis_table.columns.analysis_id]).where(
            analysis_table.columns.logic_name == db.bindparam("logic_name")
        )
        for analysis in analyses.keys():
            connection.execute(analysis_insert, {"logic_name": analyses[analysis]})
            # fetch the inserted analysis_id
            analysis_results = connection.execute(
                analysis_query, {"logic_name": analyses[analysis]}
            ).fetchall()
            analyses_id[analyses[analysis]] = analysis_results[0][0]

        daf_insert = dna_align_feature_table.insert()
        counter = 1
        daf_values = []
        logger.debug("Disable KEYS")
        connection.execute(db.text("ALTER TABLE dna_align_feature DISABLE KEYS"))
        logger.info("Loading daf stuff")
        for seq_region in daf_table.keys():
            for intron_id in daf_table[seq_region].keys():
                seq_region_data = intron_id.split(":")
                for logic_name in daf_table[seq_region][intron_id].keys():
                    if int(seq_region_data[3]) > 0:
                        hit_name = f"{counter}:canon"
                    else:
                        hit_name = f"{counter}:non canon"
                    daf_values.append(
                        {
                            "seq_region_id": slices[seq_region],
                            "seq_region_start": seq_region_data[0],
                            "seq_region_end": seq_region_data[1],
                            "seq_region_strand": seq_region_data[2],
                            "hit_name": hit_name,
                            "hit_start": 1,
                            "hit_end": (
                                int(seq_region_data[1]) - int(seq_region_data[0]) + 1
                            ),
                            "hit_strand": 1,
                            "align_type": "ensembl",
                            "analysis_id": analyses_id[logic_name],
                            "score": daf_table[seq_region][intron_id][logic_name],
                            "cigar_line": f"{int(seq_region_data[1])-int(seq_region_data[0])+1}M",
                        }
                    )
                    if (counter % batch_size) == 0:
                        connection.execute(daf_insert, daf_values)
                        daf_values = []
                        logger.debug("Loading daf stuff %d", counter)
                    counter += 1
        if daf_values:
            connection.execute(daf_insert, daf_values)
            daf_values = []
        logger.info("Stored %d daf stuff", counter)
        connection.execute(db.text("ALTER TABLE dna_align_feature ENABLE KEYS"))
        logger.debug("KEYS enabled")


def main() -> None:
    """Retrieve splice site information for database storage.

    It will retrieve the possible analyses based on the TSV file provided. Then it will
    process all the STAR junctions files in the directory provided and store the information
    in the database provided"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--junctions_dir",
        type=str,
        help="directory with the files describing junctions found by STAR and me",
    )
    parser.add_argument(
        "--tsv_file",
        type=str,
        help="Tabulated file containing the sample names and the sample ids",
    )
    parser.add_argument(
        "--species",
        type=str,
        help="Species name, used to construct the logic_name of the analysis",
    )
    parser.add_argument(
        "--intron_db", type=str, help="URI to database using 'mysql+pymysql'"
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        help="Number of values to insert at each batch",
        default=500,
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Print debugging message"
    )
    args = parser.parse_args()

    logger = logging.getLogger("star2introns")
    logger.setLevel(logging.INFO)
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    analyses = get_analyses(args.tsv_file, args.species)

    daf_table = {}
    for filename in Path(args.junctions_dir).glob("*SJ.out.tab"):
        process_file(filename, analyses, daf_table)

    if len(daf_table) > 0:
        engine = get_engine(args.intron_db)
        slices = fetch_slice_ids(engine)
        write_output(engine, analyses, slices, daf_table, args.batch_size)
        # I'm not sure the dispose() is needed
        engine.dispose()
    else:
        logger.error("Could not load any data")
        sys.exit()


# Setting up the logger
main_logger = logging.getLogger("star2introns")
console_handler = logging.StreamHandler()
console_format = logging.Formatter(
    "%(asctime)s| %(levelname)s | %(module)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
console_handler.setFormatter(console_format)
main_logger.addHandler(console_handler)

if __name__ == "__main__":
    main()
