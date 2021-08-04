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
    star2introns.py --num_cpus 12 --species salmo_salar
        --intron_db mysql+pymysql://rw_user:password@genebuild6:3306/salmo_salar_gca905237065v2_introns_104
        --csv_file /hps/nobackup/genebuild/salmon/rnaseq/salmo_salar.csv
        --junctions_dir /hps/nobackup/genebuild/salmon/rnaseq/output --verbose

"""

from typing import Dict
import sys
import logging
import argparse
from multiprocessing import Pool
from multiprocessing.pool import AsyncResult
from re import search
from pathlib import Path
from math import ceil
import csv
import sqlalchemy as db


def get_analyses(csv_file: str, species: str) -> Dict[str, str]:
    """Parse the csv_file to generate the logic names.
    If there is only one sample, it does not create
    the merged analysis"""
    with open(csv_file, newline="") as csvfile:
        samplereader = csv.reader(csvfile, delimiter="\t")
        unique_analyses = {}
        analyses = {}
        for row in samplereader:
            logic_name = "_".join([species, row[0], "rnaseq_daf"])
            unique_analyses[logic_name] = logic_name
            analyses[row[1]] = logic_name
        if len(unique_analyses.values()) > 1:
            analyses["merged"] = "_".join([species, "merged", "rnaseq_daf"])
    return analyses


def fetch_slice_ids(intron_db: str) -> Dict[str, int]:
    """ Retrieve the seq_region names and dbIDs """
    slice_ids = {}
    engine = db.create_engine(intron_db)
    with engine.connect() as connection:
        rows = connection.execute(
            db.text(
                "SELECT sr.seq_region_id, sr.name"
                + " FROM seq_region sr, seq_region_attrib sra"
                + " WHERE sr.seq_region_id = sra.seq_region_id AND sra.attrib_type_id = 6"
            )
        )
        for row in rows:
            slice_ids[row[1]] = row[0]
    engine.dispose()
    return slice_ids


def process_file(
        filename: str, analyses: Dict[str, str]
        ) -> Dict[str, Dict[str, Dict[str, int]]]:
    """Parse a junctions file and store the information in
    a dictionary structure: seq_region->position->logic_name"""
    file_id = search("^([^_]+)", filename.name)
    daf_table = {}
    if file_id:
        logic_name = analyses[file_id.group(1)]
        if "merged" in analyses:
            has_merged_analysis = analyses["merged"]
        with open(filename, newline="") as csvfile:
            intronreader = csv.reader(csvfile, delimiter="\t")
            for row in intronreader:
                seq_region = row[0]
                intron_id = ":".join(row[1:5])
                depth = int(row[6]) + ceil(int(row[7]) / 2)
                if seq_region in daf_table and intron_id in daf_table[seq_region]:
                    daf_table[seq_region][intron_id][logic_name] += depth
                    if has_merged_analysis:
                        daf_table[seq_region][intron_id][has_merged_analysis] += depth
                else:
                    if seq_region not in daf_table:
                        daf_table[seq_region] = {}
                    if has_merged_analysis:
                        daf_table[seq_region][intron_id] = {
                            logic_name: depth,
                            has_merged_analysis: depth,
                        }
                    else:
                        daf_table[seq_region][intron_id] = {logic_name: depth}
    return daf_table


def write_output(
        intron_db: str,
        analyses: Dict[str, str],
        slices: Dict[str, int],
        daf_table: Dict[str, Dict[str, Dict[str, int]]],
        batch_size: int,
    ) -> None:
    """ Stores the analyses and the junction information into an Ensembl database """
    engine = db.create_engine(intron_db)
    metadata = db.MetaData()

    analysis_table = db.Table("analysis", metadata, autoload=True, autoload_with=engine)
    dna_align_feature_table = db.Table(
        "dna_align_feature", metadata, autoload=True, autoload_with=engine
    )
    with engine.connect() as connection:
        logging.debug("Inserting analyses")
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

        logging.debug("Preparing daf stuff")
        daf_insert = dna_align_feature_table.insert()
        counter = 1
        daf_values = []
        logging.debug("Disable KEYS")
        connection.execute(db.text("ALTER TABLE dna_align_feature DISABLE KEYS"))
        logging.debug("Loading daf stuff")
        for seq_region in daf_table.keys():
            for intron_id in daf_table[seq_region].keys():
                seq_region_data = intron_id.split(":")
                for logic_name in daf_table[seq_region][intron_id].keys():
                    if seq_region_data[2] == "-1":
                        seq_region_data[2] = -1
                    else:
                        seq_region_data[2] = 1
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
                                int(seq_region_data[2]) - int(seq_region_data[1]) + 1
                            ),
                            "hit_strand": 1,
                            "align_type": "ensembl",
                            "analysis_id": analyses_id[logic_name],
                            "score": daf_table[seq_region][intron_id][logic_name],
                            "cigar_line": "{:d}M".format(
                                (int(seq_region_data[2]) - int(seq_region_data[1]) + 1)
                            ),
                        }
                    )
                    if (counter % batch_size) == 0:
                        connection.execute(daf_insert, daf_values)
                        daf_values = []
                        logging.debug("Loading daf stuff %d", counter)
                    counter += 1
        if daf_values:
            connection.execute(daf_insert, daf_values)
            daf_values = []
        logging.debug("Stored %d daf stuff", counter)
        connection.execute(db.text("ALTER TABLE dna_align_feature ENABLE KEYS"))
        logging.debug("KEYS enabled")

    engine.dispose()


def work_done(results: AsyncResult) -> None:
    """Print 'Successful X objects' when all files have been parsed successfully and
    --verbose was specified on the command line"""
    logging.debug("Successful %d objects", len(results))


def work_failed(error: BaseException) -> None:
    """ Print the error reported by one of the process and exits """
    logging.error(error)
    sys.exit()


def main() -> None:
    """It will prepare the analyses, then fetch the seq_region ids.
    Then it will create a pool of worker to prcoess as many files as
    there is cores assigned.
    Finally it will collapse the results and write the results into a database"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--junctions_dir",
        type=str,
        help="directory with the files describing junctions found by STAR and me",
    )
    parser.add_argument(
        "--csv_file",
        type=str,
        help="Tabulated file containing the sample names and the sample ids",
    )
    parser.add_argument(
        "--species",
        type=str,
        help="Species name, used to construct the logic_name of the analysis",
    )
    parser.add_argument(
        "--num_cpus", type=int, help="Number of cores to use", default=1
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

    log_format = "%(asctime)s: %(message)s"
    log_level = logging.INFO
    if args.verbose:
        log_level = logging.DEBUG
    logging.basicConfig(format=log_format, level=log_level, datefmt="%H:%M:%S")

    analyses = get_analyses(args.csv_file, args.species)
    files = []
    for filename in Path(args.junctions_dir).glob("*SJ.out.tab"):
        files.append([filename, analyses])

    if len(files) > 0:
        slices = fetch_slice_ids(args.intron_db)
    else:
        logging.error("Could not find any files")
        sys.exit()

    with Pool(processes=args.num_cpus) as executor:
        logging.debug("Starting")
        results = executor.starmap_async(process_file, files, 1, work_done, work_failed)
        executor.close()
        executor.join()
        daf_table = {}
        index = 1
        logging.debug("Process results")
        for item in results.get():
            index += 1
            for seq_region in item.keys():
                if seq_region not in daf_table:
                    daf_table[seq_region] = {}
                for intron_id in item[seq_region].keys():
                    for logic_name in item[seq_region][intron_id].keys():
                        if (
                            intron_id in daf_table[seq_region]
                            and logic_name in daf_table[seq_region][intron_id]
                        ):
                            daf_table[seq_region][intron_id][logic_name] += item[
                                seq_region
                            ][intron_id][logic_name]
                        else:
                            if intron_id in daf_table[seq_region]:
                                daf_table[seq_region][intron_id][logic_name] = item[
                                    seq_region
                                ][intron_id][logic_name]
                            else:
                                daf_table[seq_region][intron_id] = {
                                    logic_name: item[seq_region][intron_id][logic_name]
                                }

    logging.debug("Write to database")
    write_output(args.intron_db, analyses, slices, daf_table, args.batch_size)


if __name__ == "__main__":
    main()
