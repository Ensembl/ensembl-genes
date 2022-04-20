# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" Module to run Red to find repeats and store them in the given Ensembl core database """

import os
import errno
import subprocess
import tempfile
import shutil
import sqlalchemy as db
from pathlib import Path

# sqlalchemy requires MySQLdb but MySQLdb doesn't support Python 3.x
# pymysql can be imported and used instead
import pymysql

pymysql.install_as_MySQLdb()

import eHive


class Repeatmask_Red(eHive.BaseRunnable):
    """Runnable that runs Red to find repeats and store them in the target database."""

    def param_defaults(self):
        """It sets the parameters default values."""
        return {
            "logic_name": "repeatdetector",
            "target_db_url": "",  # 'driver://user:pass@host:port/dbname'
            "red_path": "",
            "red_meta_key": 0,
        }

    def fetch_input(self):
        """It fetches the input parameters and it checks that they are correct."""

        # get new temporary directory name in the default path with prefix gnm_ (eg '/scratch/gnm_nar4nry8')
        temp_gnm_dir = tempfile.TemporaryDirectory(prefix="gnm_").name

        genome_file = self.param_required("genome_file")
        gnm = self.param("gnm", temp_gnm_dir)
        msk = self.param_required("msk")
        rpt = self.param_required("rpt")
        red_path = self.param_required("red_path")
        target_db_url = self.param_required("target_db_url")
        self.param("target_db_url", f"{target_db_url}?local_infile=1")

        # make the temporary directory 'gnm' and copy the genome file into it
        # in this way we make sure that the only .fa file to be processed is the one we want
        gnm_path = Path(gnm)
        try:
            gnm_path.mkdir()
        except PermissionError:
            print(f"Could not create {gnm_path} directory.")
            raise

        # copy the genome file into the temporary directory
        # add suffix '.fa' to make sure it ends with '.fa' as required by Red
        # In this way we make sure that the only .fa file to be processed is the one we want
        genome_file_path = Path(genome_file)
        new_genome_file = gnm_path / f"{genome_file_path.name}.fa"
        try:
            os.symlink(genome_file, new_genome_file)
        except PermissionError:
            print(f"Could not create symlink to {genome_file} in directory {gnm}")
            raise

        genome_file = self.param(genome_file, new_genome_file)

        # check that the Red binary exists
        red_path_obj = Path(red_path)
        if not (red_path_obj.exists()):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), red_path)

        # connect to the target database and fetch the seq_region_ids required later
        engine = db.create_engine(self.param("target_db_url"))
        connection = engine.connect()
        if connection:
            metadata = db.MetaData()

            sr = db.Table("seq_region", metadata, autoload=True, autoload_with=engine)
            sra = db.Table(
                "seq_region_attrib", metadata, autoload=True, autoload_with=engine
            )
            at = db.Table("attrib_type", metadata, autoload=True, autoload_with=engine)

            query = db.select([sr.columns.seq_region_id, sr.columns.name])
            query = query.select_from(
                sr.join(
                    sra, sr.columns.seq_region_id == sra.columns.seq_region_id
                ).join(at, at.columns.attrib_type_id == sra.columns.attrib_type_id)
            )
            query = query.where(at.columns.code == "toplevel")

            results = connection.execute(query).fetchall()

            seq_region = {}
            for seq_region_id, name in results:
                seq_region[name] = seq_region_id
            self.param("seq_region", seq_region)

            connection.close()
        else:
            raise ValueError(
                f"Could not connect to the target database {target_db_url}."
            )

        engine.dispose()

        # make sure that the output directories exist and they are empty
        msk_path = Path(msk)
        if msk_path.exists():
            try:
                shutil.rmtree(msk)
            except OSError as e:
                print(f"Error: {msk} : {e.strerror}")

        rpt_path = Path(rpt)
        if rpt_path.exists():
            try:
                shutil.rmtree(rpt)
            except OSError as e:
                print(f"Error: {rpt} : {e.strerror}")

        try:
            msk_path.mkdir()
        except PermissionError:
            print(f"Could not create {msk} directory.")
            raise

        try:
            rpt_path.mkdir()
        except PermissionError:
            print(f"Could not create {rpt} directory.")
            raise

    def run(self):
        """It runs the Red program."""
        # output format: 1 (chrName:start-end) or 2 (chrName start end)
        # Note that chrName includes the '>' character
        cmd = [
            self.param("red_path"),
            "-frm",
            "2",
            "-gnm",
            self.param("gnm"),
            "-msk",
            self.param("msk"),
            "-rpt",
            self.param("rpt"),
        ]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as err:
            print(
                f'Could not run Red. Command: {" ".join(cmd)} Return code {str(err.returncode)}'
            )

    def write_output(self):
        """It parses the Red's program output and inserts it into
        the given Ensembl core database."""
        engine = db.create_engine(self.param("target_db_url"))
        connection = engine.connect()
        metadata = db.MetaData()

        analysis_table = db.Table(
            "analysis", metadata, autoload=True, autoload_with=engine
        )
        meta_table = db.Table("meta", metadata, autoload=True, autoload_with=engine)
        repeat_consensus_table = db.Table(
            "repeat_consensus", metadata, autoload=True, autoload_with=engine
        )
        db.Table("repeat_feature", metadata, autoload=True, autoload_with=engine)

        # insert Red analysis
        analysis_insert = (
            analysis_table.insert(None)
            .prefix_with("IGNORE")
            .values(
                {
                    "created": db.sql.func.now(),
                    "logic_name": self.param("logic_name"),
                    "program": self.param("logic_name"),
                    "program_version": "05/22/2015",
                    "program_file": self.param("red_path"),
                }
            )
        )
        connection.execute(analysis_insert)

        # fetch the inserted analysis_id
        analysis_query = db.select([analysis_table.columns.analysis_id])
        analysis_query = analysis_query.where(
            analysis_table.columns.logic_name == self.param("logic_name")
        )
        analysis_results = connection.execute(analysis_query).fetchall()
        analysis_id = analysis_results[0][0]

        # insert repeat analysis meta keys
        if self.param("red_meta_key"):
            meta_insert = (
                meta_table.insert(None)
                .prefix_with("IGNORE")
                .values(
                    {
                        "species_id": 1,
                        "meta_key": "repeat.analysis",
                        "meta_value": self.param("logic_name"),
                    }
                )
            )
            connection.execute(meta_insert)

        # insert dummy repeat consensus
        repeat_consensus_insert = repeat_consensus_table.insert(None).values(
            {
                "repeat_name": self.param("logic_name"),
                "repeat_class": self.param("logic_name"),
                "repeat_type": self.param("logic_name"),
                "repeat_consensus": "N",
            }
        )
        connection.execute(repeat_consensus_insert)

        # fetch the inserted repeat_consensus_id
        repeat_consensus_query = db.select(
            [repeat_consensus_table.columns.repeat_consensus_id]
        )
        repeat_consensus_query = repeat_consensus_query.where(
            repeat_consensus_table.columns.repeat_name == self.param("logic_name")
        )
        repeat_consensus_results = connection.execute(repeat_consensus_query).fetchall()
        repeat_consensus_id = repeat_consensus_results[0][0]

        # parse the repeats file and make a tsv file ready to load into the repeat_feature table
        repeats_file = self.parse_repeats(
            self.param("rpt"), repeat_consensus_id, analysis_id
        )

        # insert repeat features
        repeat_feature_query = (
            f'LOAD DATA LOCAL INFILE "{repeats_file}"'
            + "INTO TABLE repeat_feature \
                                 FIELDS TERMINATED BY '\\t' \
                                 LINES TERMINATED BY '\\n' \
                                 (seq_region_id,seq_region_start,seq_region_end, \
                                  repeat_start,repeat_end,repeat_consensus_id,analysis_id)"
        )
        connection.execute(repeat_feature_query)

        connection.close()
        engine.dispose()

        # delete temporary directory and its contents
        try:
            shutil.rmtree(self.param("gnm"))
        except OSError as e:
            print(f'Error: {self.param("gnm")} : {e.strerror}')

    def parse_repeats(self, rpt, repeat_consensus_id, analysis_id):
        """It parses the Red's program output and it converts it into
        a tsv file which can be loaded into an Ensembl core repeat_feature table."""
        # Required 1 file in rpt dir and it ends with .rpt
        # Red's rpt output file contains ">" which needs to be removed from each line
        # and we need to replace the seq region name with seq region id and
        # add some extra columns so it can be loaded directly
        seq_region = self.param("seq_region")
        rpt_path = Path(rpt)
        rpt_files = list(rpt_path.iterdir())
        rpt_file = Path(rpt_files[0])  # we know there is only one file
        fixed_rpt_file = Path(f"{rpt_file}.fixed")
        with rpt_file.open("r") as f_in, fixed_rpt_file.open("w") as f_out:
            for line in f_in:
                columns = line.split()

                if columns[0][0] == ">":
                    name = columns[0][1:]  # remove first character '>'
                else:
                    name = columns[0]

                seq_region_start = int(columns[1]) + 1  # Red's start is zero-based
                seq_region_end = int(columns[2]) - 1  # Red's end is exclusive
                seq_region_id = seq_region[name]
                # seq_region_id seq_region_start seq_region_end repeat_start repeat_end repeat_consensus_id analysis_id
                print(
                    str(seq_region_id)
                    + "\t"
                    + str(seq_region_start)
                    + "\t"
                    + str(seq_region_end)
                    + "\t"
                    + "1\t"
                    + str(seq_region_end - seq_region_start + 1)
                    + "\t"
                    + str(repeat_consensus_id)
                    + "\t"
                    + str(analysis_id),
                    file=f_out,
                )

        return fixed_rpt_file
