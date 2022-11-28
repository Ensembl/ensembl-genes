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

""" This module will process all junctions files created by STAR
    and store the information in an Ensembl core database
    using a tabulated file to group the data by sample names

    It expects the first column of the tabulated file to be the sample
    name and the second column of the file to be the accession/id of the
    sample. The second column value should match the following regex: (^[^_]+)
    on the junction file

"""
import logging
import pathlib
import eHive
from . import star2introns


class StarIntrons(eHive.BaseRunnable):
    """Load the introns from the Star alignments."""

    def param_defaults(self):
        """It sets the parameters default values.

        Returns:
            A dictionary.
        """
        defaults = super().param_defaults()

        defaults["batch_size"] = 500
        return defaults

    def fetch_input(self):
        """Retrieve all informations needed for processing the junction files.

        It searches all files ending with SJ.out.tab in the 'junctions_dir' directory.
        It will parse the TSV file and combine it with the sientific name to generate
        the intron logic_name for each samples. It will finally fetch the seq_region ids
        to be able to store the introns in the database.

        Raises:
            ValueError: when no junction files are found.
        """

        if self.debug > 0:
            # The name should match the logger name in the script
            logger = logging.getLogger("star2introns")
            logger.setLevel(logging.DEBUG)

        junction_files = []
        for filename in pathlib.Path(self.param_required("junctions_dir")).glob(
            "*SJ.out.tab"
        ):
            junction_files.append(filename)

        if len(junction_files) > 0:
            target_url = self.param_required("intron_db").replace(
                "mysql:", "mysql+pymysql:"
            )
            engine = star2introns.get_engine(target_url)
            self.param(
                "analyses",
                star2introns.get_analyses(
                    self.param_required("tsv_file"),
                    star2introns.fetch_species_name(engine),
                ),
            )
            self.param("slices", star2introns.fetch_slice_ids(engine))
            self.param("junction_files", junction_files)
        else:
            raise ValueError(
                f"No junctions files found in {self.param_required('junctions_dir')}"
            )

    def run(self):
        """Process all the files found in 'junctions_dir'.

        It will collapse all the introns which are from the same sample, the score representing
        the number of reads overlapping the splice site.
        """

        daf_table = {}
        for filename in self.param("junction_files"):
            star2introns.process_file(filename, self.param("analyses"), daf_table)

        self.param("output", daf_table)

    def write_output(self):
        """Write the introns into the dna_align_feature table."""

        star2introns.write_output(
            self.param("engine"),
            self.param("analyses"),
            self.param("slices"),
            self.param("output"),
            self.param("batch_size"),
        )
