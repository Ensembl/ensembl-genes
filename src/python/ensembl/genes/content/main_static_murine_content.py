#!/usr/bin/env python3
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
# pylint:disable=duplicate-code
"""Prepare static content files for main release.
This script fetches assembly information from a MySQL database and the ENA API,
and generates HTML files containing assembly and annotation content.
"""

import argparse
from pathlib import Path
import pymysql
import requests
import xmltodict


def mysql_fetch_data(query: str, database, host, port, user):
    """
    Fetch data from a MySQL database based on a provided query.

    Args:
        query (str): SQL query to execute.
        database (str): Database name to connect to.
        host (str): Host address of the MySQL server.
        port (int): Port number of the MySQL server.
        user (str): Username for the MySQL server.

    Returns:
        tuple: A tuple containing the fetched data.
    """
    try:
        conn = pymysql.connect(
            host=host, user=user, port=port, database=database.strip()
        )
        cursor = conn.cursor()
        cursor.execute(query)
        info = cursor.fetchall()
    except pymysql.Error as err:
        print(f"Database error: {err}")
        info = None
    finally:
        cursor.close()
        conn.close()
    return info


def get_assembly_info(accession):
    """
    Retrieve assembly information from the ENA API for a given accession.

    Args:
        accession (str): The assembly accession number to fetch information for.

    Returns:
        dict: A dictionary containing assembly attributes such as name, level, \
        submitter, and counts.
    """
    assembly_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    assembly_xml = requests.get(assembly_url, timeout=10)
    assembly_dict = xmltodict.parse(assembly_xml.text)
    assembly_attribs = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["ASSEMBLY_ATTRIBUTES"][
        "ASSEMBLY_ATTRIBUTE"
    ]

    return_dict = {}

    # Retrieve specific attributes with error handling
    return_dict["assembly.name"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"].get(
        "NAME", ""
    )
    return_dict["assembly.level"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"].get(
        "ASSEMBLY_LEVEL", ""
    )
    return_dict["assembly.submitter"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"][
        "IDENTIFIERS"
    ]["SUBMITTER_ID"].get("@namespace", "")

    # Fetch specific attributes based on tags
    for attrib_set in assembly_attribs:
        if attrib_set["TAG"] == "count-contig":
            return_dict["contig.count"] = attrib_set["VALUE"]
        elif attrib_set["TAG"] == "scaffold-count":
            return_dict["scaffold.count"] = attrib_set["VALUE"]
        elif attrib_set["TAG"] == "ENA-LAST-UPDATED":
            return_dict["assembly.date"] = attrib_set["VALUE"]
        elif attrib_set["TAG"] == "n50":
            return_dict["scaffold.n50"] = attrib_set["VALUE"]
        elif attrib_set["TAG"] == "contig-n50":
            return_dict["contig.n50"] = attrib_set["VALUE"]

    return return_dict


def write_content(info, out_dir, url_path, species_strain):
    """
    Write assembly and annotation content to HTML files.

    Args:
        info (dict): Assembly information to be written.
        out_dir (Path): Directory path for output files.
        url_path (str): Name used to title the output files.
        species_strain (str): Species strain name for the assembly.
    """
    with open(  # pylint: disable=unspecified-encoding
        out_dir / f"{url_path}_assembly.html", "w"
    ) as assembly_out:
        print(
            f"<p>The assembly for {species_strain} was generated as part of "
            f'<a href="https://www.mousegenomes.org/">The Mouse Genomes Project'
            f"</a>, additional species_strain can be found in "
            f'<a href="https://www.ensembl.org/Mus_musculus/Info/Strains">Ensembl</a>.</p>'
            f"<p>The assembly is on the {info['assembly.level']} "
            f"level, consisting of {info['contig.count']} contigs assembled "
            f"into {info['scaffold.count']} scaffolds. "
            f"The N50 size is the length such that 50% of the assembled genome lies in "
            f"blocks of the N50 size or longer. "
            f"The N50 length for the contigs is {info['contig.n50']} while the "
            f"scaffold N50 is {info['scaffold.n50']}.</p>",
            file=assembly_out,
        )

    with open(  # pylint: disable=unspecified-encoding
        out_dir / f"{url_path}_annotation.html", "w"
    ) as annotation_out:
        print(
            f"<p>Genome annotation was generated by mapping "  # pylint: disable=f-string-without-interpolation
            f'<a href="https://www.gencodegenes.org/mouse/release_M30.html">GENCODE M30</a> '
            f'genes and transcripts via the <a href="https://beta.ensembl.org/help/articles/human-genome-automated-annotation">Ensembl Human automated annotation system</a>, supplemented by methods from the '  # pylint: disable=line-too-long
            f'<a href="https://beta.ensembl.org/help/articles/vertebrate-genome-annotation">Ensembl vertebrate annotation pipeline</a>. '  # pylint: disable=line-too-long
            f"Mapped GENCODE structures served as the primary evidence with "
            f"gaps in the annotations filled using aligned short-read "
            f"transcriptomic data and full-length transcripts derived from PacBio IsoSeq "
            f"long-read data."
            f"</p><p>In accordance with the "
            f'<a href="https://en.wikipedia.org/wiki/Fort_Lauderdale_Agreement">Fort Lauderdale '
            f"Agreement</a>, please check the publication "
            f"status of the genome/assembly before publishing any genome-wide analyses using "
            f"these data.</p>",
            file=annotation_out,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare static content files for main release"
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Output directory for files. Uses current dir by default.",
    )
    parser.add_argument("-d", "--db_name", help="Database name", required=True)
    parser.add_argument("-s", "--host", help="Host server", required=True)
    parser.add_argument(
        "-p", "--port", type=int, help="Host server port", required=True
    )

    args = parser.parse_args()

    server_info = {
        "server": {
            "db_host": args.host,
            "db_port": args.port,
            "db_user": "ensro",
        },
    }

    db = args.db_name
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Fetch assembly accession and production name from the core database
    CORE_QUERY = (
        "SELECT meta_key, meta_value FROM meta WHERE meta_key IN "
        "('assembly.accession', 'species.url', 'species.strain');"
    )
    core_meta = mysql_fetch_data(
        CORE_QUERY,
        host=server_info["server"]["db_host"],
        user=server_info["server"]["db_user"],
        port=server_info["server"]["db_port"],
        database=db,
    )

    # Get assembly info and write content files
    if core_meta is None:
        raise RuntimeError("Failed to fetch metadata from the database.")
    core_dict = {meta_pair[0]: meta_pair[1] for meta_pair in core_meta}

    assembly_info = get_assembly_info(core_dict["assembly.accession"])
    url_name = core_dict["species.url"]
    strain = core_dict["species.strain"]

    write_content(assembly_info, output_dir, url_name, strain)
