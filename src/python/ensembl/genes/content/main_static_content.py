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
"""
Module to create static content files for main release.
This module fetches assembly information from a MySQL database and the ENA API,
and generates HTML files containing assembly and annotation content.
"""
import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple
import pymysql
import requests
import xmltodict


def mysql_fetch_data(
    query: str, database: str, host: str, port: int, user: str
) -> Optional[Tuple]:
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


def get_assembly_info(accession: str) -> Dict[str, str]:
    """
    Retrieve assembly information from the ENA API for a given accession.

    Args:
        accession (str): The assembly accession number to fetch information for.

    Returns:
        dict: A dictionary containing assembly attributes such as name, level,\
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


def write_content(info: Dict[str, str], out_dir: Path, url_path: str) -> None:
    """
    Write assembly and annotation content to HTML files.

    Args:
        info (dict): Assembly information to be written.
        output_dir (Path): Directory path for output files.
        url_path (str): Name used to title the output files.
    """
    with open(  # pylint: disable=unspecified-encoding
        out_dir / f"{url_path}_assembly.html", "w"
    ) as assembly_out:  # pylint: disable=unspecified-encoding
        print(
            f"<p>The {info['assembly.name']} assembly was submitted by "
            f"{info['assembly.submitter']} and last updated on "
            f"{info['assembly.date']}. The assembly is on the "
            f"{info['assembly.level']} "
            f"level, consisting of {info['contig.count']} contigs assembled into "
            f"{info['scaffold.count']} scaffolds. "
            f"The N50 size is the length such that 50% of the assembled genome lies in "
            f"blocks of the N50 size or longer. "
            f"The N50 length for the contigs is {info['contig.n50']} while the "
            f"scaffold N50 is {info['scaffold.n50']}.</p>",
            file=assembly_out,
        )

    with open(  # pylint: disable=unspecified-encoding
        out_dir / f"{url_path}_annotation.html", "w"
    ) as annotation_out:  # pylint: disable=unspecified-encoding
        print(
            "<p>Genome annotation was generated using the "
            '<a href="https://beta.ensembl.org/help/articles/vertebrate-genome-annotation">Ensembl vertebrate annotation pipeline</a>. '  # pylint: disable=line-too-long
            "</p><p>In accordance with the "
            '<a href="https://en.wikipedia.org/wiki/Fort_Lauderdale_Agreement">Fort Lauderdale Agreement</a>, please check the publication '  # pylint: disable=line-too-long
            "status of the genome/assembly before publishing any genome-wide analyses using these data.</p>",  # pylint: disable=line-too-long
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
        "('assembly.accession', 'species.url');"
    )
    core_meta = mysql_fetch_data(
        CORE_QUERY,
        host=server_info["server"]["db_host"],
        user=server_info["server"]["db_user"],
        port=server_info["server"]["db_port"],
        database=db,
    )

    if core_meta is None:
        raise RuntimeError("Failed to fetch metadata from the database.")
    core_dict = {meta_pair[0]: meta_pair[1] for meta_pair in core_meta}

    # Get assembly info and write content files
    assembly_info = get_assembly_info(core_dict["assembly.accession"])
    url_name = core_dict["species.url"]

    write_content(assembly_info, output_dir, url_name)
