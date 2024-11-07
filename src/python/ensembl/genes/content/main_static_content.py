import argparse
import pymysql
import requests
import xmltodict
from pathlib import Path


def mysql_fetch_data(query, database, host, port, user):
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
        conn = pymysql.connect(host=host, user=user, port=port, database=database.strip())
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
        dict: A dictionary containing assembly attributes such as name, level, submitter, and counts.
    """
    assembly_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    assembly_xml = requests.get(assembly_url)
    assembly_dict = xmltodict.parse(assembly_xml.text)
    assembly_attribs = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["ASSEMBLY_ATTRIBUTES"]["ASSEMBLY_ATTRIBUTE"]

    return_dict = {}

    # Retrieve specific attributes with error handling
    return_dict["assembly.name"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"].get("NAME", "")
    return_dict["assembly.level"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"].get("ASSEMBLY_LEVEL", "")
    return_dict["assembly.submitter"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["IDENTIFIERS"]["SUBMITTER_ID"].get("@namespace", "")

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


def write_content(assembly_info, output_dir, production_name):
    """
    Write assembly and annotation content to HTML files.

    Args:
        assembly_info (dict): Assembly information to be written.
        output_dir (Path): Directory path for output files.
        production_name (str): Name used to title the output files.
    """
    with open(output_dir / f"{production_name}_assembly.html", "w") as assembly_out:
        print(
            f"<p>The {assembly_info['assembly.name']} assembly was submitted by {assembly_info['assembly.submitter']} "
            f"and last updated on {assembly_info['assembly.date']}. The assembly is on the {assembly_info['assembly.level']} "
            f"level, consisting of {assembly_info['contig.count']} contigs assembled into {assembly_info['scaffold.count']} scaffolds. "
            f"The N50 size is the length such that 50% of the assembled genome lies in blocks of the N50 size or longer. "
            f"The N50 length for the contigs is {assembly_info['contig.n50']} while the scaffold N50 is {assembly_info['scaffold.n50']}.</p>",
            file=assembly_out
        )

    with open(output_dir / f"{production_name}_annotation.html", "w") as annotation_out:
        print(
            "<p>Genome annotation was generated using the "
            "<a href=\"https://beta.ensembl.org/help/articles/vertebrate-genome-annotation\">Ensembl vertebrate annotation pipeline</a>. "
            "</p><p>In accordance with the "
            "<a href=\"https://en.wikipedia.org/wiki/Fort_Lauderdale_Agreement\">Fort Lauderdale Agreement</a>, please check the publication "
            "status of the genome/assembly before publishing any genome-wide analyses using these data.</p>",
            file=annotation_out
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare static content files for main release")
    parser.add_argument("-o", "--output_dir", type=str, help="Output directory for files. Uses current dir by default.")
    parser.add_argument("-d", "--db_name", help="Database name", required=True)
    parser.add_argument("-s", "--host", help="Host server", required=True)
    parser.add_argument("-p", "--port", type=int, help="Host server port", required=True)

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
    core_query = (
        "SELECT meta_key, meta_value FROM meta WHERE meta_key IN "
        "('assembly.accession', 'species.production_name');"
    )
    core_meta = mysql_fetch_data(
        core_query,
        host=server_info["server"]["db_host"],
        user=server_info["server"]["db_user"],
        port=server_info["server"]["db_port"],
        database=db,
    )

    core_dict = {meta_pair[0]: meta_pair[1] for meta_pair in core_meta}
    upper_production_name = core_dict['species.production_name'].capitalize()

    # Get assembly info and write content files
    assembly_info = get_assembly_info(core_dict['assembly.accession'])
    write_content(assembly_info, output_dir, upper_production_name)
