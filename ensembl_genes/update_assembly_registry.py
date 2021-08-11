"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import mysql.connector
from mysql.connector import Error
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import argparse


def update_registry_db(query, database, host, port, user, password):
    '''Connection settings for registry database.'''
    try:
        conn = mysql.connector.connect(
            database=database, host=host, port=port, user=user, password=password
        )
        print(query)
        cursor = conn.cursor()
        cursor.execute(query)

        conn.commit()

    except Error as e:
        print(e)

    finally:
        cursor.close()
        conn.close()


def update_assembly_sheet(
    existing_sheet_records, assembly_sheet, worksheet_name, accession
):
    """This method creates a dictionary for the lists from the sheets and makes the dict key on the versioned GCA (which is unique). 
    Once the dicts are generated, it checks for the versioned GCA to update the status of the annotation
 
    Parameters:
        argument1 (list): A list of existing sheet entries.
        argument2 (str): Name of sheet.
        argument3 (str): Name of worksheet.
        argument4 (str): Assembly accession to be updated.

    Returns:
        str:Genebuild status update.
    """
    current_records = existing_sheet_records
    sheet_name = assembly_sheet
    wksht_name = worksheet_name
    genbank_accession = accession
    existing_sheet_dict = {}

    assembly_sheet_columns = [
        "GCA",
        "Clade",
        "Species name",
        "Common name",
        "Contig N50",
        "Assembly level",
        "Assembly date",
        "Assembly name",
        "RNAseq data",
        "RefSeq accession",
        "Genebuilder",
        "Status",
        "Assembly group",
        "Expected release",
        "Grant",
        "Notes",
        "Filter: Max version",
        "Filter: Genome rep",
        "Filter: N50",
        "Filter: Non-human",
    ]
    # This just makes a dict for the sheet based on the versioned GCA
    for row in current_records:
        gca = row[0]
        gca.encode("ascii", "ignore")

        if gca == "GCA":
            #do nothing
        else:
            existing_sheet_dict[gca] = row

    # Use creds to create a client to interact with the Google Drive API
    scope = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/drive",
    ]
    credentials = ServiceAccountCredentials.from_json_keyfile_name(credentials_path, scope)
    sheet_client = gspread.authorize(credentials)
    # Find a workbook by name and open the first sheet
    # Make sure you use the right name here.
    sheet_name = sheet_client.open(wksht_name).worksheet("EnsemblAssemblyRegistry")
    annotation_status = "Completed"
    sheet_row = existing_sheet_dict[genbank_accession]
    sheet_annotation_status_index = assembly_sheet_columns.index("Status")
    sheet_annotation_status_val = sheet_row[sheet_annotation_status_index]

    # Check status of the genebuild and update accordingly
    print("Updating genebuild status for: " + genbank_accession)
    row_update_index = sheet_name.find(genbank_accession).row
    update_cell_val(
        sheet_name,
        row_update_index,
        sheet_annotation_status_index,
        annotation_status,
    )


def update_cell_val(assembly_sheet, row_index, col_offset, val):
    """
    Updates assembly with genebuild status.

    Parameters:
        argument1 (str): The name of the sheet.
        argument2 (int): Sheet ssembly row index.
        argument3 (int): Sheet column index.
        argument4 (str): Genebuild status.

    Returns:
        The status of the genebuild.
    """
    assembly_sht_name = assembly_sheet
    row_id = row_index
    col_id_offset = col_offset
    value_update = val
    col_id_offset += 1
    assembly_sht_name.update_cell(row_id, col_id_offset, value_update)


def split_gca(gca):
    '''Split assembly accession to obtain chain and version values.'''
    accession = gca
    split_gca = accession.split(".")
    return split_gca


if __name__ == "__main__":
    main()


def main:
    '''Retrieve command line arguments and start the process of updating the shhet with genebuild status.'''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-as", "--assembly_accession", help="Assembly identifier", required=True
    )
    parser.add_argument(
        "-wsn",
        "--worksheet_name",
        help="The name of the Google Sheets worksheet",
        required=True,
    )
    parser.add_argument(
        "-gsc",
        "--gsheets_credentials",
        help="Path to a Google Sheets credentials JSON file for authentication",
        required=True,
    )
    parser.add_argument(
        "-ad",
        "--assembly_db_dbname",
        help="Name for assembly registry db",
        required=True,
    )
    parser.add_argument(
        "-ah", "--assembly_db_host", help="Host for assembly registry db", required=True
    )
    parser.add_argument(
        "-ap", "--assembly_db_port", help="Port for assembly registry db", required=True
    )
    parser.add_argument(
        "-au", "--assembly_db_user", help="User for assembly registry db", required=True
    )
    parser.add_argument(
        "-apw",
        "--assembly_db_password",
        help="Password for assembly registry db",
        required=True,
    )
    args = parser.parse_args()

    worksheet_name = args.worksheet_name
    credentials_path = args.gsheets_credentials
    gca = args.assembly_accession
    chain_version = split_gca(gca)

    assembly_db_database = args.assembly_db_dbname
    assembly_db_host = args.assembly_db_host
    assembly_db_port = args.assembly_db_port
    assembly_db_user = args.assembly_db_user
    assembly_db_password = args.assembly_db_password
    assembly_db_query = (
        "update assembly set annotated_status = 'completed' where chain = '"
        + chain_version[0]
        + "' and version = "
        + chain_version[1]
    )
    # use creds to create a client to interact with the Google Drive API
    scope = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/drive",
    ]

    creds = ServiceAccountCredentials.from_json_keyfile_name(credentials_path, scope)
    client = gspread.authorize(creds)
    # Find a workbook by name and open the first sheet
    # Make sure you use the right name here.
    assembly_sheet = client.open(worksheet_name).worksheet("EnsemblAssemblyRegistry")

    # Extract and print all of the values
    existing_sheet_records = assembly_sheet.get_all_values()

    update_assembly_sheet(existing_sheet_records, assembly_sheet, worksheet_name, gca)
    update_registry_db(
        assembly_db_query,
        assembly_db_database,
        assembly_db_host,
        assembly_db_port,
        assembly_db_user,
        assembly_db_password,
    )
