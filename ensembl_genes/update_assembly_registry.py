# Copyright [2018-2021] EMBL-European Bioinformatics Institute
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

import mysql.connector
from mysql.connector import Error
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import time
import argparse

# import traceback


def update_registry_db(query, database, host, port, user, password):
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
    # This method creates a dictionary for the lists from the sheets and makes the
    # dict key on the versioned GCA (which is unique). Once the dicts are generated,
    # It checks for the versioned GCA to update the status of the annotation
    existing_sheet_dict = {}

    # This ordering needs to match the ordering of the columns on the sheet
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
    for row in existing_sheet_records:
        gca = row[0]
        gca.encode("ascii", "ignore")

        if gca == "GCA":
            next
        else:
            existing_sheet_dict[gca] = row

    # This is where the majority of the work occurs. All assembly GCAs are examined to determine what
    # should be added/updated
    # Note that currently a three second sleep is needed to avoid exhausting the Sheets REST API quota
    # use creds to create a client to interact with the Google Drive API
    scope = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/drive",
    ]
    creds = ServiceAccountCredentials.from_json_keyfile_name(credentials_path, scope)
    client = gspread.authorize(creds)
    gettime = time.time()
    # Find a workbook by name and open the first sheet
    # Make sure you use the right name here.
    assembly_sheet = client.open(worksheet_name).worksheet("EnsemblAssemblyRegistry")
    annotation_status = "Completed"
    assembly_accession = accession
    sheet_row = existing_sheet_dict[assembly_accession]
    sheet_annotation_status_index = assembly_sheet_columns.index("Status")
    sheet_annotation_status_val = sheet_row[sheet_annotation_status_index]

    # Check status of the genebuild and update accordingly
    print("Updating genebuild status for: " + assembly_accession)
    row_update_index = assembly_sheet.find(assembly_accession).row
    update_cell_val(
        assembly_sheet,
        row_update_index,
        sheet_annotation_status_index,
        annotation_status,
    )


def update_cell_val(assembly_sheet, row_index, col_offset, val):
    col_offset += 1
    assembly_sheet.update_cell(row_index, col_offset, val)


def split_gca(gca):
    split_gca = gca.split(".")
    return split_gca


if __name__ == "__main__":

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
