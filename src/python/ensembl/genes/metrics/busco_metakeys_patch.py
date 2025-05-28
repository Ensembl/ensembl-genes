#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
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
import argparse
import json
import pymysql
import re
from pathlib import Path
from typing import Dict, Optional, Union


def parse_busco_file(file_path: str, db: str) -> Dict[str, Union[str, int]]:
    """
    Parses a BUSCO result file and extracts relevant data into a dictionary.

    Args:
        file_path (str): The path to the BUSCO result file.
        db(str): Core db name.

    Returns:
        Dict[str, str]: A dictionary containing parsed BUSCO data, including the dataset,
                        completeness values, and mode (proteins or genome).
    """

    # Declare the dictionary to accept str as keys and str or float as values
    data: Dict[str, Union[str, int]] = {}
    #data["core_db"] = db
    # Open and read the file
    with open(file_path, "r") as file:
        content = file.read()

    # Define regular expressions to match the relevant numbers
    version_pattern : Optional[re.Match[str]] = re.search(r"BUSCO version is: ((\d+\.\d+.\d+))", content)
    dataset_pattern : Optional[re.Match[str]] = re.search(r"The lineage dataset is: ([\w_]+)", content)
    mode_pattern: Optional[re.Match[str]] = re.search(r"BUSCO was run in mode: ([\w_]+)", content)
    completeness_pattern: Optional[re.Match[str]] = re.search(r"(\d+)\s+Complete BUSCOs \(C\)", content)
    single_copy_pattern: Optional[re.Match[str]] = re.search(
        r"(\d+)\s+Complete and single-copy BUSCOs \(S\)", content
    )
    duplicates_pattern: Optional[re.Match[str]] = re.search(
        r"(\d+)\s+Complete and duplicated BUSCOs \(D\)", content
    )
    fragmented_pattern: Optional[re.Match[str]] = re.search(r"(\d+)\s+Fragmented BUSCOs \(F\)", content)
    missing_pattern: Optional[re.Match[str]] = re.search(r"(\d+)\s+Missing BUSCOs \(M\)", content)

    # Initialize mode_match as None or str
    mode_match: Optional[str] = None

    # If match is not None, extract the group and assign it to mode_match
    if mode_pattern is not None:
        mode_match = mode_pattern.group(1)
    if mode_match in ("genome", "euk_genome_met", "euk_genome_min"):
        busco_mode = "genome"
    elif mode_match == "proteins":
        busco_mode = "protein"
    else:
        mode_match = None

    version = str(version_pattern.group(1)) if version_pattern else None
    dataset = str(dataset_pattern.group(1)) if dataset_pattern else None
    completeness = int(completeness_pattern.group(1)) if completeness_pattern else None
    single_copy = int(single_copy_pattern.group(1)) if single_copy_pattern else None
    duplicated = int(duplicates_pattern.group(1)) if duplicates_pattern else None
    fragmented = int(fragmented_pattern.group(1)) if fragmented_pattern else None
    missing = int(missing_pattern.group(1)) if missing_pattern else None

    # Extract the BUSCO summary line with completeness values
    if mode_match == "euk_genome_min":
        score_match = re.search(
            r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%.*,D:(\d+\.\d+)%\],F:(\d+\.\d+)%.*,M:(\d+\.\d+)%,n:(\d+),E:(\d+\.\d+)%",  # pylint: disable=line-too-long
            content,  # pylint: disable=line-too-long
        )
    else:
        score_match = re.search(
            r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%.*,D:(\d+\.\d+)%\],F:(\d+\.\d+)%.*,M:(\d+\.\d+)%,n:(\d+)", content
        )

    if score_match:
        score = score_match.group(0)
        total_buscos = score_match.group(6)
        if mode_match == "euk_genome_min":
            erroneus = score_match.group(7)

        if mode_match in ("genome", "euk_genome_met", "euk_genome_min"):
            # Extract the BUSCO version
            data["assembly.busco_version"] = str(version)
            # Extract the BUSCO dataset
            data["assembly.busco_dataset"] = str(dataset)
            # Store the BUSCO completeness summary with erroneous
            data["assembly.busco"] = str(score)
            data["assembly.busco_mode"] = busco_mode
            # Store the BUSCO values into individual fields
            data["assembly.busco_completeness"] = str(completeness)
            data["assembly.busco_single_copy"] = str(single_copy)
            data["assembly.busco_duplicated"] = str(duplicated)
            data["assembly.busco_fragmented"] = str(fragmented)
            data["assembly.busco_missing"] = str(missing)
            data["assembly.busco_total"] = int(total_buscos)
            if mode_match == "euk_genome_min":
                data["assembly.busco_erroneus"] = str(erroneus)

            data["assembly.busco"] = str(score)  # pylint: disable=line-too-long

        else:
            # Extract the BUSCO version
            data["genebuild.busco_version"] = str(version)
            # Extract the BUSCO dataset
            data["genebuild.busco_dataset"] = str(dataset)
            # Store the BUSCO completeness summary
            data["genebuild.busco"] = str(score)
            data["genebuild.busco_mode"] = busco_mode
            # Store the BUSCO values into individual fields
            data["genebuild.busco_completeness"] = str(completeness)
            data["genebuild.busco_single_copy"] = str(single_copy)
            data["genebuild.busco_duplicated"] = str(duplicated)
            data["genebuild.busco_fragmented"] = str(fragmented)
            data["genebuild.busco_missing"] = str(missing)
            data["genebuild.busco_total"] = int(total_buscos)
    return data


# Function to generate SQL patches
def generate_sql_patches(
        db_name: str, json_data: Dict[str, Union[str, float]], species_id: int = 1, table_name: str = "meta"
) -> str:  # pylint: disable=line-too-long
    """Creat Sql patch for database

    Args:
        db_name (str): db name
        json_data (Dict[str, Union[str, float]]): Dict of metakeys
        species_id (int, optional): species_id Defaults to 1.
        table_name (str, optional): Table name where to store data Defaults to "meta".

    Returns:
        str: list of Mysql patches
    """
    sql_statements = []
    sql_statements.append(f"USE {db_name};\n")  # Replace with your actual DB name

    # Iterate through the JSON key-value pairs
    for key, value in json_data.items():
        if value is None:
            # Skip if the value is None (or can handle it differently if needed)
            continue
        # Convert value to string and escape single quotes if necessary
        value_str = str(value).replace("'", "''")
        # Create the SQL INSERT statement
        sql_statements.append(
            f"INSERT IGNORE INTO {table_name} (species_id, meta_key, meta_value) VALUES ({species_id}, '{key}', '{value_str}');"  # pylint: disable=line-too-long
        )

    return "\n".join(sql_statements)


def process_busco_file(busco_file, db, output_dir,assembly_id=""):
    """
    Parses the BUSCO file, generates a JSON, writes it to an output file,
    and generates SQL patches.
    """
    # Parse the BUSCO file and generate the JSON
    busco_data = parse_busco_file(busco_file, db)

    # Determine the file name based on the mode (protein or genome)
    for key, value in busco_data.items():
        if key.endswith(".busco_mode"):
            busco_mode = value
            break  # Exit the loop once we find the first match

    output_file_name = f"{db}_busco_{busco_mode}_metakey.json"
    busco_data_json =busco_data.copy()
    busco_data_json['core_db']=db
    busco_data_json['assembly_id']=assembly_id
    # Convert the dictionary to a JSON object
    busco_json = json.dumps(busco_data_json, indent=4)

    # Write the JSON output to the dynamically named file
    output_path = Path(output_dir) / output_file_name
    with open(output_path, "w") as outfile:
        outfile.write(busco_json)

    # Output the JSON
    print(busco_json)

    # Generate SQL patches from the JSON
    sql_patches = generate_sql_patches(db, busco_data)

    # Return SQL patches to write them to an SQL file later
    return sql_patches

def execute_sql_patches(
    db_name: str,
    sql_statements: Dict[str, Union[str, float]],
    host: str,
    user: str,
    password: str,
    port: int
) -> str:  # pylint: disable=line-too-long
    """Create SQL patch for database and execute it

    Args:
        db_name (str): Database name
        sql_statements (Dict[str, Union[str, float]]): List of queryes
        host (str, optional): MySQL server host.
        user (str, optional): MySQL user.
        password (str, optional): MySQL password.
        port (int, optional): MySQL port.

    """
    sql_statements = sql_statements.strip().split(';\n')
    connection = None  # Initialize connection variable
    # Connect to the database and execute the SQL statements
    try:
        connection = pymysql.connect(
            host=host,
            user=user,
            password=password,
            database=db_name,
            port=int(port)
        )
        with connection.cursor() as cursor:
            for statement in sql_statements:
                print(statement.strip())
                statement=statement.strip()
                cursor.execute(statement)  # Execute each SQL statement
            connection.commit()  # Commit the changes
    except pymysql.MySQLError as e:
        print(f"Error while executing SQL: {e}")
    finally:
        connection.close()  # Close the database connection



def main():
    """
    Main function to parse a BUSCO result file and output the parsed data in JSON format.

    It expects the file path to the BUSCO result as a command-line argument.
    """

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Parse a BUSCO result file and generate JSON output.")
    parser.add_argument("-file", type=str, help="Path to the BUSCO result file")
    parser.add_argument("-db", type=str, help="Core db")
    parser.add_argument("-input_dir", type=str, help="Path for directory containing the busco output files")
    parser.add_argument("-output_dir", type=str, help="Path for output directory")
    parser.add_argument("-run_query", type=str, choices=['true', 'false'], help="Add busco metakeys to the db")
    parser.add_argument("-host", type=str, help="Server host")
    parser.add_argument("-port", type=str, help="Server port")
    parser.add_argument("-user", type=str, help="Db user with writable permission")
    parser.add_argument("-password", type=str, help="Server password")
    parser.add_argument("-assembly_id", type=str, help="Registry assembly id")
    # Parse arguments
    args = parser.parse_args()
    if args.file:
        # Process the single file and write to JSON and SQL
        sql_patches = process_busco_file(args.file, args.db, args.output_dir, args.assembly_id)
        with open(Path(args.output_dir) / f"{args.db}.sql", "a") as f:
            f.write(sql_patches)

    elif args.input_dir:
        # Process all files that end with 'busco_short_summary'
        busco_files = list(Path(args.input_dir).rglob("*busco_short_summary.txt"))

        with open(Path(args.output_dir) / f"{args.db}.sql", "a") as f:
            for file in busco_files:
                print(f"Processing file: {file}")
                sql_patches = process_busco_file(file, args.db, args.output_dir, args.assembly_id)
                f.write(sql_patches)
    if args.run_query == 'true':
        execute_sql_patches(args.db, sql_patches, args.host, args.user, args.password, int(args.port))


if __name__ == "__main__":
    main()

