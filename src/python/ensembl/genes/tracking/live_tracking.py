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
"""the live_tracking module contains functions to track live core databases on staging servers."""

import argparse
import json
from typing import List, Tuple, Any
from pathlib import Path
import pymysql

with open(  # pylint:disable=unspecified-encoding
    Path(__file__).parent / "./live_tracking_config.json", "r"
) as f:
    config = json.load(f)


def mysql_fetch_data(
    query: str, database: str, host: str, port: int, user: str
) -> List[Tuple[Any, ...]]:
    """
    Executes a given SQL query on a MySQL database and fetches the result.

    This function establishes a connection to the MySQL database \
        using provided connection details.
    It then executes the given query using a cursor obtained from \
        the connection. After executing the query,
    it fetches all the rows of the query result and returns them. \
        The function handles any errors that might
    occur during the process and ensures that the database connection is \
        closed before returning the result.

    Args:
        query (str): The SQL query to be executed.
        database (str): The name of the database to connect to.
        host (str): The host name or IP address of the MySQL server.
        port (int): The port number to use for the connection.
        user (str): The username to use for the database connection.

    Returns:
        tuple: A tuple of tuples containing the rows returned by the query execution.

    Note:
        This function does not handle database password authentication. \
            Ensure that the provided user
        has the necessary permissions and that the database is configured\
            to allow password-less connections
        from the given host.
    """
    try:
        conn = pymysql.connect(
            host=host, user=user, port=port, database=database.strip()
        )

        cursor = conn.cursor()
        cursor.execute(query)
        result = cursor.fetchall()

        cursor.close()
        conn.close()

        return list(result)

    except pymysql.Error as err:
        print(f"Error: {err}")

        try:
            cursor.close()
            conn.close()
        except:  # pylint:disable=bare-except
            pass
        return []


def check_database_on_server(
    db: str, server_key: str, config: dict  # pylint:disable=redefined-outer-name
) -> bool:
    """
    Checks if a database exists on a given server.

    Args:
        db (str): The name of the database to check.
        server_key (str): The key of the server in the config.
        server_dict (dict): Dictionary containing server connection details.

    Returns:
        bool: True if the database exists, False otherwise.
    """
    try:
        conn = pymysql.connect(
            host=config["server_details"]["staging"][server_key]["db_host"],
            user=config["server_details"]["staging"][server_key]["db_user"],
            passwd=config["server_details"]["staging"][server_key]["db_pass"],
            port=config["server_details"]["staging"][server_key]["db_port"],
        )
        with conn.cursor() as cur:
            cur.execute(
                "SELECT SCHEMA_NAME FROM information_schema.SCHEMATA\
                    WHERE SCHEMA_NAME = %s",
                (db,),
            )
            result = cur.fetchone()
            return result is not None

    except pymysql.MySQLError as e:
        print(f"Error connecting to {server_key}: {e}")
        return False
    finally:
        if "conn" in locals() and conn:
            conn.close()


def clean_server(
    config: dict, mode: str  # pylint:disable=redefined-outer-name
) -> List:
    """
    Fetches a list of core databases that are both current in the \
        rapid database and present on the MySQL server.

    This function performs the following steps:
    1. Executes a query on the rapid database to fetch all current core\
        databases.
    2. Executes a query on the MySQL server to fetch all core databases.
    3. Compares the results from both queries and compiles a list of\
        core databases that are present in both.
        
    Args:
        config (dict): Configuration dictionary containing server connection details.
        mode (str): Mode of operation, either "rapid" or "beta".    

    Returns:
        list: A list of core databases that are current in the rapid \
            database and present on the MySQL server.
    """
    rapid_live_databases = []
    beta_live_databases = []

    cores_query = "SHOW DATABASES like '%core%';"
    cores_fetch = mysql_fetch_data(
        cores_query,
        "",
        config["server_details"]["genebuild"]["prod_1"]["db_host"],
        config["server_details"]["genebuild"]["prod_1"]["db_port"],
        config["server_details"]["genebuild"]["prod_1"]["db_user"],
    )
    core_dbs = [db[0] for db in cores_fetch]

    if mode == "rapid":  # pylint:disable=no-else-return
        rapid_query = (
            "SELECT dbname FROM genome_database JOIN genome USING(genome_id) "
            "JOIN data_release USING(data_release_id) WHERE \
                is_current = 1 and dbname like '%core%';"
        )
        rapid_fetch = mysql_fetch_data(
            rapid_query,
            config["server_details"]["meta"]["rapid"]["db_name"],
            config["server_details"]["meta"]["rapid"]["db_host"],
            config["server_details"]["meta"]["rapid"]["db_port"],
            config["server_details"]["meta"]["rapid"]["db_user"],
        )
        rapid_databases = [rrdb[0] for rrdb in rapid_fetch]

        for core_db in core_dbs:
            if core_db in rapid_databases:
                for server_key in ["st5", "st6"]:
                    if check_database_on_server(core_db, server_key, config):
                        rapid_live_databases.append(core_db)
                        break
        return rapid_live_databases

    elif mode == "beta":
        beta_query = (
            "SELECT DISTINCT(dataset_source.name) AS dbname FROM genome "
            "JOIN assembly ON genome.assembly_id = assembly.assembly_id "
            "JOIN genome_dataset ON genome.genome_id = genome_dataset.genome_id "
            "JOIN dataset ON genome_dataset.dataset_id = dataset.dataset_id "
            "JOIN dataset_source ON dataset.dataset_source_id = dataset_source.dataset_source_id "
            "JOIN genome_release ON genome.genome_id = genome_release.genome_id "
            "JOIN ensembl_release ON genome_release.release_id = ensembl_release.release_id "
            "WHERE dataset.name='assembly' AND dataset.status='Released';"
        )
        beta_fetch = mysql_fetch_data(
            beta_query,
            config["server_details"]["meta"]["beta"]["db_name"],
            config["server_details"]["meta"]["beta"]["db_host"],
            config["server_details"]["meta"]["beta"]["db_port"],
            config["server_details"]["meta"]["beta"]["db_user"],
        )
        beta_databases = [rrdb[0] for rrdb in beta_fetch]

        for core_db in core_dbs:
            if core_db in beta_databases:
                for server_key in ["st5", "st6"]:
                    if check_database_on_server(core_db, server_key, config):
                        beta_live_databases.append(core_db)
                        break
        return beta_live_databases

    else:
        raise ValueError("Mode must be either 'rapid' or 'beta'")


def main():
    """the main function to run the live tracking script."""
    parser = argparse.ArgumentParser(
        description="Check rapid or beta metadata databases."
    )
    parser.add_argument(
        "--mode",
        choices=["rapid", "beta"],
        required=True,
        help="Choose either 'rapid' or 'beta'",
    )
    args = parser.parse_args()

    live_databases = clean_server(config, args.mode)

    print(
        f"Here is a list of databases on genebuild-prod-1 that\
            can also be found in the {args.mode} metadata:"
    )
    for db in live_databases:
        print(f"{db}")


if __name__ == "__main__":
    main()
