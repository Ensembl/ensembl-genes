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
MySQL Helper Module

This module provides utility functions for connecting to a MySQL database,
fetching data, and performing updates. Logging is enabled to capture
informative messages and errors during execution.

Functions:
    - mysql_fetch_data: Executes a SELECT query and returns the result as a list of dictionaries.
    - mysql_update: Executes an UPDATE/INSERT/DELETE query with optional parameters.
"""

import logging
from typing import Optional, Any
import pymysql

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("pipeline_setup.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)

__all__ = ["mysql_fetch_data", "mysql_update"]


def mysql_get_connection(
    database: str, host: str, port: int, user: str, password: str
) -> Optional[pymysql.connections.Connection]:  # pylint: disable=unsubscriptable-object
    """
    Establish a connection to the MySQL database.
    """
    try:
        conn = pymysql.connect(
            host=host,
            user=user,
            port=port,
            password=password,
            database=database.strip(),
            cursorclass=pymysql.cursors.DictCursor,
        )
        return conn  # type: ignore
    except pymysql.Error as err:
        print(f"MySQL error: {err}")
        return None


def mysql_fetch_data(  # pylint:disable=too-many-arguments
    query: str,
    database: str,
    host: str,
    port: int,
    user: str,
    password: str,
    params: Optional[tuple[Any, ...] | list[Any]] = None,
) -> list[dict[str, Any]]:
    """
    Execute a SELECT query on the specified MySQL database and return results.

    Args:
        query (str): SQL SELECT query to execute.
        database (str): Name of the database.
        host (str): Host address of the MySQL server.
        port (int): Port number of the MySQL server.
        user (str): Username to connect to the database.
        password (str): Password for the user.
        params (tuple or list, optional): Parameters to safely interpolate into the query.

    Returns:
        list[dict]: A list of dictionaries representing query results.
                    Returns an empty list on failure.
    """
    try:
        conn = pymysql.connect(
            host=host,
            user=user,
            port=int(port),
            password=password,
            database=database.strip(),
            cursorclass=pymysql.cursors.DictCursor,
        )
        with conn.cursor() as cursor:
            cursor.execute(query, params or ())
            results = cursor.fetchall()
        conn.close()
        logger.info("Query successful")
        return list(results)

    except pymysql.Error as err:
        logger.error(  # pylint:disable=logging-fstring-interpolation
            f"MySQL error during fetch: {err}"
        )
        return []


def mysql_update(  # pylint:disable=too-many-arguments
    query: str,
    database: str,
    host: str,
    port: int,
    user: str,
    password: str,
    params: Optional[tuple[Any, ...] | list[Any]] = None,
) -> bool:
    """
    Execute an UPDATE, INSERT, or DELETE query on the specified MySQL database.

    Args:
        query (str): SQL query to execute.
        database (str): Name of the database.
        host (str): Host address of the MySQL server.
        port (int): Port number of the MySQL server.
        user (str): Username to connect to the database.
        password (str): Password for the user.
        params (tuple or list, optional): Parameters to safely interpolate into the query.

    Returns:
        bool: True if the query was executed successfully, False otherwise.
    """
    try:
        conn = pymysql.connect(
            host=host,
            user=user,
            port=port,
            password=password,
            database=database.strip(),
            cursorclass=pymysql.cursors.DictCursor,
        )
        with conn.cursor() as cursor:
            cursor.execute(query, params or ())
            conn.commit()
        conn.close()
        logger.info("Update successful.")
        return True

    except pymysql.Error as err:
        logger.error(  # pylint:disable=logging-fstring-interpolation
            f"MySQL error during update: {err}"
        )
        return False
