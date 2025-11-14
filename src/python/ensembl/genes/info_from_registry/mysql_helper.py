"""
MySQL Helper Module

This module provides utility functions for connecting to a MySQL database,
fetching data, and performing updates. Logging is enabled to capture
informative messages and errors during execution.

Functions:
    - mysql_fetch_data: Executes a SELECT query and returns the result as a list of dictionaries.
    - mysql_update: Executes an UPDATE/INSERT/DELETE query with optional parameters.
"""

import pymysql
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline_setup.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def mysql_get_connection(database, host, port, user, password):
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
        return conn
    except pymysql.Error as err:
        print(f"MySQL error: {err}")
        return None


def mysql_fetch_data(query, database, host, port, user, password, params=None):
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
            port=port,
            password=password,
            database=database.strip(),
            cursorclass=pymysql.cursors.DictCursor,
        )
        with conn.cursor() as cursor:
            cursor.execute(query, params or ())
            results = cursor.fetchall()
        conn.close()
        logger.info("Query successful")
        return results

    except pymysql.Error as err:
        logger.error(f"MySQL error during fetch: {err}")
        return []


def mysql_update(query, database, host, port, user, password, params=None):
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
        logger.error(f"MySQL error during update: {err}")
        return False