import json
import pymysql
from typing import List, Tuple, Any, Dict

with open("./live_tracking_config.json", "r") as f:
        config = json.load(f)

import pymysql
from typing import List, Tuple, Any
import sys

def mysql_fetch_data(query: str, database: str, host: str, port: int, user: str) -> List[Tuple[Any, ...]]:
    """
    Executes a given SQL query on a MySQL database and fetches the result.
    
    This function establishes a connection to the MySQL database using provided connection details.
    It then executes the given query using a cursor obtained from the connection. After executing the query,
    it fetches all the rows of the query result and returns them. The function handles any errors that might
    occur during the process and ensures that the database connection is closed before returning the result.

    Args:
        query (str): The SQL query to be executed.
        database (str): The name of the database to connect to.
        host (str): The host name or IP address of the MySQL server.
        port (int): The port number to use for the connection.
        user (str): The username to use for the database connection.

    Returns:
        tuple: A tuple of tuples containing the rows returned by the query execution.

    Note:
        This function does not handle database password authentication. Ensure that the provided user
        has the necessary permissions and that the database is configured to allow password-less connections
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
        
        return result

    except pymysql.Error as err:
        print(f"Error: {err}")
        
        try:
                cursor.close()
                conn.close()
        except:
                pass
        return []

            
def clean_server():
    """
    Fetches a list of core databases that are both current in the rapid database and present on the MySQL server.
    
    This function performs the following steps:
    1. Executes a query on the rapid database to fetch all current core databases.
    2. Executes a query on the MySQL server to fetch all core databases.
    3. Compares the results from both queries and compiles a list of core databases that are present in both.
    
    Returns:
        list: A list of core databases that are current in the rapid database and present on the MySQL server.
    """
    live_databases = []
    
    rapid_query =(
        "SELECT dbname FROM genome_database JOIN genome USING(genome_id) JOIN data_release USING(data_release_id) WHERE is_current = 1 and dbname like '%core%';"
    )
    rapid_fetch = mysql_fetch_data(
        rapid_query,
        config["server_details"]["meta"]["rapid"]["db_name"],
        config["server_details"]["meta"]["rapid"]["db_host"],
        config["server_details"]["meta"]["rapid"]["db_port"],
        config["server_details"]["meta"]["rapid"]["db_user"],
    )
    rapid_databases = [rrdb[0] for rrdb in rapid_fetch]
    
    cores_query = (
        "SHOW DATABASES like '%core%';"
    )
    cores_fetch = mysql_fetch_data(
        cores_query,
        '',
        config["server_details"]["genebuild"]["prod_1"]["db_host"],
        config["server_details"]["genebuild"]["prod_1"]["db_port"],
        config["server_details"]["genebuild"]["prod_1"]["db_user"],
    )
    
    for core_db in ([db[0] for db in cores_fetch]):
        if core_db in rapid_databases:
            live_databases.append(core_db)
            
    return live_databases

def main():
    drop_databases = clean_server()
                       
    print("Here is a list of databases on genebuild-prod-1 that can also be found in the current rapid release:")
    for db in drop_databases:
        print(f"{db}")
                    
if __name__ == "__main__":
    main()
