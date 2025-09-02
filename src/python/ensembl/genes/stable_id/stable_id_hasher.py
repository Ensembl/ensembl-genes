#!/usr/bin/env python3
"""
Script to find duplicate stable_id values across multiple MySQL databases
containing 'core' in their name.
"""

import mysql.connector
import argparse
import sys
from collections import defaultdict
from itertools import combinations

# Server configuration dictionary
SERVERS = {
    'sta1': {
        'host': 'mysql-ens-sta-1',
        'user': 'ensro',
        'password': '',
        'port': 4519,
        'name': 'sta1'
    },
    'sta5': {
        'host': 'mysql-ens-sta-5',
        'user': 'ensro',
        'password': '',
        'port': 4684,
        'name': 'sta5'
    },
    'sta6': {
        'host': 'mysql-ens-sta-6',
        'user': 'ensro',
        'password': '',
        'port': 4695,
        'name': 'sta6'
    },
    'gb6': {
        'host': 'mysql-ens-genebuild-prod-6',
        'user': 'ensro',
        'password': '',
        'port': 4532,
        'name': 'gb6'
    },
    'gb5': {
        'host': 'mysql-ens-genebuild-prod-5',
        'user': 'ensro',
        'password': '',
        'port': 4531,
        'name': 'gb5'
    },
}

def get_databases_with_keywords(cursor, keywords):
    """Get all database names containing all specified keywords"""
    if isinstance(keywords, str):
        keyword_list = [k.strip() for k in keywords.split(',')]
    else:
        keyword_list = keywords
    
    # Build LIKE pattern: %keyword1%keyword2%keyword3%
    pattern = '%' + '%'.join(keyword_list) + '%'
    
    cursor.execute(f"SHOW DATABASES LIKE '{pattern}'")
    databases = cursor.fetchall()
    return [db[0] for db in databases]

def get_table_with_stable_id(cursor, database):
    """Check if the gene table exists and has stable_id column"""
    cursor.execute(f"USE `{database}`")
    
    # Check if gene table exists
    cursor.execute("SHOW TABLES")
    tables = [table[0] for table in cursor.fetchall()]
    
    if 'gene' not in tables:
        return False
    
    # Check if stable_id column exists
    cursor.execute(f"DESCRIBE `gene`")
    columns = [col[0] for col in cursor.fetchall()]
    
    return 'stable_id' in columns

def get_stable_ids(cursor, database):
    """Get all stable_id values from the gene table"""
    cursor.execute(f"USE `{database}`")
    cursor.execute("SELECT stable_id FROM gene")
    return [row[0] for row in cursor.fetchall() if row[0] is not None]

def find_duplicate_stable_ids(server_list, db_keywords='core'):
    """Main function to find duplicate stable_ids across databases matching keywords"""
    
    stable_id_map = defaultdict(list)
    
    print(f"Processing {len(server_list)} servers...")
    
    for server_config in server_list:
        try:
            # Connect to MySQL server
            connection = mysql.connector.connect(
                host=server_config['host'],
                user=server_config['user'],
                password=server_config['password'],
                port=server_config['port']
            )
            cursor = connection.cursor()
            
            print(f"Connected to MySQL server {server_config['name']} at {server_config['host']}:{server_config['port']}")
            
            # Get databases with keywords in name
            matching_databases = get_databases_with_keywords(cursor, db_keywords)
            
            print(f"Scanning {len(matching_databases)} databases on {server_config['name']} matching pattern '%{db_keywords.replace(',', '%')}%'...\n")
            
            # Process each database
            processed_dbs = 0
            total_dbs = len(matching_databases)
            
            for i, database in enumerate(matching_databases, 1):
                try:
                    # Check if table exists and has stable_id column
                    if not get_table_with_stable_id(cursor, database):
                        print(f"{i}/{total_dbs} WARNING: Skipping {server_config['name']}:{database}: table 'gene' not found or no 'stable_id' column")
                        continue
                    
                    # Get stable_ids from this database
                    stable_ids = get_stable_ids(cursor, database)
                    print(f"{i}/{total_dbs} {server_config['name']}:{database}: {len(stable_ids)} stable_ids")
                    
                    # Add to our mapping with server prefix
                    for stable_id in stable_ids:
                        stable_id_map[stable_id].append(f"{server_config['name']}:{database}")
                    
                    processed_dbs += 1
                    
                except mysql.connector.Error as e:
                    print(f"{i}/{total_dbs} ERROR: Error processing database {server_config['name']}:{database}: {e}")
                    continue
            
            print(f"{'='*60}")
            print(f"SUMMARY for {server_config['name']}: Processed {processed_dbs}/{total_dbs} databases")
            print(f"{'='*60}\n")
            
        except mysql.connector.Error as e:
            print(f"MySQL connection error for {server_config['name']}: {e}")
            continue
            
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
                print(f"MySQL connection to {server_config['name']} closed.\n")
    
    # Find duplicates (stable_ids that appear in more than one database)
    duplicates = {k: v for k, v in stable_id_map.items() if len(v) > 1}
    
    if duplicates:
        print(f"\nFOUND {len(duplicates)} DUPLICATE stable_id VALUES")
        
        # Write detailed duplicates file
        with open('duplicates_detailed.txt', 'w') as f:
            f.write(f"DUPLICATE stable_id VALUES FOUND\n")
            f.write(f"{'='*60}\n\n")
            
            for stable_id, databases in duplicates.items():
                f.write(f"stable_id: '{stable_id}'\n")
                f.write(f"  Found in {len(databases)} databases:\n")
                for db in databases:
                    f.write(f"    - {db}\n")
                f.write(f"\n")
        
        print("Detailed duplicates written to: duplicates_detailed.txt")
        
        # Generate database pairs and write pairs file
        db_pairs = set()
        for stable_id, databases in duplicates.items():
            # Generate all combinations of database pairs for this stable_id
            for pair in combinations(databases, 2):
                db_pairs.add(tuple(sorted(pair)))
        
        with open('duplicates_pairs.txt', 'w') as f:
            f.write(f"DATABASE PAIRS SHARING DUPLICATE stable_ids\n")
            f.write(f"{'='*50}\n\n")
            
            for db1, db2 in sorted(db_pairs):
                f.write(f"{db1}, {db2}\n")
        
        print(f"Database pairs written to: duplicates_pairs.txt ({len(db_pairs)} pairs)")
        
        # Print summary stats
        print(f"\n{'='*60}")
        print(f"Total unique stable_ids: {len(stable_id_map)}")
        print(f"Duplicate stable_ids: {len(duplicates)}")
        print(f"Clean stable_ids: {len(stable_id_map) - len(duplicates)}")
        print(f"{'='*60}")
        
    else:
        print("\nNO DUPLICATES FOUND!")
        print(f"All {len(stable_id_map)} stable_id values are unique across databases.")

def main():
    parser = argparse.ArgumentParser(
        description='Find duplicate stable_id values across MySQL databases containing specified keywords',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python script.py -s prod1
  python script.py -s prod1,prod5,prod6 -k erebo,core
  python script.py -s prod5 -k test,dev,core
        """
    )
    
    parser.add_argument('-s', '--servers', required=True,
                       help='Server aliases to process, comma-separated (available: ' + 
                            ', '.join(SERVERS.keys()) + ')')
    parser.add_argument('-k', '--keywords', default='core',
                       help='Keywords to search for in database names, comma-separated (default: core)')
    
    args = parser.parse_args()
    
    # Parse server aliases
    server_aliases = [s.strip() for s in args.servers.split(',')]
    server_list = []
    
    for alias in server_aliases:
        if alias not in SERVERS:
            print(f"ERROR: Unknown server alias '{alias}'. Available servers: {', '.join(SERVERS.keys())}")
            sys.exit(1)
        server_list.append(SERVERS[alias])
    
    find_duplicate_stable_ids(server_list, args.keywords)

if __name__ == "__main__":
    main()
