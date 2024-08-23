# scripts/metadata

## About

This script was used to produce the metadata sql patch for the given database, <database_name>.sql.

It will also produce a log file with any CRITICAL all WARNING messages, <database_name>_metadata.log, you should check this before applying the patch.

## Running core_meta_data.py

**python core_meta_data.py -h**

usage: core_meta_data.py [-h] [-o OUTPUT_DIR] -d DB_NAME -s HOST -p PORT

Prepare SQL updates for core dbs

optional arguments:
  -h, --help            show this help message and exit	

  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
     Path where the output and temp files will write to. Uses current dir by default

  -d DB_NAME, --db_name DB_NAME
     Database name

  -s HOST, --host HOST  Host server

  -p PORT, --port PORT  Host server port

**NOTE**

This script doesn't currently deal with collections, species_ids are hardcoded at **LN:207**. This should be updated, but for now, if you are dealing with a colleciotn db, please hardcode the species_id.