# scripts/stats

## About

This script is used to produce the stats sql patch for the given database, stats_<database_name>.sql.

It will also produce a log file with any error messages, <production_name>.log, you should check this before applying the patch.

## Running generate_species_homepage_stats.pl

**perl generate_species_homepage_stats.pl -dbname DB_NAME -host HOST -port PORT -production_name PRODUCTION_NAME -output_dir OUTPUT_DIR**

**NOTE**

This script will fail if the assembly name in the meta table does not match exactly the assembly name in the assembly report

e.g. for database , the value for assembly.name in the meta table is "PGSBv2.0_Norin61", but in the [assembly report](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/904/066/035/GCA_904066035.1_10wheat_assembly_norin61/GCA_904066035.1_10wheat_assembly_norin61_assembly_report.txt) it is "10wheat_assembly_norin61". The path to the report cannot be contructed with the incorrect name, so the script will fail.

**Solution**

If you know in advanced that your assembly name is different, or you find out after trying to run the script and it failing, add it to the core_meta_updates/scripts/stats/assembly_names.txt file in the format:

triticum_aestivum_norin61_core_59_112_1 10wheat_assembly_norin61 #thats a tab separator!


