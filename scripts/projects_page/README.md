# Projects Pages

Code for updating the ensembl projects pages

## Running write_yaml.py

**python write_yaml.py -h**

usage: write_yaml.py [-h] -f DB_FILE -p {aquafaang,bovreg,dtol,geneswitch,vgp}

Create species.yaml file for a given project page.

optional arguments:
  -h, --help            show this help message and exit
  -f DB_FILE, --db_file DB_FILE
                        Name for file containing list of VGP databases on the Rapid Release or Main server
  -p {aquafaang,bovreg,dtol,geneswitch,vgp}, --project {aquafaang,bovreg,dtol,geneswitch,vgp}
                        Name of the project this set of database belongs to

## Updating the DB_FILE

- The DB_FILEs in this repo (aquafaang_dbs.txt, bovreg_dbs.txt, dtol_dbs.txt, geneswitch_dbs.txt) should contain all the core databases for everything that is represented on the corresponding project page. 
- These databases can be found on the mirror server (mysql-ens-mirror-1) if on the main site or on the staging server (mysql-ens-sta-5) if on the rapid site.
- To add new annotations to a project page you must add the corresponding core database names to the apropriate DB_FILE.
- The database version must match the current live release version. If you run this script after a live release update you will need to update the versions of the exitsing core databases in the DB_FILE you are using (otherwise they will not be found on the servers).

NOTE:
Don't forget to push updated DB_FILEs to this repo.

## Newly created species.yaml files

- The species.yaml files in this repo (aquafaang_species.yaml, bovreg_species.yaml, dtol_species.yaml, geneswitch_species.yaml) should contain all the core databases for everything that is represented on the corresponding project page.
- When you run write_yaml.py you will create a species.yaml file for the given project - this will overwrite the existing species.yaml file.
- You should check this file, add it to [Ensembl projects repo](https://github.com/Ensembl/projects.ensembl.org/tree/master/_data) in the appropriate data folder and create a Pull Request for the new species.yaml file to be merged.

NOTE:
Don't forget to push updated species.yaml to this repo. 
