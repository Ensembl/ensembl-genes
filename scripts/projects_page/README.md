# Projects Pages

Code for updating the ensembl projects pages

## GENESWICTH PROJECTS PAGE CURRENTLY HAS A CUSTOM FORMAT, DO NOT RUN THE SRIPT FOR THIS PROJECT - SEE LEANNE FOR DETAILS!

## Running write_yaml.py

**python write_yaml.py -h**

usage: write_yaml.py [-h] -f DB_FILE -p {aquafaang,bge,bovreg,cbp,dtol,erga,geneswitch,vgp}

Create species.yaml file for a given project page.

optional arguments:
  -h, --help            show this help message and exit
  -f DB_FILE, --db_file DB_FILE
                        Name for file containing list of VGP databases on the Rapid Release or Main server
  -p {aquafaang,bge,bovreg,cbp,dtol,erga,geneswitch,vgp}, --project {aquafaang,bge,bovreg,cbp,dtol,erga,geneswitch,vgp}
                        Name of the project this set of database belongs to

## Creating the DB_FILE
- The DB_FILE should contain a list of all databases that should be included on the project page
- You can create a list of database using BioPorject IDs or Taxon IDs, using [this script](https://github.com/Ensembl/ensembl-genes/blob/main/src/python/ensembl/genes/tracking/bioproject_tracking.py)

## Newly created species.yaml files

- When you run write_yaml.py you will create a species.yaml file for the given project.
- You should check this file, add it to [Ensembl projects repo](https://github.com/Ensembl/projects.ensembl.org/tree/master/_data) in the appropriate data folder and create a Pull Request for the new species.yaml file to be merged.

## Adding new image icons

- Images used for many of the project pages are stored in the [Ensembl projects repo](https://github.com/Ensembl/projects.ensembl.org/tree/master/) at /img/vgp/ 
- The file, icons.txt, in this repo contains a list of species classifications that have already been added to the pages
- If you are adding a species that does not fall into one of these classifications, then you will need to update the file 
   - check the species.classification keys in the meta table in the core db and choose a classification level that is appropriate
   - add the classification and image name to the icons.txt file
   - make sure the image exists in the [Ensembl projects repo](https://github.com/Ensembl/projects.ensembl.org/tree/master/) at /img/vgp/ - if it does not, you will need to add it (you should take the image from the rapid release, if it exists, otherwise, you should talk to someone from the Web team about getting a new appropriate image)

## Things to note in write_yaml.py

There are certain species for which the standard url format will not work and so custom urls must be created
```
      # there are species on main for which the upper case production name is used in the url instead of the upper case species name
      prod_url_list = ["bos_taurus_hybrid", "bos_indicus_hybrid"]
```

For some of the projects pages, there is an expectation that certain species/assemblies are listed first, so some ordering had to be added
```
    if project == "aquafaang":
        # move danio rerio reference dbs to end of the list
        for db in sorted_db_list:
            if "danio_rerio_core" in db:
                sorted_db_list.append(sorted_db_list.pop(sorted_db_list.index(db)))

    if project == "bovreg":
        # move bos taurus reference dbs to top of the list
        for db in sorted_db_list:
            if "bos_taurus_core" in db:
                sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

    if project == "geneswitch":
        # move sus scrofa or gallus gallus reference dbs to top of the list
        for db in sorted_db_list:
            if "sus_scrofa_core" in db or "gallus_gallus_core" in db:
                sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))
```

For the chicken reference (GENE-SWitCH project) the submitter information cannot be pulled from the assembly report so it has to be set
```
    # for the chicken reference the submitter info is missing, set it manually
    if info_dict["assembly.accession"] is 'GCA_000002315.5':
        submitter = "Genome Reference Consortium"
```

The GENE-SWitCH consortium want a release freeze on the data so all links must point to the 102 release
```
        if project == "geneswitch":
            release = "release-102"
            release_number = 102
```
NOTE: this will need to be changed when new annotations are added to this page e.g. the new VGP chickens (relese-102 data won't exists for these)