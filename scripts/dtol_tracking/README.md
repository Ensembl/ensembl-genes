# DToL Tracking script

A script to check the status of Ensembl annotation on the DToL assemblies. It can also produce a list of databases that need to be added to the DToL projects page, https://projects.ensembl.org/darwin-tree-of-life/.

## To run

**python3 dtol_tracking.py -h**
\usage: dtol_tracking.py [-h] [--summary SUMMARY] [--in_progress IN_PROGRESS] [--unannotated UNANNOTATED] [--projects_missing PROJECTS_MISSING]

Track status of DToL assemblies in Ensembl.

optional arguments:
> -h, --help	show this help message and exit
> --summary SUMMARY      Provide a summary of the status of DToL assemblies in Ensembl
> --in_progress IN_PROGRESS      Provide a summary of the status of DToL assemblies in Ensembl
> --unannotated UNANNOTATED      Provide a list of GCAs for the DToL assemblies that have yet to be annotated.
> --projects_missing PROJECTS_MISSING    Provide a list of the core dbs for the Ensembl annotated DToL assemblies that have been release on rapid.ensembl.org
