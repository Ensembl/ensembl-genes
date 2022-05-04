# GSC eHive pipelines

eHive pipelines to run a [GSC](https://github.com/Ensembl/gene_symbol_classifier) classifier on an Ensembl core database, or Ensembl Rapid Release assemblies.


## standalone pipeline

initialize the pipeline
```
init_pipeline.pl GeneSymbolClassifier_conf --pipe_db_server <pipeline MySQL server hostname> --pipe_db_port <pipeline MySQL server port> --user <username> --password <password> --user_r ensro --pipe_db_name <pipeline database name> --core_db_server_host <core database server hostname> --core_db_server_port <core database server port> --core_db_name <core database name> --annotation_data_directory <annotation data directory> --singularity_image <Singularity image path> --classifier_directory <saved classifier directory> --classifier_filename <saved classifier filename> --scientific_name <assembly scientific name> --ehive_singularity_image <eHive pipeline Singularity image path> --loading_threshold <loading symbols probability threshold>
```


## generate Singularity image

build Docker image
```
docker image build --tag williamebi/gene_symbol_classifier_ehive_pipeline:<image version> --file Dockerfile .
```

upload Docker image to Docker Hub
```
docker push williamebi/gene_symbol_classifier_ehive_pipeline:<image version>
```

generate Singularity image from the Docker image
```
singularity pull docker://williamebi/gene_symbol_classifier_ehive_pipeline:<image version>
```


## run specific tasks with Singularity

filter symbol assignments
```
singularity run --bind <assignments directory path>:/app/data <Singularity image path> --symbol_assignments /app/data/<symbol assignments filename> --threshold <assignment probability threshold>
```
