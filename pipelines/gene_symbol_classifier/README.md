# Gene Symbol Classifier eHive pipeline

eHive pipeline to run a trained [Gene Symbol Classifier](https://github.com/Ensembl/gene_symbol_classifier) network for an Ensembl core database.


## run the pipeline

initialize the pipeline
```
init_pipeline.pl gene_symbol_classifier_conf --pipe_db_server <pipeline MySQL server hostname> --pipe_db_port <pipeline MySQL server port> --user <username> --password <password> --user_r ensro --pipe_db_name <pipeline database name> --core_db_server_host <core database server hostname> --core_db_server_port <core database server port> --core_db_name <core database name> --annotation_data_directory <annotation data directory> --singularity_image <singularity image path> --classifier_directory <classifier checkpoint directory> --classifier_filename <classifier checkpoint filename> --scientific_name <assembly scientific name> --ehive_singularity_image <eHive pipeline singularity image path> --loading_threshold <loading symbols threshold probability>
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
SINGULARITY_IMAGE=<Singularity image path>; ASSIGNMENTS_DIRECTORY=<assignments directory path>; ASSIGNMENTS_CSV=<symbol assignments filename>; THRESHOLD=<assignment probability threshold>; singularity run --bind "$ASSIGNMENTS_DIRECTORY":/app/data "$SINGULARITY_IMAGE" --symbol_assignments "/app/data/$ASSIGNMENTS_CSV" --threshold "$THRESHOLD"
```

generate a symbol assignments CSV with gene description from an existing assignments CSV
```
SINGULARITY_IMAGE=<Singularity image path>; ASSIGNMENTS_DIRECTORY=<assignments directory path>; ASSIGNMENTS_CSV=<symbol assignments filename>; singularity run --bind "$ASSIGNMENTS_DIRECTORY":/app/data "$SINGULARITY_IMAGE" --symbol_assignments "/app/data/$ASSIGNMENTS_CSV"
```
