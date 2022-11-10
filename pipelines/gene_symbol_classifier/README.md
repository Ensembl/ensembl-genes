# GSC eHive pipelines

eHive pipelines to run a [GSC](https://github.com/Ensembl/gene_symbol_classifier) classifier on an Ensembl core database, or Ensembl Rapid Release assemblies.


## assign symbols with eHive pipeline

initialize pipeline
```
init_pipeline.pl GeneSymbolClassifier_conf \
    --pipe_db_server <pipeline db server hostname> --pipe_db_port <pipeline db server port> \
    --user <username> --password <password> --user_r ensro \
    --pipe_db_name <pipeline database name> \
    --core_db_server_host <core db server hostname> --core_db_server_port <core db server port> \
    --core_db_name <core database name> \
    --annotation_data_directory <annotation data directory> \
    --singularity_image <Singularity image path> \
    --classifier_directory <saved classifier directory> \
    --classifier_filename <saved classifier filename> \
    --scientific_name <assembly scientific name> \
    --ehive_singularity_image <eHive pipeline Singularity image path> \
    --loading_threshold <loading symbols probability threshold>
```

run pipeline
```
beekeeper.pl --url $EHIVE_URL --loop
```


## generate pipeline Singularity image

build Docker image
```
docker image build --tag ensemblorg/gsc_pipeline:<image version> --file Dockerfile .
```

upload Docker image to Docker Hub
```
docker push ensemblorg/gsc_pipeline:<image version>
```

generate Singularity image from the Docker image
```
singularity pull docker://ensemblorg/gsc_pipeline:<image version>
```


## filter symbol assignments directly with Singularity

```
singularity run \
    --bind <assignments directory path>:/app/data \
    <Singularity image path> \
    --symbol_assignments /app/data/<symbol assignments filename> \
    --threshold <assignment probability threshold>
```
