# Gene Symbol Classifier eHive pipeline

eHive pipeline to run the Gene Symbol Classifier on an Ensembl core database.


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


## run tasks with Singularity

filter symbol assignments
```
SINGULARITY_IMAGE=<Singularity image path>; ASSIGNMENTS_DIRECTORY=<assignments directory path>; ASSIGNMENTS_CSV=<symbol assignments filename>; THRESHOLD=<assignment probability threshold>; singularity run --bind "$ASSIGNMENTS_DIRECTORY":/app/data "$SINGULARITY_IMAGE" --symbol_assignments "/app/data/$ASSIGNMENTS_CSV" --threshold "$THRESHOLD"
```

generate a symbol assignments CSV with gene description from an existing assignments CSV
```
SINGULARITY_IMAGE=<Singularity image path>; ASSIGNMENTS_DIRECTORY=<assignments directory path>; ASSIGNMENTS_CSV=<symbol assignments filename>; singularity run --bind "$ASSIGNMENTS_DIRECTORY":/app/data "$SINGULARITY_IMAGE" --symbol_assignments "/app/data/$ASSIGNMENTS_CSV"
```
