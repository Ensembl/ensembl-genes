# Gene Symbol Classifier eHive pipeline

eHive pipeline to run the Gene Symbol Classifier on an Ensembl core database.


## generate Docker / Singularity image

build Docker image
```
docker image build --tag williamebi/gene_symbol_classifier_ehive_pipeline:<image version> --file Dockerfile .
```

assign gene symbols with a Docker container
```
ASSIGNMENTS_DIRECTORY=<assignments directory path>; ASSIGNMENTS_CSV=<symbol assignments filename>; THRESHOLD=<assignment probability threshold>; IMAGE_VERSION=<container image version>; docker run --read-only --volume="$ASSIGNMENTS_DIRECTORY":/app/data williamebi/gene_symbol_classifier_ehive_pipeline:"$IMAGE_VERSION" --symbol_assignments "/app/data/$ASSIGNMENTS_CSV" --threshold "$THRESHOLD"
```

upload the Docker image to Docker Hub
```
docker push williamebi/gene_symbol_classifier_ehive_pipeline:<image version>
```

generate a Singularity image from a Docker image at Docker Hub
```
singularity pull docker://williamebi/gene_symbol_classifier_ehive_pipeline:<image version>
```

assign gene symbols with a Singularity image
```
SINGULARITY_IMAGE=<Singularity image path>; ASSIGNMENTS_DIRECTORY=<assignments directory path>; ASSIGNMENTS_CSV=<symbol assignments filename>; THRESHOLD=<assignment probability threshold>; singularity run --bind "$ASSIGNMENTS_DIRECTORY":/app/data "$SINGULARITY_IMAGE" --symbol_assignments "/app/data/$ASSIGNMENTS_CSV" --threshold "$THRESHOLD"
```
