# AGAT Docker image

Suite of tools to handle gene annotations in any GTF/GFF format.
Dockerfile and commands to create the AGAT Docker image and run it as a Singularity.

Links:
[source code](https://github.com/NBISweden/AGAT)

Docker image at Docker Hub:
https://hub.docker.com/r/ftricomi/agat


## Build and run the Docker image

Build the Docker image locally, using the Dockerfile
```
docker image build -t ftricomi/agat:latest .
```
or pull the existing image from DockerHub 
```
docker pull ftricomi/agat
```

Run the image (see [source code](https://github.com/NBISweden/AGAT) to choose the tool)
```
docker -v <working_dir> run ftricomi/agat:latest 
```

## Run the Docker image as a  Singularity 

Pull the Docker image from DockerHub in the `/hps/software/users/ensembl/genebuild/$USER/singularity` directory, it will be saved as a Singularity image
```
singularity pull docker://ftricomi/agat:latest
```

Submit a job for running the Singularity image on LFS
```
bsub -M 4000 -R "rusage[mem=40000]" "singularity exec --bind <working_dir_absolute_path>:/data:rw /hps/software/users/ensembl/genebuild/$USER/singularity/agat.simg agat_sp_extract_sequences.pl --gff input.gtf -f /data/genome.fa -p  -o /data/output.fa"
```
