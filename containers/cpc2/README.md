# CPC2 Docker image

Dockerfile and commands to create the CPC2 Docker image and run it as a Singularity (tag 1.0.1).

Links:
[source code](https://github.com/gao-lab/CPC2_standalone),
[publication](https://academic.oup.com/nar/article/45/W1/W12/3831091)

Docker image at Docker Hub:
https://hub.docker.com/r/ftricomi/cpc2


## Build and run the Docker image

Build the Docker image locally, using the Dockerfile
```
docker image build -t ftricomi/cpc2:latest .
```
or pull the existing image from DockerHub 
```
docker pull ftricomi/cpc2
```

Run the image (see [source code](https://github.com/gao-lab/CPC2_standalone) to set up the command for CPC2)
```
docker -v <working_dir> run ftricomi/cpc2:latest python3 /CPC2_standalone-1.0.1/bin/CPC2.py -i input.fa --ORF -o out.txt
```

## Run the Docker image as a  Singularity 

Pull the Docker image from DockerHub in the `/hps/software/users/ensembl/genebuild/$USER/singularity` directory, it will be saved as a Singularity image
```
singularity pull docker://ftricomi/cpc2:latest
```

Submit a job for running the Singularity image on LFS
```
bsub -M 4000 -R "rusage[mem=40000]" "singularity exec --bind <working_dir_absolute_path>:/app:rw /hps/software/users/ensembl/genebuild/$USER/singularity/cpc2.sif python3 /CPC2_standalone-1.0.1/bin/CPC2.py -i input.fa --ORF -o out.txt"
```
