# Red Docker image

Dockerfile for generating a Docker image for Red, "an intelligent, rapid, accurate tool for detecting repeats de-novo on the genomic scale".

Links:
[source code](https://github.com/BioinformaticsToolsmith/Red)
[publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0654-5)


Docker image at Docker Hub:
https://hub.docker.com/r/williamebi/red


## build, publish, and run a Docker image

build a Docker image
```
docker image build -t williamebi/red:2.0 .
```

run the Docker image locally
```
# print the help message
docker run williamebi/red:2.0

# run repeat masking
docker run williamebi/red:2.0 -gnm <genome_directory> -msk <output_directory>
```

upload the Docker image to Docker Hub (you need an account and to log in with `docker login`)
```
docker push williamebi/red:2.0
```

run the Docker image directly with Singularity on a compute or login node
```
singularity run docker://williamebi/red:2.0 -gnm <genome_directory> -msk <output_directory>"
```

generate a Singularity image from the Docker image at Docker Hub, which will be saved in the `/hps/nobackup2/singularity/$USER` directory
```
singularity pull docker://williamebi/red:2.0
```

run the Singularity image on a compute or login node
```
singularity run "/hps/nobackup2/singularity/$USER/red-2.0.simg"
```

submit a job for running the Singularity image on LFS
```
bsub -I -M10000 -R"select[mem>10000] rusage[mem=10000]" "singularity run /hps/nobackup2/singularity/$USER/red-2.0.simg -gnm <genome_directory> -msk <output_directory>"
```
