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

run the Docker image
```
docker run williamebi/red:2.0
```

upload the Docker image to Docker Hub (you need an account and to log in with `docker login`)
```
docker push williamebi/red:2.0
```

run the Docker image directly with Singularity on a compute or login node
```
singularity run docker://williamebi/red:2.0
```

generate a Singularity image from the Docker image at Docker Hub
```
singularity pull docker://williamebi/red:2.0
```

run the Singularity image on a compute or login node
```
singularity run /hps/nobackup2/singularity/william/red-2.0.simg
```

submit a job for running the Singularity image on LFS
```
bsub "singularity run /hps/nobackup2/singularity/william/red-2.0.simg"
```
