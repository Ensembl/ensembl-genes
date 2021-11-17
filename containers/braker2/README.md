# BRAKER2 Docker image

BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database
Dockerfile and commands to create the BRAKER2 Docker image and run it as a Singularity (tag 2.1.5).

Links:
[source code](https://github.com/Gaius-Augustus/BRAKER),
[publication](https://academic.oup.com/nargab/article/3/1/lqaa108/6066535?login=true)

### Notes
The Docker image is not publicly available because a private licence is needed to run GeneMark. 
Download GeneMark software `gmes_linux_64.tar` and license `gm_key_64.gz` from ['GeneMark'](http://topaz.gatech.edu/GeneMark/license_download.cgi) and store them in the same directory with the Dockerfile. 

## Build and run the Docker image

Build the Docker image locally, using the Dockerfile
```
docker image build -t braker2:latest .
```
Upload the Docker image to DockerHub (you need an account and to log in with `docker login`) 
```
docker push <DockerHub username>/braker2:latest
```

Run the image (see [source code](https://github.com/Gaius-Augustus/BRAKER) to set up the command for BRAKER2)
```
docker run -v <working_dir> braker2:latest 
```

## Run the Docker image as a  Singularity 

Pull the Docker image from DockerHub in the `/hps/software/users/ensembl/genebuild/$USER/singularity` directory, it will be saved as a Singularity image
```
singularity pull docker://<DockerHub username>/braker2:latest
```

Submit a job to run the Singularity image on LFS
```
bsub -N 10 -M 4000 -R "rusage[mem=40000]" "singularity exec --bind <working_dir_absolute_path>:/app:rw /hps/software/users/ensembl/genebuild/$USER/singularity/braker2.sif braker.pl --genome=genome.fa --softmasking --esmode --core=10"
```
