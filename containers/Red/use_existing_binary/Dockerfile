# Dockerfile for Red, "an intelligent, rapid, accurate tool for detecting repeats
# de-novo on the genomic scale"
# https://github.com/BioinformaticsToolsmith/Red

# build stage
################################################################################
# use the latest Ubuntu LTS version as the base image
# https://hub.docker.com/_/ubuntu
FROM ubuntu:20.04 as builder

# re-synchronize package index files from their sources, install packages,
# and clean up the apt cache
RUN apt-get update && apt-get install -y \
    wget \
    && rm -rf /var/lib/apt/lists/*

# create and set a working directory
WORKDIR /app

# download the published Red binary
# http://toolsmith.ens.utulsa.edu/
# NOTE
# no HTTPS source available, this is a security issue
RUN wget http://toolsmith.ens.utulsa.edu/red/data/DataSet2Unix64.tar.gz
RUN tar -xvf DataSet2Unix64.tar.gz
################################################################################

# production image generation stage
################################################################################
FROM ubuntu:20.04

# copy Red binary from the builder stage to the production image
COPY --from=builder /app/redUnix64/Red /usr/local/bin

# specify a command to always run when running the image
ENTRYPOINT ["Red"]
################################################################################
