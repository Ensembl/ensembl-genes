# Dockerfile for Red, "an intelligent, rapid, accurate tool for detecting repeats
# de-novo on the genomic scale"
# https://github.com/BioinformaticsToolsmith/Red

# base image to build upon
# using the latest Ubuntu LTS version
# https://hub.docker.com/_/ubuntu
FROM ubuntu:20.04

# re-synchronize package index files from their sources, install packages,
# and clean up the apt cache
RUN apt-get update && apt-get install -y \
    wget \
    && rm -rf /var/lib/apt/lists/*

# create and set a working directory
WORKDIR /app

# set up Red
################################################################################
# download and use the published binary
# http://toolsmith.ens.utulsa.edu/
# NOTE
# no HTTPS source available, this is a security issue
RUN wget http://toolsmith.ens.utulsa.edu/red/data/DataSet2Unix64.tar.gz
RUN tar -xvf DataSet2Unix64.tar.gz
RUN cp /app/redUnix64/Red /usr/local/bin
RUN rm -rf /app/*
################################################################################

# specify a command to always run when running the image
ENTRYPOINT ["Red"]
