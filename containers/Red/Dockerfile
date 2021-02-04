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
    git \
    make \
    g++-8 \
    && rm -rf /var/lib/apt/lists/*

# create and set a working directory
WORKDIR /app

# compile a Red binary from the source code
RUN git clone https://github.com/BioinformaticsToolsmith/Red.git
WORKDIR /app/Red/src_2.0
RUN make bin
RUN make
RUN cp /app/Red/bin/Red /usr/local/bin
RUN rm -rf /app/*

# specify a command to always run when running the image
ENTRYPOINT ["Red"]
