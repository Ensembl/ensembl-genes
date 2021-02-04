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
################################################################################

# production image generation stage
################################################################################
FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# copy compiled binary from the builder stage to the production image
COPY --from=builder /app/Red/bin/Red /usr/local/bin

# specify a command to always run when running the image
ENTRYPOINT ["Red"]
################################################################################
