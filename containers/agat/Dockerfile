FROM ubuntu:16.04
LABEL maintainer "ftricomi@ebi.ac.uk"
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/London"

# Install dependencies
RUN apt-get update && apt-get install -y cmake build-essential git wget gzip perl cpanminus autoconf r-base
    
RUN cpanm install Clone Graph::Directed LWP::UserAgent Statistics::R JSON Carp Sort::Naturally File::Share File::ShareDir::Install Moose

RUN apt install -y libbio-perl-perl 
 		
RUN git clone https://github.com/NBISweden/AGAT.git 

RUN cd /AGAT/ && perl Makefile.PL && make  && make install && cd ..

RUN mkdir /data
VOLUME ["/data"]

WORKDIR /data
