FROM ubuntu:20.04
LABEL maintainer "ftricomi@ebi.ac.uk"
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/London"

# Install dependencies
RUN apt-get update && apt-get upgrade -y -q
RUN apt-get update && apt-get install -y \
    software-properties-common
RUN add-apt-repository universe
RUN apt-get update && apt-get install -y \
    apache2 \
    curl \
    git \
    python3.4 \
    python3-pip
RUN apt-get install -y     build-essential     ca-certificates     gcc     git     libpq-dev     make \     
         ssh     && apt-get autoremove     && apt-get clean    
RUN apt-get autoremove && apt-get clean   
RUN pip install --no-cache-dir 'biopython' 
RUN pip install --no-cache-dir 'numpy' 

#Download CPC2 repo
RUN wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz
RUN gzip -dc v1.0.1.tar.gz | tar xf -
RUN cd /CPC2_standalone-1.0.1/libs/libsvm/ && gzip -dc libsvm-3.18.tar.gz | tar xf - && cd libsvm-3.18 && make clean && make 

ENV PATH="/CPC2_standalone-1.0.1/bin/:${PATH}"
ENV CPC_HOME="/CPC2_standalone-1.0.1/bin/"

RUN cd $CPC_HOME     && sed -i "s/if not os.path.isfile(options.fasta):/if (not os.path.isfile(options.fasta) and os.path.islink(options.fasta)):/" CPC2.py

RUN mkdir /data
VOLUME ["/data"]

WORKDIR /data
 
ENTRYPOINT ["python3"] 
