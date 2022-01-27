FROM ubuntu:20.04
LABEL maintainer "ftricomi@ebi.ac.uk"
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/London"

#Install dependencies
RUN apt-get update
RUN apt-get install -y -q \
    software-properties-common \
    libboost-iostreams-dev libboost-system-dev libboost-filesystem-dev \
    zlibc gcc-multilib apt-utils zlib1g-dev python \
    cmake tcsh build-essential g++ git wget gzip perl cpanminus autoconf

#GENE MARK
COPY gmes_linux_64.tar /
RUN tar -xf gmes_linux_64.tar   		

#AUGUSTUS 
RUN wget https://github.com/Gaius-Augustus/Augustus/archive/v3.4.0.tar.gz
RUN tar -xzf v3.4.0.tar.gz

RUN apt-get install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
RUN apt-get install -y libsqlite3-dev libmysql++-dev

# Install dependencies for bam2hints and filterBam 
RUN apt-get install -y libbamtools-dev

# Install additional dependencies for bam2wig
RUN apt-get install -y samtools libhts-dev

# Install additional dependencies for scripts
RUN apt-get install -y cdbfasta diamond-aligner libfile-which-perl libparallel-forkmanager-perl libyaml-perl libdbd-mysql-perl
RUN apt-get install -y --no-install-recommends python3-biopython

# Build AUGUSTUS
RUN cd /Augustus-3.4.0/ && make clean && make && make install && cd ..

#BAMTOOLS
RUN apt-get install bamtools -y -q

#DIAMOND
RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.0.6/diamond-linux64.tar.gz
RUN tar -xzf diamond-linux64.tar.gz

#PROTHINT
RUN wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.4.0/ProtHint-2.4.0.tar.gz
RUN tar -xzf ProtHint-2.4.0.tar.gz

#GENOME THREADER
RUN wget https://genomethreader.org/distributions/gth-1.7.3-Linux_x86_64-64bit.tar.gz
RUN tar -zxf gth-1.7.3-Linux_x86_64-64bit.tar.gz

RUN mkdir /ProtHint-2.4.0/bin/gth
RUN cp -r /gth-1.7.3-Linux_x86_64-64bit/bin/. /ProtHint-2.4.0/bin/gth/

#BRAKER2
RUN wget https://github.com/Gaius-Augustus/BRAKER/archive/v2.1.5.tar.gz
RUN tar -xzf v2.1.5.tar.gz

RUN cpanm File::Spec::Functions Hash::Merge List::Util Logger::Simple Module::Load::Conditional Parallel::ForkManager POSIX Scalar::Util::Numeric Exporter File::Basename Math::Utils YAML Cwd File::Glob File::Temp IO::File Test::More MCE::Mutex threads File::HomeDir Thread::Queue


RUN mkdir /data
VOLUME ["/data"]

WORKDIR /data

COPY gm_key_64.gz /
RUN  gunzip /gm_key_64.gz 
RUN  cp /gm_key_64 /.gm_key
RUN  cp /.gm_key /data/.gm_key
RUN  cp /.gm_key ~/.gm_key
RUN  echo ~/
ENV  HOME /home/
RUN  cp /.gm_key /$HOME/.gm_key   
RUN  chmod a+rwx -R /Augustus-3.4.0/config

RUN  cp -r /Augustus-3.4.0/auxprogs/utrrnaseq /Augustus-3.4.0/bin/utrrnaseq

ENV PATH="/BRAKER-2.1.5/scripts:${PATH}"
ENV PATH="/Augustus-3.4.0/bin:${PATH}"
ENV PATH="/Augustus-3.4.0/config:${PATH}"
ENV PATH="/Augustus-3.4.0/scripts:${PATH}"
ENV PATH="/gmes_linux_64/:${PATH}"
ENV PATH="/usr/bin:${PATH}"
ENV PATH="/ProtHint-2.4.0/bin:${PATH}"
ENV PATH="/Augustus-3.4.0/bin/utrrnaseq:${PATH}"
ENV PATH="/gth-1.7.3-Linux_x86_64-64bit/bin/:${PATH}"

ENV AUGUSTUS_CONFIG_PATH="/Augustus-3.4.0/config"
ENV AUGUSTUS_SCRIPTS_PATH="/Augustus-3.4.0/scripts"
ENV AUGUSTUS_BIN_PATH="/Augustus-3.4.0/bin/"
ENV GENEMARK_PATH="/gmes_linux_64/"
ENV BAMTOOLS_PATH="/usr/bin/"
ENV DIAMOND_PATH="/"
ENV PROTHINT_PATH="/ProtHint-2.4.0/bin"
ENV ALIGNMENT_TOOL_PATH="/gth-1.7.3-Linux_x86_64-64bit/bin/"
RUN echo $HOME  
RUN echo $PATH
RUN echo $WORKDIR

