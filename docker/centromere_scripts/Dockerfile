FROM ubuntu:20.04
MAINTAINER Mira Mastoras, mmastora@ucsc.edu

# 1. Set base ubuntu env
RUN apt update --fix-missing && \
    DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
    apt upgrade -y && \
    apt install -y build-essential=12.8ubuntu1 default-jre=2:1.11-72 bedtools=2.27.1+dfsg-4ubuntu1 samtools=1.10-3 r-base=3.6.3-2 && \
    apt install -y vim htop git wget pigz

# update and install dependencies
RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ vim sudo build-essential bzip2 zlib1g-dev libbz2-dev \
        libcurl4-gnutls-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev

# Install Python 3 and necessary dependencies
RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN pip install pysam

#Install python, biopython
RUN apt-get -y update && \
 apt-get install -y python3-numpy python3-scipy python3-dev && \
 pip install biopython && \
 apt-get install -y libatlas-base-dev gfortran && \
 pip install svgwrite

RUN pip install lxml
RUN apt-get install -y libgl1
RUN pip install PyQt5=5.15.7

RUN apt-get install -y python-six python3-matplotlib && \
 pip install --upgrade ete3

# Set Python 3 as the default
RUN ln -s /usr/bin/python3 /usr/bin/python

# Verify the Python version
RUN python --version

RUN pip install pandas

# Get centromere scripts
WORKDIR /opt
RUN git clone https://github.com/miramastoras/centromere-scripts.git

# get samtools
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
    rm -r /opt/samtools/samtools-1.9.tar.bz2 && \
    cd samtools-1.9/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure && \
    make && \
    mkdir -p /root/bin/samtools_1.9 && \
    cp samtools /root/bin/samtools_1.9

# install mummer
ARG MUMMER_VER="4.0.0"

RUN wget https://github.com/mummer4/mummer/releases/download/v${MUMMER_VER}rc1/mummer-${MUMMER_VER}rc1.tar.gz && \
  tar -xvf mummer-${MUMMER_VER}rc1.tar.gz && \
  rm mummer-${MUMMER_VER}rc1.tar.gz && \
  cd mummer-${MUMMER_VER}rc1 && \
  ./configure && \
  make && \
  make install && \
  ldconfig && \
  cd .. && \
  mkdir /data

WORKDIR /data
