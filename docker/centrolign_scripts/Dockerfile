FROM ubuntu:22.04
MAINTAINER Mira Mastoras, mmastora@ucsc.edu

# Install Python 3 and necessary dependencies
RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Set Python 3 as the default
RUN ln -s /usr/bin/python3 /usr/bin/python

# Verify the Python version
RUN python --version

# Get centromere scripts
WORKDIR /opt
RUN git clone https://github.com/jeizenga/centromere-scripts.git

WORKDIR /data
