FROM ubuntu:22.04
MAINTAINER Mira Mastoras, mmastora@ucsc.edu

# Install required packages
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y \
        gcc \
        g++ \
        cmake \
        git \
    && rm -rf /var/lib/apt/lists/*

# Verify installations
RUN gcc --version && cmake --version && git --version

# Get centrolign
WORKDIR /opt

RUN git clone https://github.com/jeizenga/centrolign.git
WORKDIR centrolign
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN make -j 8

ENV PATH="/opt/centrolign/build/":$PATH

# Set the default working directory
WORKDIR /data
