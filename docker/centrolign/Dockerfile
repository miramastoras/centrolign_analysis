FROM ubuntu:22.04
MAINTAINER Mira Mastoras, mmastora@ucsc.edu

# Set environment variables to prevent interactive prompts

# Update packages and install prerequisites
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    gcc \
    g++ \
    git \
    && rm -rf /var/lib/apt/lists/*

# Verify CMake and GCC versions
RUN cmake --version && gcc --version

# Get centrolign
WORKDIR /opt

RUN git clone https://github.com/jeizenga/centrolign.git && \
    cd centrolign && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 8 && \
    make install

ENV PATH=/opt/centrolign/build/centrolign:$PATH

# Set the default working directory
WORKDIR /data
