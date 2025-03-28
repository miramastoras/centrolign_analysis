FROM ubuntu:22.04
MAINTAINER Mira Mastoras, mmastora@ucsc.edu

WORKDIR /app

# Update and install dependencies: GCC, CMake, Git
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    gcc \
    g++ \
    cmake \
    curl \
    pkg-config \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --recursive https://github.com/metaphysicser/RaMA.git

# Set the working directory to the cloned RaMA folder
WORKDIR RaMA
RUN mkdir build
WORKDIR build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release
RUN make -j$(nproc)

ENV PATH="app/RaMA/build/":$PATH
