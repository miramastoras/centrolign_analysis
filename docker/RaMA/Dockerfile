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
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --recursive https://github.com/metaphysicser/RaMA.git /app

# Set the working directory to the cloned RaMA folder
WORKDIR /app
WORKDIR RaMA
WORKDIR build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release
RUN make -j$(nproc)
RUN ./RaMA -h
