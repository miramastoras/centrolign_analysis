FROM ubuntu:22.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
# ENV MPLCONFIGDIR=/tmp/matplotlib

# Before importing matplotlib in your script
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p /tmp/matplotlib && chmod 777 /tmp/matplotlib

# Fontconfig cache
ENV FONTCONFIG_CACHE=/tmp/fontconfig
RUN mkdir -p /tmp/fontconfig && chmod 777 /tmp/fontconfig

# Update system and install essentials
RUN apt-get update && apt-get install -y \
    wget bzip2 ca-certificates git build-essential \
    libssl-dev libbz2-dev libreadline-dev libsqlite3-dev \
    libffi-dev zlib1g-dev curl && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Pull environment.yml from GitHub (raw URL)
RUN wget -O /tmp/environment.yml https://raw.githubusercontent.com/miramastoras/centrolign_analysis/main/conda/tree_heatmap/environment-linux.yml

# Accept conda ToS
RUN conda config --set always_yes yes --set changeps1 no && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create the environment
RUN conda env create -f /tmp/environment.yml && rm /tmp/environment.yml

# Activate environment by default
# Activate environment by default
ENTRYPOINT ["conda", "run", "-n", "tree_python", "--no-capture-output"]

# Set working directory
WORKDIR /workspace
