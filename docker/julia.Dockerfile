# This Dockerfile is meant to build the image for working and running code related to my MSc thesis.
# Run the following command to run the container and mount the correct folder.
# docker run -it --rm -v $(pwd):/app/vu-msc-thesis jaccusaurelius/vu-msc-thesis:latest
# jaccusaurelius/vu-msc-thesis:kube can be used for docker with k8s
FROM julia:1.8-bullseye
ENV DEBIAN_FRONTEND=noninteractive
ENV PROJ_DIR=/app/vu-msc-thesis
ENV ZHU_DIR=/app/admm_src_zhu
ENV JULIA_DIR=/app/juliaenv

# Install necessary packages
RUN apt-get update && \
    apt-get install -y wget && \
    apt-get install -y uuid-runtime && \
    apt-get install -y locales && \
    sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen

ENV LC_ALL en_US.UTF-8 
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en

# Install build tools for R packages
RUN apt-get update && \
    apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

# Clone Git repository into container
RUN apt-get update && \
    apt-get install -y git
#     git clone https://github.com/jaccobroere/vu-msc-thesis.git /app

# Then fix other dependencies
COPY juliaenv /app/juliaenv
RUN julia --project=$JULIA_DIR -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
RUN julia --project=$JULIA_DIR -e 'using Distributed, LoopVectorization, Tables, LinearAlgebra, Random, CSV, StatsBase, DataFrames, Distributions, SparseArrays, Statistics, Graphs, GraphIO, EzXML, MatrixMarket'


# Specify command to run when container starts
CMD ["bash"]
