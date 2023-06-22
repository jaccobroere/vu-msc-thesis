# This Dockerfile is meant to build the image for working and running code related to my MSc thesis.
# Run the following command to run the container and mount the correct folder.
# docker run -it --rm -v $(pwd):/app/vu-msc-thesis jaccusaurelius/vu-msc-thesis:workspace
# jaccusaurelius/vu-msc-thesis:kube can be used for docker with k8s
FROM julia:1.9-bullseye
ENV DEBIAN_FRONTEND=noninteractive
ENV PROJ_DIR=/app/vu-msc-thesis
ENV ZHU_DIR=/app/admm_src_zhu
ENV JULIA_DIR=/app/juliaenv

# Install necessary packages
RUN apt-get update && \
    apt-get install -y wget gnupg uuid-runtime locales && \
    sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen

ENV LC_ALL en_US.UTF-8 
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en

# Install Python and necessary packages
RUN apt-get update && \
    apt-get install -y python3 python3-pip

# Install R and necessary packages
RUN apt-get update && \
    apt-get install -y r-base r-base-dev

# Packages to get faster linear algebra
RUN apt-get install -y libatlas3-base libopenblas-base
# Debian binaries for R packages
RUN apt-get update && \ 
    apt-get install -y r-cran-data.table \
    r-cran-matrixstats \ 
    r-cran-glmnet \ 
    r-cran-igraph \
    r-cran-matrix

# Install build tools for R packages
RUN apt-get update && \
    apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

# Install R packages not available as binaries
COPY admm_src_zhu /app/admm_src_zhu
COPY splash_1.0.tar.gz /app/splash_1.0.tar.gz
COPY r-requirements.R /app/r-requirements.R
RUN cd /app && Rscript /app/r-requirements.R && cd ..
RUN cd /app/admm_src_zhu && \
    make clean && \
    make 

# Clone Git repository into container
RUN apt-get update && \
    apt-get install -y git
# git clone https://github.com/jaccobroere/vu-msc-thesis.git /app

# Then fix other dependencies
COPY python-requirements.txt /app/python-requirements.txt
COPY juliaenv /app/juliaenv
# # Install Python, R and Julia packages
RUN pip3 install -r /app/python-requirements.txt
RUN julia --project=$JULIA_DIR -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
RUN julia --project=$JULIA_DIR -e 'using Distributed, LoopVectorization, Tables, LinearAlgebra, Random, CSV, StatsBase, DataFrames, Distributions, SparseArrays, Statistics, Graphs, GraphIO, EzXML, MatrixMarket'
RUN ln -s /usr/bin/python3 /usr/bin/python

# Copy code to run
COPY src /app/vu-msc-thesis/src
COPY scripts /app/vu-msc-thesis/scripts

# Make output and app directory inside the container
RUN mkdir /app/vu-msc-thesis/data
RUN mkdir /app/vu-msc-thesis/out

# Set working directory
WORKDIR /app/vu-msc-thesis

# Specify command to run when container starts
CMD ["bash"]
