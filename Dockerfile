# This Dockerfile is meant to build the image for working and running code related to my MSc thesis.
# Run the following command to run the container and mount the correct folder.
# docker run -it --rm -v /path/to/local/vu-msc-thesis:/app/vu-msc-thesis my-image
FROM ubuntu:latest

# Install necessary packages
RUN apt-get update && \
    apt-get install -y wget

# Install Python and necessary packages
RUN apt-get update && \
    apt-get install -y python3 python3-pip

# Install R and necessary packages
RUN apt-get update && \
    apt-get install -y r-base r-base-dev

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz && \
    tar zxvf julia-1.8.5-linux-x86_64.tar.gz && \
    MV julia-1.8.5 /usr/local/julia-1.8.5 && \
    ln -s /usr/local/julia-1.8.5/bin/julia /usr/local/bin/julia

# Add Julia to system PATH
ENV PATH="/usr/local/julia-1.8.5/bin:${PATH}"

# Install build tools for R packages
RUN apt-get update && \
    apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

# Clone Git repository into container
# RUN apt-get update && \
#     apt-get install -y git && \
#     git clone https://github.com/jaccobroere/vu-msc-thesis.git /app

VOLUME /app/vu-msc-thesis

# Install Python, R and Julia packages
RUN pip3 install -r /app/vu-msc-thesis/python-requirements.txt
RUN Rscript /app/vu-msc-thesis/r-requirements.R
RUN julia --project=/app/vu-msc-thesis/juliaenv -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Install C dependencies for R package from Zhu et al. (2015)
RUN cd /app/vu-msc-thesis/admm_src_zhu && \
    make

# Set working directory
WORKDIR /app/vu-msc-thesis

# Specify command to run when container starts
CMD ["bash"]
