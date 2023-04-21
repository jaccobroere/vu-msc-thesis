# This Dockerfile is meant to build the image for working and running code related to my MSc thesis.
# Run the following command to run the container and mount the correct folder.
# ocker run -it --rm -v $(pwd):/app/vu-msc-thesis jaccusaurelius/vu-msc-thesis:bullseye
FROM julia:1.8-bullseye
ENV DEBIAN_FRONTEND=noninteractive
ENV PROJ_DIR=/app/vu-msc-thesis
ENV ZHU_DIR=/app/admm_src_zhu

# Install necessary packages
RUN apt-get update && \
    apt-get install -y wget && \
    apt-get install -y locales && \
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

# Install build tools for R packages
RUN apt-get update && \
    apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

# Clone Git repository into container
RUN apt-get update && \
    apt-get install -y git
#     git clone https://github.com/jaccobroere/vu-msc-thesis.git /app

COPY python-requirements.txt /app/python-requirements.txt
COPY r-requirements.R /app/r-requirements.R
COPY juliaenv /app/juliaenv
COPY admm_src_zhu /app/admm_src_zhu

# Install Python, R and Julia packages
RUN pip3 install -r /app/python-requirements.txt
RUN Rscript /app/r-requirements.R
RUN julia --project=/app/juliaenv -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Install C dependencies for R package from Zhu et al. (2015)
RUN cd /app/admm_src_zhu && \
    make clean && \
    make

# Set working directory
WORKDIR /app/vu-msc-thesis

# Specify command to run when container starts
CMD ["bash"]
