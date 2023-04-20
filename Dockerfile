FROM ubuntu:latest

# Install Python and necessary packages
RUN apt-get update && \
    apt-get install -y python3 python3-pip python3-dev

# Install R and necessary packages
RUN apt-get update && \
    apt-get install -y r-base r-base-dev

# Install Julia
RUN apt-get update && \
    apt-get install -y julia

# Install build tools for R packages
RUN apt-get update && \
    apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

# Clone Git repository into container
RUN apt-get update && \
    apt-get install -y git && \
    git clone https://github.com/your-username/your-repo.git /app

# Set working directory
WORKDIR /app

# Compile C code for R package using Makefile
RUN make myfile.o

# Specify command to run when container starts
CMD ["bash"]
