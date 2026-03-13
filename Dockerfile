# Use Python 3.9.16 base image
FROM python:3.9.16-slim

# Set environment variables to take away user interaction in install
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Clone the repository
RUN git clone --depth 1 --branch preprint_keras https://github.com/nvk23/CNV-Finder.git /opt/cnv-finder

# Set working directory
WORKDIR /opt/cnv-finder

# Install Python requirements
RUN pip install --no-cache-dir -r requirements.txt

# Set default shell entry point
ENTRYPOINT ["/bin/sh"]