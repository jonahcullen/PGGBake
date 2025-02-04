FROM ubuntu:20.04

# labels
LABEL version="2.3"
LABEL description="Mash 2.3 image for sequence analysis"
LABEL org.opencontainers.image.ref.name="ubuntu"
LABEL org.opencontainers.image.version="20.04"

# install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    wget \
    zlib1g-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# define mash install directory
ENV MASH_INSTALL_DIR=/opt/mash

# download and install mash v2.3
WORKDIR /tmp
RUN wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar && \
    tar -xf mash-Linux64-v2.3.tar && \
    mkdir -p ${MASH_INSTALL_DIR}/bin && \
    mv mash-Linux64-v2.3/mash ${MASH_INSTALL_DIR}/bin/ && \
    rm -rf mash-Linux64-v2.3 mash-Linux64-v2.3.tar

# add mash to PATH
ENV PATH="${MASH_INSTALL_DIR}/bin:${PATH}"

# set working directory
WORKDIR /data

# default command
ENTRYPOINT ["mash"]
