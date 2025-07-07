FROM debian:bullseye-slim
LABEL maintainer="Oliver Kleinmann okleinmann@uchicago.edu"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    libhdf5-dev \
    libeigen3-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    python3-distutils \
    procps \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Anaconda to /opt/conda
RUN wget https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh && \
    bash Anaconda3-2023.09-0-Linux-x86_64.sh -b -p /opt/conda && \
    rm Anaconda3-2023.09-0-Linux-x86_64.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy
ENV PATH /opt/conda/bin:$PATH

# Create and activate a virtual environment
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    numpy \
    scipy \
    tables \
    pandas \
    h5py \
    pymbar \
    setuptools

# See mdtraj issue #1476 not yet updated in PyPi
RUN pip install --no-cache-dir git+https://github.com/mdtraj/mdtraj@master
# prody requires python3 --version <= 3.11
RUN pip install --no-cache-dir prody 
    
# Clone the Upside2-md repository
COPY . /upside2-md

# Specify Anaconda and eigen3 path
WORKDIR /upside2-md
RUN sed -i '3c\export MY_PYTHON=/opt/conda' source_sh
RUN sed -i '4c\export EIGEN_HOME=/usr/include/eigen3' source_sh

# Build Upside using the provided install.sh script
RUN chmod +x install.sh && ./install.sh

# Add the compiled Upside executable directory to PATH
ENV PATH="/upside2-md/obj:${PATH}"

CMD ["/bin/bash"]
