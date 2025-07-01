# .devcontainer/Dockerfile

FROM continuumio/miniconda3
LABEL maintainer="sosnicklab"

WORKDIR /upside2-md

RUN useradd -ms /bin/bash user
RUN apt-get update && apt-get install -y sudo && \
    echo 'user ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Install essential build tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    libhdf5-dev \
    libeigen3-dev \
    python3 \
    python3-dev \
    python3-pip \
    procps \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create a symlink to fix the hardcoded Eigen3 path in the build scripts.
# This makes the scripts' expected path point to the actual installation location.
RUN mkdir -p /home/user && \
    ln -s /usr/include/eigen3 /home/user/Eigen-3.3.5 && \
    ln -s /opt/conda/ /home/user/Anaconda3

# Clone the project repository into the working directory.
# This creates a self-contained image with the project code.
RUN git clone https://github.com/sosnicklab/upside2-md .

# Copy environment definition file into the container.
# This will overwrite the cloned version with the local one, allowing for easy customization.
COPY .devcontainer/environment.yml /tmp/environment.yml

# Add conda to the PATH
ENV PATH="/opt/conda/bin:$PATH"


# Configure conda for faster dependency resolution
RUN conda config --set channel_priority strict && \
    conda config --set solver libmamba

# Install upside from git. Developers may overwrite this in post-create.sh.
RUN ./install.sh

# Environment variables for the build process
ENV MY_PYTHON="/opt/conda"
ENV EIGEN_HOME="/usr/include/eigen3"
ENV UP_PATH="/upside2-md"
ENV UP_HOME="/upside2-md"

# Update conda and create the environment with a retry loop to handle network issues.
RUN conda update -n base -c defaults conda && \
    ( for i in 1 2 3; do \
        conda env create --name upside2-env -f /tmp/environment.yml && break; \
        echo "Conda env create failed, retrying ($i/3)..."; \
        sleep 5; \
      done; \
      if [ $? -ne 0 ]; then echo "Failed to create conda environment after 3 attempts."; exit 1; fi )

# Add conda init to the user's bashrc to ensure it's available in interactive shells.
RUN { \
    echo ''; \
    echo '# >>> conda initialize >>>'; \
    echo '# !! Contents within this block are managed by '\''conda init'\'' !!'; \
    echo '__conda_setup="$('\'/opt/conda/bin/conda\'' '\''shell.bash'\'' '\''hook'\'' 2> /dev/null)"'; \
    echo 'if [ $? -eq 0 ]; then'; \
    echo '    eval "$__conda_setup"'; \
    echo 'else'; \
    echo '    if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then'; \
    echo '        . "/opt/conda/etc/profile.d/conda.sh"'; \
    echo '    else'; \
    echo '        export PATH="/opt/conda/bin:$PATH"'; \
    echo '    fi'; \
    echo 'fi'; \
    echo 'unset __conda_setup'; \
    echo '# <<< conda initialize <<<'; \
    echo ''; \
    echo '# Activate the default conda environment.'; \
    echo 'conda activate upside2-env'; \
    echo ''; \
    echo '# Set project-specific environment variables.'; \
    echo 'export UPSIDE_HOME="/upside2-md"'; \
    echo 'export PATH="$UPSIDE_HOME/py:$UPSIDE_HOME/obj:$PATH"'; \
    echo 'export PYTHONPATH="$UPSIDE_HOME/py:$PYTHONPATH"'; \
} >> /home/user/.bashrc

# Change ownership of the workspace to the non-root user.
RUN chown -R user:user /upside2-md
# The only reason this is here is to not recompile for the nth time
COPY .devcontainer/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Create a directory for user-provided scripts and set ownership
RUN mkdir /upside2-md/scripts && chown user:user /upside2-md/scripts

# Switch to the non-root user. This will be the default for `docker run`.
USER user
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
