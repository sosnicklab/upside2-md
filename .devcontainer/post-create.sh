#!/bin/bash
# .devcontainer/post-create.sh

# Exit immediately if a command exits with a non-zero status.
set -e

echo "--- Starting Post-Create Setup ---"

# --- Build the C++ Code ---
echo "Setting up build environment for C++ compilation..."
export EIGEN_HOME=/usr/include/eigen3

echo "Building Upside C++ code..."
sudo ./install.sh

echo "--- Post-Create Setup Complete ---"

# --- Configure Shell for Interactive Use ---
echo "Configuring .bashrc for interactive shells..."
echo '' >> ~/.bashrc
echo '# >>> conda initialize >>>' >> ~/.bashrc
echo '# !! Contents within this block are managed by '\''conda init'\'' !!' >> ~/.bashrc
echo '__conda_setup="$('\'/opt/conda/bin/conda\'' '\''shell.bash'\'' '\''hook'\'' 2> /dev/null)"' >> ~/.bashrc
echo 'if [ $? -eq 0 ]; then' >> ~/.bashrc
echo '    eval "$__conda_setup"' >> ~/.bashrc
echo 'else' >> ~/.bashrc
echo '    if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then' >> ~/.bashrc
echo '        . "/opt/conda/etc/profile.d/conda.sh"' >> ~/.bashrc
echo '    else' >> ~/.bashrc
echo '        export PATH="/opt/conda/bin:$PATH"' >> ~/.bashrc
echo '    fi' >> ~/.bashrc
echo 'fi' >> ~/.bashrc
echo 'unset __conda_setup' >> ~/.bashrc
echo '# <<< conda initialize <<<' >> ~/.bashrc
echo '' >> ~/.bashrc
echo '# Activate the default conda environment' >> ~/.bashrc
echo 'conda activate upside2-env' >> ~/.bashrc
