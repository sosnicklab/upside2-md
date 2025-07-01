# Using Upside2 with Docker

This guide provides instructions for running the Upside2 application using Docker. This is the recommended method for most users as it provides a consistent, isolated environment without requiring manual dependency installation.

## Prerequisites

Before you begin, make sure you have Docker installed on your system. You can download it from the official [Docker website](https://www.docker.com/get-started).

## Docker Setup

To install and manage Docker manually, you can download it from the [official Docker website](https://www.docker.com/get-started). Once Docker is installed, follow the steps below.

### Getting the Docker Image

#### Option 1: Pull from Docker Hub (Recommended)

For most users, the easiest way to get started is to pull the pre-built image from Docker Hub:

```bash
docker pull oliverkleinmann/upside2-md
```

This image is regularly updated and contains all the necessary dependencies to run Upside2 simulations.

### Option 2: Build from Source

Building the image from source is only recommended if you need to add new dependencies that are not included in the default image. If you do add dependencies, please consider contributing by committing your changes to the `environment.yml` or `Dockerfile` and submitting a pull request on GitHub.

To build the image from source:
```bash
docker build -t upside2 https://github.com/sosnicklab/upside2-md.git
```

## Running the Container

### Non-Interactive Mode

To run a script inside the container and see the output, you can use non-interactive mode. For example, to run one of the example scripts:

```bash
docker run --rm oliverkleinmann/upside2-md python example/01.GettingStarted/0.run.py
```
The container will execute the command and then exit.

### Interactive Mode

For an interactive session inside the container, which is useful for exploration and running multiple commands:

```bash
docker run -it --rm oliverkleinmann/upside2-md
```
This command does the following:
*   `docker run`: Starts a new container.
*   `-it`: Opens an interactive terminal session.
*   `--rm`: Automatically removes the container when you exit.
*   `oliverkleinmann/upside2-md`: The name of the image to use.

### Mounting Custom Scripts

You can mount a directory of your own custom C++ scripts or other files into the container. This is useful for running your own c++ code within the container's pre-configured environment.

```bash
docker run -it --rm \
  -v "/path/to/your/scripts":/scripts \
  oliverkleinmann/upside2-md
```
Replace `/path/to/your/scripts` with the actual path to your scripts directory on your local machine. Your scripts will then be available inside the container at the `/scripts` directory.

## C++ Development

For C++ development, the expected IDE is VS Code with the Dev Containers extension.

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/sosnicklab/upside2-md.git
    cd upside2-md
    ```

2.  **Open in VS Code and Reopen in Container:**
    Open the `upside2-md` folder in VS Code. It will detect the `.devcontainer` configuration and prompt you to "Reopen in Container". This will build the development environment and mount your local repository files into it.

3.  **Configuring Other IDEs:**
    Support for other IDEs can be configured by modifying the `.devcontainer/devcontainer.json` file.

## Using Upside2 with Singularity on a Slurm Cluster

For users on HPC clusters that use the Slurm workload manager, you can run Upside2 jobs inside a Singularity container using manual Slurm job submission.

### Prerequisites

-   Singularity is installed on the cluster.
-   You have access to a Slurm cluster.

### Getting the Singularity Image

You can pull the image directly from Docker Hub:

```bash
singularity pull docker://oliverkleinmann/upside2-md
```

This command downloads the latest version of the Docker image and converts it into a Singularity Image File (`.sif`). The resulting file will be named `upside2-md_latest.sif`. You should rename it for convenience:

```bash
mv upside2-md_latest.sif upside2-md.sif
```

### Manual Job Submission

To submit jobs manually, create your own Slurm batch script. Here's an example template:

```bash
#!/bin/bash
#SBATCH --job-name=upside2_job
#SBATCH --time=4:00:00
#SBATCH --partition=cpu
#SBATCH --mem=16gb
#SBATCH --account=your-account

# Run your command inside the Singularity container
singularity exec upside2-md.sif python /path/to/your/script.py
```

Save this as a `.sbatch` file and submit it with:

```bash
sbatch your_job.sbatch
```


## Frequently Asked Questions (FAQ)

**Q: My project has a dependency that isn't included in the Docker/Singularity image. What should I do?**

**A:** You have two main options:

1.  **Extend the Image (Recommended for personal use):**
    You can create a new `Dockerfile` or Singularity definition file that uses the `oliverkleinmann/upside2-md` image as its base and adds your required dependencies. This is the best approach for project-specific needs.

    *   **For Docker:** Create a `Dockerfile` like this:
        ```Dockerfile
        FROM oliverkleinmann/upside2-md
        RUN conda install -y your-dependency
        ```
        Then build it: `docker build -t my-custom-upside2 .`

    *   **For Singularity:** Create a `.def` file to add your packages.

2.  **Contribute to the Main Image:**
    If you believe the dependency would be beneficial for many other users, please consider contributing to the project. You can do this by:
    *   Forking the [upside2-md repository](https://github.com/sosnicklab/upside2-md).
    *   Adding the dependency to the `environment.yml` file.
    *   Submitting a pull request with your changes.
