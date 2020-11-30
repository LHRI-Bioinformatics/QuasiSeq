# Quasi-Seq

## Description

Quasi-Seq is a pipeline for the deconvolution of closely related quasispecies from PacBio sequencing data.

All code for running Quasi-Seq is available for download.  However, given the many system and software dependencies of the Quasi-Seq pipeline, we have decided to create a [Docker] image in which all dependencies are resolved and example data sets are included.  The user should be able to run Quasi-Seq on any system/operating system with a [Docker] host environment.

## Installation

Follow the [Docker installation instructions] for your system in order to use [Docker]

After installing [Docker], change to a directory on the [Docker] host for which the user has read/write access.  This directory will eventually contain Quasi-Seq logs and pipeline output.  

Next, log in to the gitlab registry: `docker login registry.gitlab.com`

Next, run the following command which will automatically pull down the Quasi-Seq image and run a container with the Quasi-Seq pipeline, interface, and example datasets:

`docker run --rm -i -t -v .:/data -p 9000:9000 registry.gitlab.com/lhri/lhri_bioinformatics/quasi-seq`

The  .:/data portion of the command mounts the /data directory within the Quasi-Seq docker container to the current host directory(.).  The user may specify any directory to which they have read/write access.

The -p 9000:9000 portion of the command maps port 9000 on the host machine to port 9000 in the Quasi-Seq docker container which is used for access to the django web interface for the pipeline.  The host port may be changed to any unused, open port (i.e. 8080:9000).


## Acknowledgments

Tom, Hiromi, Cliff, PacBio



## License
[Matlab Runtime] is needed to run the clustering code.  By downloading and running Quasi-Seq, you are accepting Matlab's license.

[DAVID]: https://david.ncifcrf.gov/
[LHRI_GIT]: https://gitlab.com/LHRI/
[Matlab Runtime]: https://www.mathworks.com/products/compiler/matlab-runtime.html
[Docker]: https://www.docker.com/
[Docker installation instructions]: https://docs.docker.com/engine/installation/
