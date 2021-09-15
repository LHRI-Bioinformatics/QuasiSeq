# QuasiSeq

## Description

QuasiSeq is a pipeline for the deconvolution of closely related quasispecies from PacBio sequencing data.

All code for running QuasiSeq is available for download.  However, given the many system and software dependencies of the QuasiSeq pipeline, we have decided to create a [Docker] image in which all dependencies are resolved and example data sets are included.  The user should be able to run QuasiSeq on any system/operating system with a [Docker] host environment.

## Installation

Follow the [Docker installation instructions] for your system in order to use [Docker]

After installing [Docker], change to a directory on the [Docker] host for which the user has read/write access.  You will clone the repository to this directory and it will also contain QuasiSeq logs and pipeline output once the application is run.  

Once in the directory, clone the repository: `git clone https://github.com/LHRI-Bioinformatics/QuasiSeq.git`

Change to the resulting QuasiSeq directory and run the following script: `./runQuasiSeq.sh`

The script will prompt for a port to use for the QuasiSeq application (default port 8080).

Next, the script will prompt for a number of processors on the system to use for the QuasiSeq application (default is the total detected system processors).

The script will then prompt for the amount of memory in GB to use for the QuasiSeq application (default is the total detected available system memory).

Finally, the script will prompt for the path where you would like logging and output to be saved (default is the current directory and within a subdirectory name "output")

You should now be able to access the QuasiSeq application interface in a browser at `[Docker] host ip address:port entered` (e.g. 127.0.0.1:8080).
  
<!-- Next, log in to the gitlab registry: `docker login registry.gitlab.com`

Next, run the following command which will automatically pull down the QuasiSeq image and run a container with the QuasiSeq pipeline, interface, and example datasets:

`docker run --rm -i -t -v .:/data -p 9000:9000 registry.gitlab.com/lhri/lhri_bioinformatics/quasi-seq`

The  .:/data portion of the command mounts the /data directory within the QuasiSeq docker container to the current host directory(.).  The user may specify any directory to which they have read/write access.

The -p 9000:9000 portion of the command maps port 9000 on the host machine to port 9000 in the QuasiSeq docker container which is used for access to the django web interface for the pipeline.  The host port may be changed to any unused, open port (i.e. 8080:9000). -->


## Acknowledgments

Tom, Hiromi, Cliff, PacBio



## License
[Matlab Runtime] is needed to run the clustering code.  By downloading and running QuasiSeq, you are accepting Matlab's license.

[DAVID]: https://david.ncifcrf.gov/
[LHRI_GIT]: https://gitlab.com/LHRI/
[Matlab Runtime]: https://www.mathworks.com/products/compiler/matlab-runtime.html
[Docker]: https://www.docker.com/
[Docker installation instructions]: https://docs.docker.com/engine/installation/
