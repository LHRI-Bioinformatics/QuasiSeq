#!/bin/bash
#Reference the docker image to be used
# docker_image="lhri/lhri_bioinformatics/quasi-seq:1.3"
docker_image="ghcr.io/lhri-bioinformatics/quasiseq/quasiseq:1.3"

#Prompt for a port to use
myport=8080;
read -p "Enter an open port [${myport}]: " myportPrompted
myportPrompted=${myportPrompted:-${myport}}
myport=${myportPrompted}


#Prompt for the number of processes to use
ncores=$(getconf _NPROCESSORS_ONLN);
read -p "Enter the number of processes [${ncores}]: " ncoresPrompted
ncoresPrompted=${ncoresPrompted:-${ncores}}
ncores=${ncoresPrompted}

njobs=$(( ncores > 1 ? (ncores - 1) : 1 ))

#Prompt for the amount of memory to use in GB
memFree=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
memFree=$((memFree / 1024 / 1024))
read -p "Enter the amount of memory to use in GB [${memFree} GB]: " memFreePrompted
memFreePrompted=${memFreePrompted:-${memFree}}
memFree=${memFreePrompted}

#Prompt for the output path
outputDir="$(pwd)/output"
read -p "Enter the path where you would like logging and output to be saved [${outputDir}]: " outputDirPrompted
outputDirPrompted=${outputDirPrompted:-${outputDir}}
outputDir=${outputDirPrompted}

cp ./quasi_seq/cluster_interface/quasiseq/tasksMulti.py ./quasi_seq/cluster_interface/quasiseq/tasks.py

docker run --cpuset-cpus="0-${njobs}" --cpus "${ncores}" --memory "${memFree}"g --env "my_threads=${ncores}" --rm -i -t -v ${outputDir}:/data -v $(pwd)/quasi_seq:/opt/quasi-seq -p ${myport}:9000 ${docker_image}
