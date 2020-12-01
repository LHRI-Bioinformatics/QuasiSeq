#change celery concurrrency in init.sh to number cpus desired

#data_dir=/mnt/LHRI/Bioinformatics/Projects/Quasispecies/Docker_runs/multi59_16cpu_64GB
#mkdir -p $data_dir
ncores=$(getconf _NPROCESSORS_ONLN); njobs=$(( ncores > 1 ? (ncores - 1) : 1 ))

cp ./quasi_seq/cluster_interface/quasiseq/tasksMulti.py ./quasi_seq/cluster_interface/quasiseq/tasks.py

docker run --cpuset-cpus="0-${njobs}" --cpus "${ncores}" --memory 64g --rm -i -t -v $(pwd)/output:/data -v $(pwd)/quasi_seq:/opt/quasi-seq -p 8080:9000 lhri/lhri_bioinformatics/quasi-seq:1.3

#in docker container:
	#1) run /data/init.sh

	 
