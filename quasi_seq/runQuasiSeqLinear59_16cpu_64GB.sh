#change celery concurrrency in init.sh to number cpus desired

data_dir=/mnt/LHRI/Bioinformatics/Projects/Quasispecies/Docker_runs/linear59_16cpu_64GB
mkdir -p $data_dir
cp ./cluster_interface/quasiseq/tasksLinear.py ./cluster_interface/quasiseq/tasks.py

docker run --cpuset-cpus="0-15" --cpus 16 --memory 64g --rm -i -t -v $data_dir:/data -v $(pwd):/opt/quasi-seq -p 8080:9000 lhri/lhri_bioinformatics/quasi-seq:1.3

#in docker container:
	#1) run /data/init.sh

	 
