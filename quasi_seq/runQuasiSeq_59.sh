data_dir=/mnt/LHRI/Bioinformatics/Projects/Quasispecies/Docker_runs/multi59
mkdir -p $data_dir

docker run --cpus 4 --memory 16g --rm -i -t -v $data_dir:/data -v $(pwd):/opt/quasi-seq -p 8080:9000 lhri/lhri_bioinformatics/quasi-seq:1.3

#in docker container:
	#1) run /data/init.sh

	 
