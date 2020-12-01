ncores=$(getconf _NPROCESSORS_ONLN); njobs=$(( ncores > 1 ? (ncores - 1) : 1 ))

cp ./quasi_seq/cluster_interface/quasiseq/tasksMulti.py ./quasi_seq/cluster_interface/quasiseq/tasks.py

docker run --cpuset-cpus="0-${njobs}" --cpus "${ncores}" --memory 64g --env my_threads="${ncores}" --rm -i -t -v $(pwd)/output:/data -v $(pwd)/quasi_seq:/opt/quasi-seq -p 8080:9000 lhri/lhri_bioinformatics/quasi-seq:1.3 
	 
