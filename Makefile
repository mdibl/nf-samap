visualize:
	# make docker
	# docker run --rm ryansonder/samap:latest visualize_samap.py -i data/samap_obj.pkl
	python scripts/visualize_samap.py -i data/samap_obj.pkl

run:
	make docker
	make clean
	nextflow run main.nf --with-docker

docker:
	docker build -t ryansonder/samap:latest .

clean:
	rm -rf results/
	rm -rf .vscode/nextflow.config
	rm -rf .vscode/
	rm -rf .nextflow/
	rm -rf work/
	rm -f .nextflow.log*
	rm -f nextflow.log
	rm -f nextflow.trace