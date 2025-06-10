visualize:
	docker run --rm -it \
		-v $(PWD):/workspace \
		-w /workspace \
		ryansonder/samap:latest \
		python scripts/visualize_samap.py -i data/samap_obj.pkl

run:
	make docker
	nextflow run main.nf --with-docker

docker:
	docker build -t ryansonder/samap:latest .

docker-shell:
	docker run --rm -it \
		-v $(PWD):/workspace \
		-w /workspace \
		--entrypoint /bin/bash \
		ryansonder/samap:latest


clean:
	rm -f .nextflow.log*
	rm -f nextflow.log
	rm -f nextflow.trace
	rm -rf results/
	rm -rf .vscode/nextflow.config
	rm -rf .vscode/
	rm -rf .nextflow/
	rm -rf work/