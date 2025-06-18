run:
	nextflow run main.nf --with-docker

run-noblast:
	nextflow run main.nf --use_precomputed_blast --with-docker

docker:
	docker build -f Dockerfile.samap -t pipeline/samap:latest .
	docker build -f Dockerfile.blast -t pipeline/samap-blast:latest .


rebuild:
	docker rmi pipeline/samap:latest || true
	docker build -t pipeline/samap:latest .


docker-shell-samap:
	docker run --rm -it \
		-v $(PWD):/workspace \
		-w /workspace \
		--entrypoint /bin/bash \
		pipeline/samap:latest

docker-shell-blast:
	docker run --rm -it \
		-v $(PWD):/workspace \
		-w /workspace \
		--entrypoint /bin/bash \
		pipeline/samap-blast:latest


clean-nextflow:
	rm -rf work/*
	rm -rf .nextflow*

clean-results:
	rm -rf results/*

clean-docker:
	docker rmi pipeline/samap:latest || true

clean-all: clean-nextflow clean-results