run:
	nextflow run main.nf --with-docker


docker:
	docker build -f Dockerfile.samap -t ryansonder/samap:latest .
	docker build -f Dockerfile.blast -t ryansonder/samap-blast:latest .


rebuild:
	docker rmi ryansonder/samap:latest || true
	docker build -t ryansonder/samap:latest .


docker-shell-samap:
	docker run --rm -it \
		-v $(PWD):/workspace \
		-w /workspace \
		--entrypoint /bin/bash \
		ryansonder/samap:latest

docker-shell-blast:
	docker run --rm -it \
		-v $(PWD):/workspace \
		-w /workspace \
		--entrypoint /bin/bash \
		ryansonder/samap-blast:latest


clean-nextflow:
	rm -rf work/*
	rm -rf .nextflow*

clean-results:
	rm -rf results/*

clean-docker:
	docker rmi ryansonder/samap:latest || true

clean-all: clean-nextflow clean-results clean-docker