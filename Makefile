run:
	nextflow run main.nf 

docker:
	docker build -f Dockerfile.samap -t pipeline/samap:latest .
	docker build -f Dockerfile.blast -t pipeline/samap-blast:latest .

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

clean:
	rm -rf work/*
	rm -rf .nextflow*
	rm -rf results/*