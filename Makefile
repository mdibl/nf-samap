run:
	nextflow run main.nf 

docker:
	docker build -f Dockerfile.samap -t docker.io/mdiblbiocore/samap:latest .
	docker build --platform=linux/amd64 -f Dockerfile.blast -t docker.io/mdiblbiocore/samap-blast:latest .
	docker push docker.io/mdiblbiocore/samap:latest
	docker push docker.io/mdiblbiocore/samap-blast:latest  

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