run:
	nextflow run main.nf --with-docker


docker:
	docker build -t ryansonder/samap:latest .


rebuild:
	docker rmi ryansonder/samap:latest || true
	docker build -t ryansonder/samap:latest .


docker-shell:
	docker run --rm -it \
		-v $(PWD):/workspace \
		-w /workspace \
		--entrypoint /bin/bash \
		ryansonder/samap:latest


clean-nextflow:
	rm -rf work/*
	rm -rf .nextflow*

clean-results:
	rm -rf results/*

clean-docker:
	docker rmi ryansonder/samap:latest || true

clean-all: clean-nextflow clean-results clean-docker