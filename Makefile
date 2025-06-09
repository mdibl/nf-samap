docker:
	docker build -t ryansonder:latest .

clean:
	rm -rf results/
	rm -rf .vscode/nextflow.config
	rm -rf .vscode/
	rm -rf .nextflow/
	rm -rf work/
	rm -f .nextflow.log*
	rm -f nextflow.log
	rm -f nextflow.trace