build: Dockerfile
	docker build -t hhiiro/photoemission1d:1.0 .

build_clean: Dockerfile
	docker build -t hhiiro/photoemission1d:1.0 --no-cache .

run:
	docker run -it --rm -v D:\Win10\Documents\Docker:/work/output hhiiro/photoemission1d:1.0 

br:
	make build && make run