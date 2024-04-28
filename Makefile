build: Dockerfile
	docker build -t hhiiro/photoemission1d:1.0 .

build_clean: Dockerfile
	docker build -t hhiiro/photoemission1d:1.0 --no-cache .

run:
	docker run -it --rm hhiiro/photoemission1d:1.0 