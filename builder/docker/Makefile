


all: basis-uedge



basis-uedge:
	docker buildx build -f ./Dockerfile-basisuedge -t fusionbillm/basisuedge .
basis-uedge-push:
	docker push fusionbillm/basisuedgde
clean:
	docker system prune
realclean:
	docker system prune -a
