docker-imgs:
	docker build --platform linux/amd64 -t medbioinf/msfragger -f docker/msfragger/Dockerfile docker/msfragger/.
	docker build --platform linux/amd64 -t medbioinf/oktoberfest:latest -f docker/oktoberfest/Dockerfile docker/oktoberfest
