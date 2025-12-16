docker-imgs:
	docker pull ghcr.io/medbioinf/mspepid:latest
	
	docker pull proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.25073-842baef
	docker pull quay.io/medbioinf/tdf2mzml:0.4
	docker pull quay.io/medbioinf/openms:3.4.1
	docker pull quay.io/medbioinf/fdrbench-nightly:146f77

	docker pull quay.io/medbioinf/comet-ms:v2024.01.0
	docker pull quay.io/medbioinf/maxquant:2.6.3.0
	docker pull quay.io/medbioinf/msamanda:3.0.22.071
	docker pull quay.io/medbioinf/msgfplus:v2024.03.26
	docker pull quay.io/medbioinf/mzid-merger:1.4.26
	docker pull quay.io/medbioinf/sage:v0.15.0-beta.1
	docker pull quay.io/medbioinf/xtandem:2017.2.1.4
	
	docker pull ghcr.io/percolator/percolator:branch-3-08

	docker build --platform linux/amd64 -t medbioinf/msfragger -f docker/msfragger/Dockerfile docker/msfragger/.

	docker build --platform linux/amd64 -t medbioinf/oktoberfest:latest -f docker/oktoberfest/Dockerfile docker/oktoberfest
