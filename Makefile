all: docopt
	mkdir -p build
	cd build && cmake ../src && make

docopt:
	cat src/bmrfcpu/bmrfcpu.docopt | python src/bmrfcpu/docopt-c/docopt_c.py > src/bmrfcpu/docopt.c
	
clean:
	rm -rf build; rm -f src/bmrfcpu/docopt.c; rm -rf src/bmrfcpu/docopt-c/docopt.pyc

# Simple makefile showing the vagrant commands needed
setup:
	# build box
	cd devops/ubuntu-ec2; make
	# make sure vagrant aws plugin is installed
	vagrant plugin install vagrant-aws
	# add box to vagrant
	vagrant box add ubuntu devops/ubuntu.box --provider=aws

up:
	vagrant up --provider=aws

down:
	vagrant destroy
