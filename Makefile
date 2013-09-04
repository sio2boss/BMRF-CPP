all:
	mkdir -p build
	cd build && cmake ../src && make

docopt:
	cat src/bmrfcpu/bmrfcpu.docopt | python src/bmrfcpu/docopt-c/docopt_c.py > src/bmrfcpu/docopt.c
	
clean:
	rm -rf build; rm -f src/bmrfcpu/docopt.c

