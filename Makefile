all:build
build:
	cd main/build && cmake .. && make && make test && ./bin/examples
