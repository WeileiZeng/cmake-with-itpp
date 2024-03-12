all:build run
build:
	cd main/build && cmake .. && make
#	cd main/build && cmake .. && make && make test 
run:
	make build
	./main/build/bin/MainProject
#	./bin/examples
#	cd main/build && cmake .. && make && make test && cd ../.. && ./main/build/bin/test_lib
srun:
	make build
	srun -p debug -n 1 -c 32 --time=6:00:00 --cpus-per-task=32 ./main/build/bin/MainProject
clean:
	rm -rf main/build && mkdir main/build
