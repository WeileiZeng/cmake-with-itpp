# cmake-with-itpp
dated: Mar 12, 2024

## Introduction
This is a C++ project build with Cmake. It runs a complete search to generate CSS codes at give size n, where n is the number of physical qubits.


## How to run
see `Makefile`. All necessary cmds is put in `Makefile`. One can simply run `make`. Specifically, `make build` will build the project, `make run` or `make srun` will run the script.

## Package usedic
This project utilize the `itpp` package, which is for Information Technology. It has nice functions with binary calculations. All vector, matrices can be saved in binary format. All matrix algebra and gaussian ellimination are modula 2.

doc page: https://itpp.sourceforge.net/4.3.1/installation.html

recommended installation is to install it as a statice library, then put the libitpp.so file into `lib` folder. The header files were already in the `include` folder.


Steps to Compile and install IT++ library as static library:

```
# download the package and untar it to folder, e.g. $HOME/itpp-4.3.0
cd $HOME/itpp-4.3.0
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/it++4.3.0 -DBLA_STATIC=on -DITPP_SHARED_LIB=off
make && make install

# then find the installed file, and copy it to the folder we want to use.
```
mv $HOME/it++4.3.0/lib/libitpp.so folder/of/your/workspace
# or the following cmd if there are .so files with multiple version names
mv $HOME/it++4.3.0/lib/* folder/of/your/workspace/lib
```

```


# TODO
- [ ] class CSSCode
- [ ] mmio to itfiles
- [ ] test for bvec and GF2mat

