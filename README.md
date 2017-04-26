# USAGE INSTRUCTIONS

## DOWNLOAD

Either download an Elegent release (http://www.hepforge.org/downloads/elegent) or checkout the code from git:
```
	git clone https://github.com/jan-kaspar/elegent.git
```


## BUILD

The build is based on the CMake system (http://www.cmake.org/). While you may use any functionality it supports, the recommended build is as follows:
```
mkdir build
cd build
cmake .. <options>
make
```

CMake will automatically search for GSL, HepMC and ROOT installations. If it fails or if you wish to use a given installation, apply the following options:
 * `-DHEPMC_PREFIX:STRING=<your chosen directory>`
 * `-DROOT_PREFIX:STRING=<your chosen directory>`
 * `-DGSL_CONFIG:STRING=<path to gsl-config>`

The program files (executables, library and header files) can be installed by
```
make install
```

The installation destination can be set by option
```
-DCMAKE_INSTALL_PREFIX=<your chosen directory>
```
passed to cmake at configure time.


## BUILD TEST

To verify the success of the recommended build, you can do the following:
```
cd ../test/build
./package_test ../../build
```


## PROGRAM EXECUTION

All executable programs can be run with `-h` option in oder to get help. Alternatively, you can consult the Wiki pages: http://elegent.hepforge.org/trac/wiki .
