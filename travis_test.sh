mkdir build
cd build
cmake ..
make

cd ../test/build
./package_test ../../build
