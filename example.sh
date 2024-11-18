cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_TEST=0 -DUSE_MAIN=0 -DUSE_EXAMPLE=1
cmake --build build
./build/example/Example
