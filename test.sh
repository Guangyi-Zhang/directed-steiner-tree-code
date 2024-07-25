set -e

# cmake -S standalone -B build/standalone -DCMAKE_BUILD_TYPE=Release
# cmake -S standalone -B build/standalone -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_CXX_FLAGS=-pg
# cmake -S standalone -B build/standalone -DCMAKE_BUILD_TYPE=Debug
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_TEST=1 -DUSE_MAIN=0
cmake --build build
# ./build/test/Test -tc=crossing
if [ $# -eq 0 ]; then
    ./build/test/Test
else
    ./build/test/Test -tc=$1
fi
