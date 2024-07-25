set -e

cmake -S ./ -B build -DCMAKE_BUILD_TYPE=Debug -DUSE_MAIN=1 -DUSE_TEST=0
cmake --build build

LD_PRELOAD=/usr/local/lib/libprofiler.so CPUPROFILE=main.prof ./build/main/Main -b Debug -v 0 -r 1 -s 1721650591 -k 10 -m fast_level2 -d random_graph_wTrue_nghb10_n1000.csv

#google-pprof --text ./build/main/Main main.prof
google-pprof --pdf ./build/main/Main main.prof > main.pdf
#google-pprof --lines ./build/main/Main main.prof
