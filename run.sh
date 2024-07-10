set -e

#buildtype=Debug
#buildtype=RelWithDebInfo

cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_MAIN=1 -DUSE_TEST=0
cmake --build build

./build/main/Main -b RelWithDebInfo -m level3 -a 0.5 -d random_graph_wFalse_p01_1000.csv
./build/main/Main -b RelWithDebInfo -m level3 -a 0.5 -d random_graph_wTrue_p01_1000.csv
./build/main/Main -b RelWithDebInfo -m level3 -a 0.5 -d random_graph_wFalse_p002_1000.csv
./build/main/Main -b RelWithDebInfo -m level3 -a 0.5 -d random_graph_wTrue_p002_1000.csv
#./build/main/Main -b RelWithDebInfo -m level3 -a 0.5 -d soc-Epinions1.txt
