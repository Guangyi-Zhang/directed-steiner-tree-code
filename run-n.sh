set -e

#buildtype=Debug
buildtype=RelWithDebInfo

version=101 # @20240721: first try
version=201 # @20240722: common seed
version=202 # @20240722: bug in preprocessing
version=203 # @20240723: bug in fast_level2 denlb
version=204 # @20240724: flac, weighted graphs
version=205 # @20240725: flac reachable terms
version=206 # @20240725: level3 takes a lot memory
rep=1
seed=$(date +%s)

cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0
cmake --build build

for n in 1000 10000 100000 1000000; do
    for method in fast_level2 level2 level1 adaptive_level1; do

        nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" -m $method  &

    done
done
sleep 5

for n in 1000 10000 100000; do
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a 0.5 -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
done
sleep 5

for n in 1000 10000; do
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m level3 -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
done

