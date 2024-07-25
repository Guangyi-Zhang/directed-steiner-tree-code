set -e

#buildtype=Debug
buildtype=RelWithDebInfo

version=501 # @20240725:
version=502 # @20240725:
version=503 # @20240725:
rep=1
seed=$(date +%s)

n=10000
alpha=1

for flag in 0 1; do
    cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0 -DNO_FAST2_LB=$flag -DNO_FAST3_LB=0 -DNO_FAST3_PRUNE_U=0 -DNO_FAST3_PRUNE_SUBTREE=0
    cmake --build build
    cmake -S ./ -B build -U NO_*

    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level2 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
done

for flag in 0 1; do
    cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0 -DNO_FAST2_LB=0 -DNO_FAST3_LB=$flag -DNO_FAST3_PRUNE_U=0 -DNO_FAST3_PRUNE_SUBTREE=0
    cmake --build build
    cmake -S ./ -B build -U NO_*

    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
done

flag=1
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0 -DNO_FAST2_LB=0 -DNO_FAST3_LB=0 -DNO_FAST3_PRUNE_U=$flag -DNO_FAST3_PRUNE_SUBTREE=0
cmake --build build
cmake -S ./ -B build -U NO_*

nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &

flag=1
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0 -DNO_FAST2_LB=0 -DNO_FAST3_LB=0 -DNO_FAST3_PRUNE_U=0 -DNO_FAST3_PRUNE_SUBTREE=$flag
cmake --build build
cmake -S ./ -B build -U NO_*

nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
