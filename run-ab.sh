set -e

#buildtype=Debug
buildtype=RelWithDebInfo

version=501 # @20240725:
version=502 # @20240725:
version=503 # @20240725: 10000 too big
version=504 # @20240725:
version=505 # @20240725:
version=506 # @20240725: should only turn one off
version=507 # @20240727: alpha=0.5
version=508 # @20240727: alpha=0.5 n=50
rep=1
seed=$(date +%s)

n=10000
alpha=0.5

#for flag in 0 1; do
#    cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0 -DNO_FAST2_LB=$flag -DNO_FAST3_LB=0 -DNO_FAST3_PRUNE_U=0 -DNO_FAST3_PRUNE_SUBTREE=0
#    cmake --build build
#    cmake -S ./ -B build -U NO_*
#
#    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level2 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
#done

n=1000
for flag in 0 1; do
    cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0 -DNO_FAST2_LB=0 -DNO_FAST3_LB=$flag -DNO_FAST3_PRUNE_U=1 -DNO_FAST3_PRUNE_SUBTREE=1
    cmake --build build
    cmake -S ./ -B build -U NO_*

    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
done

flag=0
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0 -DNO_FAST2_LB=0 -DNO_FAST3_LB=1 -DNO_FAST3_PRUNE_U=$flag -DNO_FAST3_PRUNE_SUBTREE=1
cmake --build build
cmake -S ./ -B build -U NO_*

nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &

flag=0
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0 -DNO_FAST2_LB=0 -DNO_FAST3_LB=1 -DNO_FAST3_PRUNE_U=1 -DNO_FAST3_PRUNE_SUBTREE=$flag
cmake --build build
cmake -S ./ -B build -U NO_*

nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
