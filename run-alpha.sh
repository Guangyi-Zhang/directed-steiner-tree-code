set -e

#buildtype=Debug
buildtype=RelWithDebInfo

version=101 # @20240721: trying alpha
version=102 # @20240722:
version=103 # @20240722: common seed
version=104 # @20240722: bug in preprocessing
version=105 # @20240725: +nthresholds
version=106 # @20240725: n=1000
rep=1
seed=$(date +%s)

cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0
cmake --build build

n=1000
for alpha in 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1; do
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a $alpha -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
done

alpha=1
for nthresholds in 1 50 100; do # 10 included above
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a $alpha -n $nthresholds -d random_graph_wTrue_nghb10_n${n}.csv -k 10 -t "note" &
done
