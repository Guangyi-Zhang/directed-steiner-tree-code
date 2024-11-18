set -e

#buildtype=Debug
buildtype=RelWithDebInfo

version=301 # @20240722:
version=302 # @20240724: flac
version=303 # @20240725:
version=304 # @20240725:
version=305 # @20240725:
version=306 # @20240725:
rep=1
seed=$(date +%s)

cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0
cmake --build build

n=10000


seed=1721906676
for k in 10 100; do
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m level3 -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
done

exit

seed=1721906769
for k in 10 100; do
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m level3 -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
done

exit






for k in 10 100 500 1000; do
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m level1 -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
    sleep 3 # takes time to write many terminals into log
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m adaptive_level1 -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
    sleep 3
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level2 -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
    sleep 3
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m level2 -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
done

for k in 10 100; do
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a 0.5 -n 50 -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
done

