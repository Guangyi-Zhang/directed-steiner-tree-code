set -e

#buildtype=Debug
buildtype=RelWithDebInfo

version=301 # @20240722:
version=302 # @20240724: flac
version=303 # @20240725:
version=304 # @20240725:
version=305 # @20240725:
version=306 # @20240725:
version=307 # @20250201: large k
rep=2
seed=$(date +%s)

cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0
cmake --build build

n=100000


for k in 10 100 1000 10000; do
    for method in fast_level2 level2 level1 adaptive_level1; do
        if [ $k -eq 10000 ] && [ "$method" = "level2" ]; then # use up to 50G memory
            continue
        fi  
        nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m $method -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
        sleep 3 # takes time to write many terminals into log
    done
done

for k in 10 100; do
    #nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a 0.5 -n 50 -d random_graph_wTrue_nghb10_n${n}.csv -k $k -t "note" &
    echo
done

