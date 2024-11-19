set -e

#buildtype=Debug
buildtype=RelWithDebInfo

version=401 # @20240722:
version=402 # @20240722: bug in inconsecutive vertex id
version=403 # @20240722: random root
version=404 # @20240722: drop SFRoad
version=405 # @20240723: bug in fast_level2 denlb
version=406 # @20240723: bug in fast_level2 <
version=407 # @20240724: flac
version=408 # @20240725: flac reachable terms
version=409 # @20240727: n=50
version=410 # @20241117: road-CA
version=411 # @20241117: token_transfers.csv
version=413 # @20241117: token_transfers.csv
version=414 # @20241118: test road-CA
version=415 # @20241118: test wiki-topcats.txt
version=416 # @20241118: reserve wiki-topcats.txt
version=417 # @20241119: reserve soc data
version=418 # @20241119: reserve soc data
seed=$(date +%s)

cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0
cmake --build build

k=10

# too dis-connected
#for dataset in cit-HepPh.txt; do
#for dataset in SFRoad; do
#for dataset in Hongkong.road-d; do
#for dataset in rec-libimseti-dir.edges; do

# weird transaction network
#for dataset in token_transfers.csv; do

#for dataset in soc-Epinions1.txt web-Google.txt token_transfers.csv soc-pokec-relationships.txt advogato.edges soc-LiveJournal1.txt; do

rep=6
#for dataset in advogato.edges; do
#for dataset in soc-Epinions1.txt; do
#for dataset in web-Google.txt; do
#for dataset in soc-pokec-relationships.txt; do
#for dataset in soc-LiveJournal1.txt; do
for dataset in wiki-topcats.txt; do
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m level2 -d $dataset -k $k -t "note" &
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m level1 -d $dataset -k $k -t "note" &
    nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m adaptive_level1 -d $dataset -k $k -t "note" &
    #nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level2 -d $dataset -k $k -t "note" &

    #if [ "$dataset" != "soc-LiveJournal1.txt" ] && [ "$dataset" != "web-Google.txt" ]; then
    if [ "$dataset" != "soc-LiveJournal1.txt" ]; then
        :
        #nohup ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level3 -a 0.5 -n 50 -d $dataset -k $k -t "note" &
    fi

    ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -m fast_level2 -d $dataset -k $k -t "note" > tmp

    cd DSTAlgoEvaluation/src
    javac -cp .:gson-2.11.0.jar Main.java

    while read -r line
    do
        echo "$line"
        java -cp .:gson-2.11.0.jar Main $line
    done < <(python3 ../readinput.py $version $rep $dataset )
done


