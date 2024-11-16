set -e

#buildtype=Debug
buildtype=RelWithDebInfo

version=601 # @20241116
version=602 # @20241116: k=100
version=603 # @20241116: k=1000
rep=1
seed=$(date +%s)

cmake -S ./ -B build -DCMAKE_BUILD_TYPE=$buildtype -DUSE_MAIN=1 -DUSE_TEST=0
cmake --build build

n=1000

for nbr in 10 100 1000 ; do
    for method in fast_level2 ; do

        #./build/main/Main -b $buildtype -v $version -r $rep -s $seed -d random_graph_wTrue_nghb${nbr}_n${n}.csv -k 10 -t "note" -m $method
        ./build/main/Main -b $buildtype -v $version -r $rep -s $seed -d random_graph_wTrue_nghb${nbr}_n${n}.csv -k 1000 -t "note" -m $method
        :

    done
done


# remember to modify DSTAlgoEvaluation/readinput.py to set method=fast_level2
cd DSTAlgoEvaluation/src
javac -cp .:gson-2.11.0.jar Main.java

while read -r line
do
    echo "$line"
    java -cp .:gson-2.11.0.jar Main $line
done < <(python3 ../readinput.py $version $rep )

