set -e

version=408
rep=1

cd DSTAlgoEvaluation/src
javac -cp .:gson-2.11.0.jar Main.java

line="407 1 soc-Epinions1.txt 1721803096 75879 508837 10 74411 71449,20661,18045,34514,58357,68162,44340,38050,20044,65175"
#echo $line
#java -cp .:gson-2.11.0.jar Main $line
#exit

while read -r line
do
    echo "$line"
    java -cp .:gson-2.11.0.jar Main $line
#done < <(python3 ../readinput.py $version $rep | grep "Live" )
done < <(python3 ../readinput.py $version $rep )
#done < "$input"

