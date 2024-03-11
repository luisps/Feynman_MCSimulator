#!/bin/bash
appname="Feynman"
jobname="F""_$1_$2_S$6"
filename="$appname""_$1_$2_$3_$4_T$5_S$6"
RL_label=""
if  [ $9 ]; then
	filename+="_L$9"
	RL_label+=" -l $9"
fi
if  [ $7 ] && [ $8 ]; then
	RL_label+=" -r $7 $8"
fi
sbatch <<EOT
#!/bin/bash

#SBATCH --job-name=$jobname
#SBATCH --time=1:40:0
#SBATCH --partition=fct
#SBATCH --qos=cpca095372023
#SBATCH --output=$filename.o
#SBATCH --error=$filename.e
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=$5

./Feynman -f ~/QCircuits_BenchTest/circuits/$1/circuit_$1 -a $2 -i $3 -o $4 -t $5 -s $6 $RL_label

if [ -s "$filename".e ]; then
	mv "$filename".e ~/QCircuits_BenchTest/circuits/$1
else
	rm "$filename".e 
fi

if [ -s "$filename".o ]; then
	mv "$filename".o ~/QCircuits_BenchTest/circuits/$1

fi

exit 0
EOT
