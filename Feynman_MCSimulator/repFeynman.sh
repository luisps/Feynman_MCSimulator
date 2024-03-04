#!/bin/bash
jobname="Feynman"
filename="$jobname""_$1_$2_$3_$4_T$5_S$6"
R_label=""
if  [ $7 ] && [ $8 ]; then
	R_label+=" -r $7 $8"
fi
for run in {1..3}
do
	filename+="_L${run}"
	L_label+=" -l ${run}"
sbatch <<EOT
#!/bin/bash

#SBATCH --job-name=Feynman
#SBATCH --time=0:10:0
#SBATCH --partition=fct
#SBATCH --qos=cpca095372023
#SBATCH --output=$filename.o
#SBATCH --error=$filename.e
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=$5

./Feynman -f ~/QCircuits_BenchTest/circuits/$1/circuit_$1 -a $2 -i $3 -o $4 -t $5 -s $6 $R_label $L_label

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
done
