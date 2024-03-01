#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH --job-name=Feynman
#SBATCH --time=0:10:0
#SBATCH --partition=fct
#SBATCH --qos=cpca095372023
#SBATCH --output=%x_$1_$2_$3_$4_T$5_S$6.o
#SBATCH --error=%x_$1_$2_$3_$4_T$5_S$6.e
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=$5

./Feynman ~/QCircuits_BenchTest/circuits/$1/circuit_$1 $2 $3 $4 $5 $6 $7 $8

mv Feynman_$1_$2_$3_$4_T$5_S$6.? ~/QCircuits_BenchTest/circuits/$1

exit 0
EOT
