#!/bin/bash

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p steinlab
#SBATCH --mem=50g
#SBATCH -t 01-00:00:00

#module add matlab/2019a
module add anaconda
module add gcc

SAMPLES="NF1R2M9R"
STAGE="register"

for sample in $SAMPLES
do
  matlab -nodesktop -nodisplay -nosplash -singleCompThread -r "myCluster = parcluster; myCluster.NumWorkers = 4; parpool(myCluster); NM_config('$STAGE','$sample',true);"
done
