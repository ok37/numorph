#!/bin/bash

SAMPLES="NF1R2M1R"	#or specify multiple like "TEST1 TEST2..."
STAGE="process"

for sample in $SAMPLES
do 
  matlab -nodesktop -nodisplay -nosplash -r "try NM_config('$STAGE','$sample',true); catch fprintf('Error running sample %s\n','$sample'); end; exit;"
done 

echo "Completed All Samples"
