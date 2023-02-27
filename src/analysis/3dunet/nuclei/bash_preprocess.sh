#! /bin/bash

PART="075"
matlab -nodesktop -nodisplay -nosplash -r "part='$PART'; run preprocess_3dunet(part).m; exit;"