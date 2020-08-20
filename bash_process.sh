#! /bin/bash

SAMPLE1="WT7R"
SAMPLE2="WT11L"
SAMPLE3="WT1L"
SAMPLE4="WT8R"
SAMPLE5="TOP16R"
SAMPLE6="TOP14R"
SAMPLE7="TOP11L"
SAMPLE8="TOP110R"

if [ -n "${SAMPLE1}" ]; then

#  matlab -nodesktop -nodisplay -nosplash -r "sample='$SAMPLE1'; run ./templates/NMp_template.m; exit;"
#  matlab -nodesktop -nodisplay -nosplash -r "sample='$SAMPLE2'; run ./templates/NMp_template.m; exit;"
#  matlab -nodesktop -nodisplay -nosplash -r "sample='$SAMPLE3'; run ./templates/NMp_template.m; exit;"
#  matlab -nodesktop -nodisplay -nosplash -r "sample='$SAMPLE4'; run ./templates/NMp_template.m; exit;"
  matlab -nodesktop -nodisplay -nosplash -r "sample='$SAMPLE5'; run ./templates/NMp_template.m; exit;"
#  matlab -nodesktop -nodisplay -nosplash -r "sample='$SAMPLE6'; run ./templates/NMp_template.m; exit;"
#  matlab -nodesktop -nodisplay -nosplash -r "sample='$SAMPLE7'; run ./templates/NMp_template.m; exit;"
#  matlab -nodesktop -nodisplay -nosplash -r "sample='$SAMPLE8'; run ./templates/NMp_template.m; exit;"

else
  echo "First parameter not supplied."
fi

echo "Completed All Samples"
