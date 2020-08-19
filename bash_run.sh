#! /bin/bash

SAMPLE1="WT8R"
#SAMPLE1="WT1L"
#SAMPLE3="TOP16R"
#SAMPLE1="TOP11L"


if [ -n "${SAMPLE1}" ]; then
  matlab -nodesktop -nodisplay -nosplash -r "run TC_template_"$SAMPLE1".m; exit;"
  #matlab -nodesktop -nodisplay -nosplash -r "run TC_template_"$SAMPLE1"_2.m; exit;"
  #matlab -nodesktop -nodisplay -nosplash -r "run TC_template_"$SAMPLE3".m; exit;"
  #matlab -nodesktop -nodisplay -nosplash -r "run TC_template_"$SAMPLE4".m; exit;"
else
  echo "First parameter not supplied."
fi

echo "Completed Sample "$SAMPLE1""
