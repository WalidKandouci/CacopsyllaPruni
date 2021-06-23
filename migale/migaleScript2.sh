#!/bin/sh

## The following takes the qsub job number and prints it to the qsub output file in ~/
## Note. This works in BASH shells but not C shells.
## By default migale uses a C shell
## So use a -S /bin/sh argument when using qsub

qsubID=$JOB_ID
iModel=$1

echo "qsub process number is" $qsubID
echo "iModel is" $iModel

## Now call the R script
echo "Now source the R script..."
echo 'source("/projet/extern/save/wkandouci/CacopsyllaPruni/src/fitModel.R");' | R --vanilla --quiet --args $iModel $qsubID
