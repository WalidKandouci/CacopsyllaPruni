#!/bin/sh

## The following takes the qsub job number and prints it to the qsub output file in ~/
## Note. This works in BASH shells but not C shells.
## By default migale uses a C shell
## So use a -S /bin/sh argument when using qsub

qsubID=$JOB_ID
echo "qsub process number is" $qsubID

## Now call the R script
echo "Now source the R script..."
echo 'source("/home/dpleydell/work/CacopsyllaPruni/migale/test.R");' | R --vanilla --quiet --args 
