for iModel in $(seq 1 1 5) ## 50
do
    echo "From bashScript1.sh"
    echo $iModel
    qsub -q long.q -S /bin/sh -V /projet/extern/save/wkandouci/CacopsyllaPruni/bashScript2.sh $iModel
done

## This first script
