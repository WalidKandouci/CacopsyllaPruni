for iModel in $(seq 1 1 6) ## 50
do
    echo "From bashScript1.sh"
    echo $iModel
    qsub -cwd -q long.q -S /bin/sh -V /projet/extern/work/dpleydell/CacopsyllaPruni/migaleScript2.sh $iModel
done

## long.q    has a 5 day (120 hours) run time limit
## infinit.q has no limit, but only 3 nodes

## https://migale.inrae.fr/cluster
## https://migale.inrae.fr/sge
