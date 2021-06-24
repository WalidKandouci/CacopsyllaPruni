for iModel in $(seq 1 1 6) ## 50
do
    echo "From bashScript1.sh"
    echo $iModel
    # qsub -cwd -N model$iModel -q short.q@@recentNodes -S /bin/sh -V /home/dpleydell/work/CacopsyllaPruni/migale/migaleScript2.sh $iModel
    # qsub -cwd -N model$iModel -q long.q@@recentNodes -S /bin/sh -V /home/dpleydell/work/CacopsyllaPruni/migale/migaleScript2.sh $iModel
    qsub -cwd -N model$iModel -q infinit.q@@recentNodes -S /bin/sh -V /home/dpleydell/work/CacopsyllaPruni/migale/migaleScript2.sh $iModel
done

## long.q    has a 5 day (120 hours) run time limit
## infinit.q has no limit, but only 3 nodes

## https://migale.inrae.fr/cluster
## https://migale.inrae.fr/sge
