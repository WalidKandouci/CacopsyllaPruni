echo "From migaleTest.sh"
qsub -cwd -q short.q -S /bin/sh -V   /home/dpleydell/work/CacopsyllaPruni/migale/migaleTest2.sh $iModel


## long.q    has a 5 day (120 hours) run time limit
## infinit.q has no limit, but only 3 nodes

## https://migale.inrae.fr/cluster
## https://migale.inrae.fr/sge
