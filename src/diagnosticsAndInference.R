#################
## Plot stages ##
#################

# The idea I have is to plot at first each stage on it's own, and maybe have the curv of the tree that we want to see in a difrent color ?
# I'm not sure of what and how to plot because "stage" is a big array and I don't know if we should plot everything ? and if we don't plot everything what should we really plot ?
par(mfrow=c(4,4))
for (ii in 1:nstages){
  for(tt in 1:ntrees){
    plot(states[tt, ii, substage],type='l', col="grey", xlab="time", ylab="dev")
    points(states[tt, ii, substage],type='l', col="grey")
  }
}
