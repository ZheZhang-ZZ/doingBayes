setwd("C:\\Zhangzhe\\2-learning\\statistical genomics and bioinformatics\\doin bayes")
library(BRugs)
modelstring = "
model{
  for(i in 1:N1){y1[i] ~ dbern(theta1)}
  for(i in 1:N2){y2[i] ~ dbern(theta2)}
  theta1 ~ dbeta(3, 3)
  theta2 ~ dbeta(3, 3)
}
"

.temp = file("model.txt", "w")
writeLines(modelstring, con=.temp)
close(.temp)
modelCheck("model.txt")

datalist = list(
  N1 = 7,
  y1 = c(1,1,1,1,1,0,0),
  N2 = 7,
  y2 = c(1,1,0,0,0,0,0)
)

modelData(bugsData(datalist))

modelCompile()
modelGenInits()

samplesSet(c("theta1", "theta2"))
chainlength = 10000
modelUpdate(chainlength)

theta1Sample = samplesSample("theta1")
theta2Sample = samplesSample("theta2")

windows()
par( pty="s" )
plot( theta1Sample[(chainlength-500):chainlength] ,
      theta2Sample[(chainlength-500):chainlength] , type = "o" ,
      xlim = c(0,1) , xlab = bquote(theta[1]) , ylim = c(0,1) ,
      ylab = bquote(theta[2]) , main="BUGS Result" )
# Display means in plot.
theta1mean = mean(theta1Sample)
theta2mean = mean(theta2Sample)
if (theta1mean > .5) { xpos = 0.0 ; xadj = 0.0
} else { xpos = 1.0 ; xadj = 1.0 }
if (theta2mean > .5) { ypos = 0.0 ; yadj = 0.0
} else { ypos = 1.0 ; yadj = 1.0 }
text( xpos , ypos ,
      bquote(
        "M=" * .(signif(theta1mean,3)) * "," * .(signif(theta2mean,3))
      ) ,adj=c(xadj,yadj) ,cex=1.5  )
#dev.copy2eps(file="BernTwoBugs.eps")

# Plot a histogram of the posterior differences of theta values.
thetaDiff = theta1Sample - theta2Sample
source("plotPost.R")
windows(7,4)
plotPost( thetaDiff , xlab=expression(theta[1]-theta[2]) , compVal=0.0 ,
          breaks=30 )
#dev.copy2eps(file="BernTwoBugsDiff.eps")

# For Exercise 8.5:
# Posterior prediction. For each step in the chain, use the posterior thetas 
# to flip the coins.
chainLength = length( theta1Sample )
# Create matrix to hold results of simulated flips:
yPred = matrix( NA , nrow=2 , ncol=chainLength ) 
for ( stepIdx in 1:chainLength ) { # step through the chain
  # flip the first coin:
  pHead1 = theta1Sample[stepIdx]
  yPred[1,stepIdx] = sample( x=c(0,1), prob=c(1-pHead1,pHead1), size=1 )
  # flip the second coin:
  pHead2 = theta2Sample[stepIdx]
  yPred[2,stepIdx] = sample( x=c(0,1), prob=c(1-pHead2,pHead2), size=1 )
}
# Now determine the proportion of times that y1==1 and y2==0
pY1eq1andY2eq0 = sum( yPred[1,]==1 & yPred[2,]==0 ) / chainLength