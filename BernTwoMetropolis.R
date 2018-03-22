library(MASS)

likelihood = function(theta){
    z1 = 5; N1 = 7; z2 = 2; N2 = 7
    likelihood = (theta[1]^z1*(1-theta[1])^(N1-z1)
                 * theta[2]^z2*(1-theta[2])^(N2-z2))
    return(likelihood)
    }

prior = function(theta){
    a1 = 3; b1 = 3; a2 = 3; b2 = 3
    prior = dbeta(theta[1], a1, b1) * dbeta(theta[2], a2, b2)
    return(prior)
    }

targetRelProb = function(theta){
    if(all(theta>=0) & all(theta<=1.0)){
        targetRelProbVal = likelihood(theta) * prior(theta)
        }else{
          targetRelProbVal = 0.0
        }
    return(targetRelProbVal)
    }

trajLength = ceiling(1000 / .9)
trajectory = matrix(0, nrow=trajLength, ncol=2)
trajectory[1,] = c(0.5, 0.5)
burnIn = ceiling(.1 * trajLength)
nAccepted = 0
nrejected= 0

set.seed(47405)
nDim = 2; sd1 = 0.2; sd2 = 0.2;
covarMat = matrix(c(sd1^2,0,0,sd2^2), nrow=nDim, ncol=nDim)
for(stepIdx in 1:(trajLength-1)){
  currentPosition = trajectory[stepIdx,]
  proposedJump = mvrnorm(n=1, mu=rep(0,nDim), Sigma=covarMat)
  probAccept = min(1, targetRelProb(currentPosition + proposedJump)
                   /targetRelProb(currentPosition))
  if(runif(1) < probAccept){
    trajectory[stepIdx+1, ] = currentPosition + proposedJump
    if(stepIdx > burnIn){nAccepted = nAccepted + 1}
    }else{
    trajectory[stepIdx+1, ] = currentPosition
    if(stepIdx > burnIn){nrejected = nrejected + 1}
    }
}

acceptedTraj = trajectory[(burnIn+1):dim(trajectory)[1], ]
meanTraj = apply(acceptedTraj, 2, mean)
sdTraj = apply(acceptedTraj, 2, sd)

par( pty="s" ) # makes plots in square axes.
plot( acceptedTraj , type = "o" , xlim = c(0,1) , xlab = bquote(theta[1]) ,
      ylim = c(0,1) , ylab = bquote(theta[2]) , col="skyblue" )
# Display means and rejected/accepted ratio in plot.
if ( meanTraj[1] > .5 ) { xpos = 0.0 ; xadj = 0.0
} else { xpos = 1.0 ; xadj = 1.0 }
if ( meanTraj[2] > .5 ) { ypos = 0.0 ; yadj = 0.0
} else { ypos = 1.0 ; yadj = 1.0 }
text( xpos , ypos ,	bquote(
  "M=" * .(signif(meanTraj[1],3)) * "," * .(signif(meanTraj[2],3))
  * "; " * N[pro] * "=" * .(dim(acceptedTraj)[1])
  * ", " * frac(N[acc],N[pro]) * "=" 
  * .(signif(nAccepted/dim(acceptedTraj)[1],3))
) , adj=c(xadj,yadj) , cex=1.5  )


