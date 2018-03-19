setwd("C:\\Zhangzhe\\2-learning\\统计基因组与生物信息学\\doin bayes")
myData = c(rep(1,11),rep(0,3))

# Define the Bernoulli likelihood function
likelihood <- function(theta, data){
    z = sum(data == 1)
    N = length(data)
    pDataGivenTheta = theta^z*(1-theta)^(N-z)
    pDataGivenTheta[theta > 1 | theta < 0] = 0
    return(pDataGivenTheta)
    }

# Define the prior
# prior <- function(theta){
    # prior <- rep(1, length(theta)) # uniform density over [0,1]
    # #prior = dbeta()pmin(2*theta, 2*(1-theta), 2, 2)
    # prior[theta>1 | theta<0] = 0
    # return(prior)
    # }

## exercise 7.3 C
prior <- function(theta){
    prior <- (cos(4*pi*theta)+1)^2
    #prior = dbeta()pmin(2*theta, 2*(1-theta), 2, 2)
    prior[theta>1 | theta<0] = 0
    return(prior)
    }

# Define the relative probability of the target distribution
targetRelProb <- function(theta, data){
    targetRelProb = likelihood(theta, data) * prior(theta)
    return(targetRelProb)
    }

# Specify the length of trajectory
trajLength = 10000
trajectory = rep(0, trajLength)
trajectory[1] = 0.5
burnIn = ceiling(.01 * trajLength)
nAccepted = 0
nRejected = 0
set.seed(47405)

for(t in 1:(trajLength-1)){
    currentPosition = trajectory[t]
    proposedJump = rnorm(1, mean=0, sd=0.001)
    probAccept = min(1,
                targetRelProb(currentPosition + proposedJump, myData) / targetRelProb(currentPosition, myData)
                )
    if(runif(1) < probAccept){
        trajectory[t+1] = currentPosition + proposedJump
        if(t > burnIn){
            nAccepted = nAccepted + 1
            }
        }else{
        trajectory[t+1] = currentPosition
        if(t > burnIn){
            nRejected = nRejected + 1
            }
        }
    }

acceptedTraj = trajectory[(burnIn+1):length(trajectory)]

source("plotPost.r")
histInfo <- plotPost(acceptedTraj, xlim=c(0,1), breaks = 30)

densMax = max( density( acceptedTraj )$y )
meanTraj = mean( acceptedTraj )
sdTraj = sd( acceptedTraj )
if ( meanTraj > .5 ) {
  xpos = 0.0 ; xadj = 0.0
} else {
  xpos = 1.0 ; xadj = 1.0
}
text( xpos , 0.75*densMax ,
    bquote(N[pro] * "=" * .(length(acceptedTraj)) * "  " *
    frac(N[acc],N[pro]) * "=" * .(signif( nAccepted/length(acceptedTraj) , 3 ))
    ) , adj=c(xadj,0))
    

# Evidence for the model, p(D)
# a = meanTraj * ((meanTraj*(1-meanTraj)/sdTraj^2) - 1)
# b = (1-meanTraj) * ((meanTraj*(1-meanTraj)/sdTraj^2) - 1)

a=10
b=10
wtdEvid <- dbeta(acceptedTraj, a, b) / (
           likelihood(acceptedTraj, myData) * prior(acceptedTraj))
pData = 1 / mean(wtdEvid)

# Display p(D) in the graph
if ( meanTraj > .5 ) { xpos = 0.0 ; xadj = 0.0
} else { xpos = 1.0 ; xadj = 1.0 }
text( xpos , 0.9*densMax , bquote( p(D)==.( signif(pData,3) ) ) ,
      adj=c(xadj,0) , cex=1.5 )

windows()
plot(wtdEvid,type="l")

## exercise 7.3 D
library(BRugs)

modelString = "
model {
    for(i in 1:nFlips){
        y[i] ~ dbern(theta)
        }
    theta ~ dbeta(priorA, priorB)
    priorA<-1
    priorB<-1
    theta ~ (cos(4*3.14*theta)+1)^2
    }
"
writeLines(modelString,con="C:\\Users\\au588655\\Desktop\\qtl\\model.txt")
modelCheck( "C:\\Users\\au588655\\Desktop\\qtl\\model.txt" )
