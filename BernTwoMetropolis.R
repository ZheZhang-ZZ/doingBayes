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
        targetRelProb = 0.0
        }
    return(targetRelProb)
    }

trajLength = ceiling(1000 / .9)
trajectory = matrix(0, nrow=trajLength, ncol=2)
trajectory[1,] = c(0.5, 0.5)
burnIn = ceiling(.1 * trajLength)
nAccepted = 0
nrejected= 0

set.seed(47405)
nDim = 2; sd1 = 0.2; sd2 = 0.2;