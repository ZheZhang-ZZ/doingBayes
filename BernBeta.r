setwd("C:\\Zhangzhe\\2-learning\\统计基因组与生物信息学\\doin bayes")

BernBeta <- function(priorShape, dataVec, credMass=0.95, saveGraph=F){
    
    # check for errors in input arguments:
    if(length(priorShape) != 2){
        stop("priorShape must have two components.")
        }
    if(any(priorShape<=0)){
        stop("priorShape must be positive.")
        }
    if(credMass<=0 | credMass>=1.0){
        stop("credMass must be between 0 and 1.")
        }
    
    # Rename the prior shape parameters, for convenience:
    a = priorShape[1]
    b = priorShape[2]
    
    #Create summary values of the data
    z = sum(dataVec==1)
    N = length(dataVec)
    
    # Compute the posterior shape parameters:
    postShape = c(a+z, b+N-z)
    
    #Compute the evidence, p(D):
    pData = beta(z+a, N-z+b)/beta(a, b)
    
    #Determine the limits of the highest density interval
    source("HDIofICDF.r")
    hpdLim = HDIofICDF(qbeta, shape1=postShape[1], shape2=postShape[2])
    
    # Now plot everything:
    # Construct grid of theta values, used for graphing
    binwidth = 0.005
    Theta = seq(from=binwidth/2, to=1-(binwidth/2), by=binwidth)
    
    #Compute the prior at each value of theta
    pTheta = dbeta(Theta, a, b)
    # Compute the likelihood of the data at each value of theta
    pDataGivenTheta = Theta^z*(1-Theta^(N-z))
    # Compute the posterior at each value of theta
    pThetaGivenData = dbeta(Theta, z+a, b+N-z)
    
    # Open a window with three panels
    windows(7,10)
    layout(matrix(c(1,2,3), nrow=3, ncol=1, byrow=F))
    par(mar=c(3,3,1,0), mgp=c(2,1,0), mai=c(0.5,0.5,0.3,0.1))
    maxY = max(c(pTheta, pThetaGivenData))
    
    #Plot the prior
    plot(Theta, pTheta, type="l", lwd=3,xlim=c(0,1), ylim=c(0,maxY), cex.axis=1.2, xlab=bquote(theta), ylab=bquote(p(theta)), cex.lab=1.5,main="Prior", cex.main=1.5)
    if(a>b){textx=0; textadj=c(0,1)}
    else{textx=1; textadj=c(1,1)}
    
    text(textx, 1.0*max(pThetaGivenData), bquote("beta(" * theta * "|" * .(a) * "," * .(b) * ")"), cex=2.0, adj=textadj)
    
    # Plot the likelihood: p(data|theta)
    plot( Theta , pDataGivenTheta , type="l" , lwd=3 ,
      xlim=c(0,1) , cex.axis=1.2 , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pDataGivenTheta)) ,
      ylab=bquote( "p(D|" * theta * ")" ) ,
      cex.lab=1.5 , main="Likelihood" , cex.main=1.5 , col="skyblue" )
      
    if ( z > .5*N ) { textx = 0 ; textadj = c(0,1) }
    else { textx = 1 ; textadj = c(1,1) }
    text( textx , 1.0*max(pDataGivenTheta) , cex=2.0 ,
          bquote( "Data: z=" * .(z) * ",N=" * .(N) ) ,adj=textadj )
    # Plot the posterior.
    plot( Theta , pThetaGivenData  ,type="l" , lwd=3 ,
          xlim=c(0,1) , ylim=c(0,maxY) , cex.axis=1.2 ,
          xlab=bquote(theta) , ylab=bquote( "p(" * theta * "|D)" ) ,
          cex.lab=1.5 , main="Posterior" , cex.main=1.5 , col="skyblue" )
    if ( a+z > b+N-z ) { textx = 0 ; textadj = c(0,1) }
    else { textx = 1 ; textadj = c(1,1) }
    text( textx , 1.00*max(pThetaGivenData) , cex=2.0 ,
          bquote( "beta(" * theta * "|" * .(a+z) * "," * .(b+N-z) * ")"  ) ,
          adj=textadj )
    text( textx , 0.75*max(pThetaGivenData) , cex=2.0 ,
          bquote( "p(D)=" * .(signif(pData,3)) ) , adj=textadj )
    # Mark the HDI in the posterior.
    hpdHt = mean( c( dbeta(hpdLim[1],a+z,b+N-z) , dbeta(hpdLim[2],a+z,b+N-z) ) )
    lines( c(hpdLim[1],hpdLim[1]) , c(-0.5,hpdHt) , type="l" , lty=2 , lwd=1.5 )
    lines( c(hpdLim[2],hpdLim[2]) , c(-0.5,hpdHt) , type="l" , lty=2 , lwd=1.5 )
    lines( hpdLim , c(hpdHt,hpdHt) , type="l" , lwd=2 )
    text( mean(hpdLim) , hpdHt , bquote( .(100*credMass) * "% HDI" ) ,
          adj=c(0.5,-1.0) , cex=2.0 )
    text( hpdLim[1] , hpdHt , bquote(.(round(hpdLim[1],3))) ,
          adj=c(1.1,-0.1) , cex=1.2 )
    text( hpdLim[2] , hpdHt , bquote(.(round(hpdLim[2],3))) ,
          adj=c(-0.1,-0.1) , cex=1.2 )
          
    # Construct filename for saved graph, and save the graph.
    if( saveGraph ) {
      filename = paste( "BernBeta_",a,"_",b,"_",z,"_",N ,sep="")
      saveGraph( file = filename , type="eps" )
        }
    return( postShape )
    }

dataVec<-sample(c(0,1),10,replace=T)
priorShape<-c(1,1)
BernBeta(priorShape, dataVec)

dataVec<-c(rep(1,4),0)
priorShape<-c(0.5,0.5)
BernBeta(priorShape, dataVec)

priorShape1 <- c(100,1)
priorShape2 <- c(1,100)

dataVec<-c(rep(1,8),rep(0,4))

post1<-BernBeta(priorShape1, dataVec)
post2<-BernBeta(priorShape2, dataVec)


