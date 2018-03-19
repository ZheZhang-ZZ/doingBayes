BernGrid<-function(Theta, pTheta, Data, credib=.95, nToPlot=length(Theta)){
    z <- sum(Data==1)
    N <- length(Data)
    # Compute likelihood
    pDataGivenTheta <- Theta^z * (1-Theta)^(N-z)
    # Compute the evidence and the posterior
    pData = sum(pDataGivenTheta * pTheta)
    pThetaGivenData = pDataGivenTheta * pTheta / pData
    
    # Plot the results
    layout( matrix( c( 1,2,3 ) ,nrow=3 ,ncol=1 ,byrow=FALSE ) ) # 3x1 panels
    par( mar=c(3,3,1,0) , mgp=c(2,1,0) , mai=c(0.5,0.5,0.3,0.1) ) # margin settings
    dotsize = 5 # how big to make the plotted dots
    # If the comb has a zillion teeth, it's too many to plot, so plot only a
    # thinned out subset of the teeth.
    nteeth = length(Theta)
    if ( nteeth > nToPlot ) {
      thinIdx = seq( 1, nteeth , round( nteeth / nToPlot ) )
      if ( length(thinIdx) < length(Theta) ) {
        thinIdx = c( thinIdx , nteeth ) # makes sure last tooth is included
      }
    } else { thinIdx = 1:nteeth }
    # Plot the prior.
    meanTheta = sum( Theta * pTheta ) # mean of prior, for plotting
    plot( Theta[thinIdx] , pTheta[thinIdx] , type="p" , pch="." , cex=dotsize ,
          xlim=c(0,1) , ylim=c(0,1.1*max(pThetaGivenData)) , cex.axis=1.2 ,
          xlab=bquote(theta) , ylab=bquote(p(theta)) , cex.lab=1.5 ,
          main="Prior" , cex.main=1.5 , col="skyblue" )
    if ( meanTheta > .5 ) {
       textx = 0 ; textadj = c(0,1)
    } else {
      textx = 1 ; textadj = c(1,1)
    }
    text( textx , 1.0*max(pThetaGivenData) ,
          bquote( "mean(" * theta * ")=" * .(signif(meanTheta,3)) ) ,
          cex=2.0 , adj=textadj )
    # Plot the likelihood: p(Data|Theta)
    plot(Theta[thinIdx] ,pDataGivenTheta[thinIdx] ,type="p" ,pch="." ,cex=dotsize
        ,xlim=c(0,1) ,cex.axis=1.2 ,xlab=bquote(theta) 
        ,ylim=c(0,1.1*max(pDataGivenTheta)) 
        ,ylab=bquote( "p(D|" * theta * ")" )  
        ,cex.lab=1.5 ,main="Likelihood" ,cex.main=1.5 , col="skyblue" )
    if ( z > .5*N ) { textx = 0 ; textadj = c(0,1) }
    else { textx = 1 ; textadj = c(1,1) }
    text( textx ,1.0*max(pDataGivenTheta) ,cex=2.0
        ,bquote( "Data: z=" * .(z) * ",N=" * .(N) ) ,adj=textadj )
    # Plot the posterior.
    meanThetaGivenData = sum( Theta * pThetaGivenData )
    plot(Theta[thinIdx] ,pThetaGivenData[thinIdx] ,type="p" ,pch="." ,cex=dotsize
        ,xlim=c(0,1) ,ylim=c(0,1.1*max(pThetaGivenData)) ,cex.axis=1.2 
        ,xlab=bquote(theta) ,ylab=bquote( "p(" * theta * "|D)" )
        ,cex.lab=1.5 ,main="Posterior" ,cex.main=1.5 , col="skyblue" )
    if ( meanThetaGivenData > .5 ) { textx = 0 ; textadj = c(0,1) } 
    else { textx = 1 ; textadj = c(1,1) }
    text(textx ,1.00*max(pThetaGivenData) ,cex=2.0
        ,bquote( "mean(" * theta * "|D)=" * .(signif(meanThetaGivenData,3)) ) 
        ,adj=textadj )
    text(textx ,0.75*max(pThetaGivenData) ,cex=2.0
        ,bquote( "p(D)=" * .(signif(pData,3)) ) ,adj=textadj )
    # Mark the highest density interval. HDI points are not thinned in the plot.
    source("HDIofGrid.R")
    HDIinfo = HDIofGrid( pThetaGivenData , credMass=credib )
    points( Theta[ HDIinfo$indices ] ,
           rep( HDIinfo$height , length( HDIinfo$indices ) ) , pch="-" , cex=1.0 )
    text( mean( Theta[ HDIinfo$indices ] ) , HDIinfo$height ,
             bquote( .(100*signif(HDIinfo$mass,3)) * "% HDI" ) ,
             adj=c(0.5,-1.5) , cex=1.5 )
    # Mark the left and right ends of the waterline. 
    # Find indices at ends of sub-intervals:
    inLim = HDIinfo$indices[1] # first point
    for ( idx in 2:(length(HDIinfo$indices)-1) ) {
      if ( ( HDIinfo$indices[idx] != HDIinfo$indices[idx-1]+1 ) | # jumps on left, OR
        ( HDIinfo$indices[idx] != HDIinfo$indices[idx+1]-1 ) ) { # jumps on right
        inLim = c(inLim,HDIinfo$indices[idx]) # include idx
      }
    }
    inLim = c(inLim,HDIinfo$indices[length(HDIinfo$indices)]) # last point
    # Mark vertical lines at ends of sub-intervals:
    for ( idx in inLim ) {
      lines( c(Theta[idx],Theta[idx]) , c(-0.5,HDIinfo$height) , type="l" , lty=2 , 
             lwd=1.5 )
      text( Theta[idx] , HDIinfo$height , bquote(.(round(Theta[idx],3))) ,
            adj=c(0.5,-0.1) , cex=1.2 )
    }

    return( pThetaGivenData )
    } # end of function

setwd("C:\\Zhangzhe\\2-learning\\统计基因组与生物信息学\\doin bayes")

# Create vector of theta values
binwidth = 1/1000
thetagrid = seq(from=binwidth/2, to=1-binwidth/2, by=binwidth)

# Specify probability mass at each theta value
relprob = pmin(thetagrid,1-thetagrid)
prior = relprob / sum(relprob)
datavec = c( rep(1,3) , rep(0,1) )
windows(7,10)
posterior = BernGrid( Theta=thetagrid , pTheta=prior , Data=datavec )


### exercise 6.2
pTheta = c( 50:1 , rep(1,50) , 1:50 , 50:1,rep(1,50) , 1:50)
pTheta = pTheta / sum( pTheta )
width = 1 / length(pTheta)
Theta = seq( from = width/2 , to = 1-width/2 , by = width )
Data<-c(rep(1,15),rep(0,5))
windows(7,10)
posterior = BernGrid( Theta=Theta , pTheta=pTheta , Data=Data )

### exercise 6.3
pTheta = c( 50:1 , rep(1,50) , 1:50 , 50:1,rep(1,50) , 1:50)
pTheta = pTheta / sum( pTheta )
width = 1 / length(pTheta)
Theta = seq( from = width/2 , to = 1-width/2 , by = width )
Data<-c(rep(1,3),rep(0,1))
windows(7,10)
posterior = BernGrid( Theta=Theta , pTheta=pTheta , Data=Data )
Data<-c(rep(1,12),rep(0,4))
windows(7,10)
posterior = BernGrid( Theta=Theta , pTheta=posterior , Data=Data )

### exercise 6.4
pTheta = rep(1,100)
pTheta = pTheta / sum( pTheta )
width = 1 / length(pTheta)
Theta = seq( from = width/2 , to = 1-width/2 , by = width )
Data<-c(rep(1,58),rep(0,42))
windows(7,10)
posterior = BernGrid( Theta=Theta , pTheta=pTheta , Data=Data )
Data<-c(rep(1,57),rep(0,43))
windows(7,10)
posterior = BernGrid( Theta=Theta , pTheta=posterior , Data=Data )

## exercise 6.6
binwidth = 1/1000
thetagrid = seq(from=binwidth/2, to=1-binwidth/2, by=binwidth)

relprob = thetagrid^2
prior = relprob / sum(relprob)
datavec = c( rep(1,2) , rep(0,2) )
windows(7,10)
posterior = BernGrid( Theta=thetagrid , pTheta=prior , Data=datavec )

## exercise 7.3 B
binwidth = 1/1000
thetagrid = seq(from=binwidth/2, to=1-binwidth/2, by=binwidth)
pTheta = (cos(4*pi*thetagrid)+1)^2
pTheta = pTheta / sum( pTheta )
Data<-c(rep(1,8),rep(0,4))
windows(7,10)
posterior = BernGrid( Theta=thetagrid , pTheta=pTheta , Data=Data )
