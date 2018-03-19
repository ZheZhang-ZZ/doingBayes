# Figure 4.2
nThetaVals <- 63
Theta <- seq(from=1/(nThetaVals+1),to=nThetaVals/(nThetaVals+1),by=1/(nThetaVals+1))
pTheta <- pmin(Theta,1-Theta)
pTheta <- pTheta / sum(pTheta)
Data = c(1,1,1,0,0,0,0,0,0,0,0,0)
nHeads = sum(Data==1)
nTails = sum(Data==0)

#likelihood
pDataGivenTheta = Theta^nHeads*(1-Theta)^nTails

#posterior
pData = sum(pDataGivenTheta * pTheta)
pThetaGivenData = pDataGivenTheta * pTheta / pData

#plot
windows(7,10)
layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=F))
par(mar=c(3,3,1,0))
par(mgp=c(2,1,0))
par(mai=c(0.5,0.5,0.3,0.1))

##prior plot
plot(Theta, pTheta, type="h", lwd=3, main="Prior", xlim=c(0,1), xlab=bquote(theta), ylim=c(0,1.1*max(pThetaGivenData)), ylab=bquote(p(theta)), cex.axis=1.2, cex.lab=1.5, cex.main=1.5)

##likelihood plot
plot(Theta, pDataGivenTheta, type="h", lwd=3, main="Likelihood", xlim=c(0,1), xlab=bquote(theta), ylim=c(0,1.1*max(pDataGivenTheta)), ylab=bquote(paste("p(D|)", theta, ")")), cex.axis=1.2, cex.lab=1.5, cex.main=1.5)
text(0.55,0.85*max(pDataGivenTheta), cex=2.0, bquote("D=" * .(nHeads) * "H, " * .(nTails) * "T"), adj=c(0,0.5))

##posteior plot
plot(Theta, pThetaGivenData, type="h", lwd=3, main="Posterior", xlim=c(0,1), xlab=bquote(theta), ylim=c(0,1.1*max(pThetaGivenData)), ylab=bquote(paste("p(",theta,"|D)")), cex.axis=1.2, cex.lab=1.5, cex.main=1.5)
text(0.55,0.85*max(pThetaGivenData), cex=2.0, bquote("D=" * .(signif(pData,3))), adj=c(0,0.5))
