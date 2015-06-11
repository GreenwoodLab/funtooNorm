#########  
# TO CREATE another function studying optional number of components in the PLS
# the creation of quantiles stays the same
# the creation of the control covariate matrix stays the same
# then these commands are run
ncmp.large <- 5  # set this larger than the default
fit2cv.red <- plsr(cbind(quantilesA.red, quantilesB.red) ~ ctl.covmat, ncomp=ncmp.large, validation="CV")
fit2cv.grn <- plsr(cbind(quantilesA.grn, quantilesB.grn) ~ ctl.covmat, ncomp=ncmp.large, validation="CV")
fit2cv.II <- plsr(cbind(quantilesA.II, quantilesB.II) ~ ctl.covmat, ncomp=ncmp.large, validation="CV")

# can you figure out a good way of visualizing the results?  There is one result per quantile for each value of ncomp


#########
# some other visualizations

# PLOT to look at loadings
par(mfrow = c(3,1), omi=c(0,0,0,0), mar=c(2,2,0,0), mgp=c(1,0,0))
plot(fit2cv.red, "loadings", comps=1:ncmp.large)
plot(fit2cv.grn, "loadings", comps=1:ncmp.large)
plot(fit2cv.II, "loadings", comps=1:ncmp.large)


# PLOT to look at fits with the chosen numbmer of components ncmp
# these commands are in the main function
# needs generalizing to any number of cell types or tissue types
fit2.red <- plsr(cbind(quantilesA.red, quantilesB.red) ~ ctl.covmat, ncomp=3)
fit2.grn <- plsr(cbind(quantilesA.grn, quantilesB.grn) ~ ctl.covmat, ncomp=3)
fit2.II <- plsr(cbind(quantilesA.II, quantilesB.II) ~ ctl.covmat, ncomp=3)

par(mfrow = c(3,2), omi=c(0,0,0,0), mar=c(2,2,0,0), mgp=c(1,0,0))
plot(c(0,nqnt), c(0,1), type='n', xlab='quantiles')
ct1 <- fit2.red$fitted.values[ct.CP=='MONO', ,1]
ct2 <- fit2.red$fitted.values[ct.CP=='TCELL', ,1]
ct3 <- fit2.red$fitted.values[ct.CP=='BCELL', ,1]
lines(1:nqnt, apply(ct1,2, mean, na.rm=TRUE), col=1)
lines(1:nqnt, apply(ct2,2, mean, na.rm=TRUE), col=2)
lines(1:nqnt, apply(ct3,2, mean, na.rm=TRUE), col=3)
mtext(side=3, line=-1, cex=0.8, "I, Red; Comp #1")

plot(c(0,nqnt), c(0,1), type='n', xlab='quantiles')
ct1 <- fit2.red$fitted.values[ct.CP=='MONO', ,2]
ct2 <- fit2.red$fitted.values[ct.CP=='TCELL', ,2]
ct3 <- fit2.red$fitted.values[ct.CP=='BCELL', ,2]
lines(1:nqnt, apply(ct1,2, mean, na.rm=TRUE), col=1)
lines(1:nqnt, apply(ct2,2, mean, na.rm=TRUE), col=2)
lines(1:nqnt, apply(ct3,2, mean, na.rm=TRUE), col=3)
mtext(side=3, line=-1, cex=0.8, "I, Red; Comp #2")

par(mfrow = c(3,2), omi=c(0,0,0,0), mar=c(2,2,0,0), mgp=c(1,0,0))
plot(c(0,nqnt), c(0,1), type='n', xlab='quantiles')
ct1 <- fit2.grn$fitted.values[ct.CP=='MONO', ,1]
ct2 <- fit2.grn$fitted.values[ct.CP=='TCELL', ,1]
ct3 <- fit2.grn$fitted.values[ct.CP=='BCELL', ,1]
lines(1:nqnt, apply(ct1,2, mean, na.rm=TRUE), col=1)
lines(1:nqnt, apply(ct2,2, mean, na.rm=TRUE), col=2)
lines(1:nqnt, apply(ct3,2, mean, na.rm=TRUE), col=3)
mtext(side=3, line=-1, cex=0.8, "I, Green; Comp #1")

plot(c(0,nqnt), c(0,1), type='n', xlab='quantiles')
ct1 <- fit2.grn$fitted.values[ct.CP=='MONO', ,2]
ct2 <- fit2.grn$fitted.values[ct.CP=='TCELL', ,2]
ct3 <- fit2.grn$fitted.values[ct.CP=='BCELL', ,2]
lines(1:nqnt, apply(ct1,2, mean, na.rm=TRUE), col=1)
lines(1:nqnt, apply(ct2,2, mean, na.rm=TRUE), col=2)
lines(1:nqnt, apply(ct3,2, mean, na.rm=TRUE), col=3)
mtext(side=3, line=-1, cex=0.8, "I, Green; Comp #2")

par(mfrow = c(3,2), omi=c(0,0,0,0), mar=c(2,2,0,0), mgp=c(1,0,0))
plot(c(0,nqnt), c(0,1), type='n', xlab='quantiles')
ct1 <- fit2.II$fitted.values[ct.CP=='MONO', ,1]
ct2 <- fit2.II$fitted.values[ct.CP=='TCELL', ,1]
ct3 <- fit2.II$fitted.values[ct.CP=='BCELL', ,1]
lines(1:nqnt, apply(ct1,2, mean, na.rm=TRUE), col=1)
lines(1:nqnt, apply(ct2,2, mean, na.rm=TRUE), col=2)
lines(1:nqnt, apply(ct3,2, mean, na.rm=TRUE), col=3)
mtext(side=3, line=-1, cex=0.8, "II; Comp #1")

plot(c(0,nqnt), c(0,1), type='n', xlab='quantiles')
ct1 <- fit2.II$fitted.values[ct.CP=='MONO', ,2]
ct2 <- fit2.II$fitted.values[ct.CP=='TCELL', ,2]
ct3 <- fit2.II$fitted.values[ct.CP=='BCELL', ,2]
lines(1:nqnt, apply(ct1,2, mean, na.rm=TRUE), col=1)
lines(1:nqnt, apply(ct2,2, mean, na.rm=TRUE), col=2)
lines(1:nqnt, apply(ct3,2, mean, na.rm=TRUE), col=3)
mtext(side=3, line=-1, cex=0.8, "II; Comp #2")

cores=1
pls.options(parallel = cores)
ncmp.large <- 5   # set this larger than the default
fit2cvA.red <- plsr(quantilesA.red ~ ctl.covmat, ncomp=ncmp.large, validation="CV")
fit2cvB.red <- plsr(quantilesB.red ~ ctl.covmat, ncomp=ncmp.large, validation="CV")
fit2cvA.grn <- plsr(quantilesA.grn ~ ctl.covmat, ncomp=ncmp.large, validation="CV")
fit2cvB.grn <- plsr(quantilesB.grn ~ ctl.covmat, ncomp=ncmp.large, validation="CV")
fit2cvA.II <- plsr(quantilesA.II ~ ctl.covmat, ncomp=ncmp.large, validation="CV")
fit2cvB.II <- plsr(quantilesB.II ~ ctl.covmat, ncomp=ncmp.large, validation="CV")



op <- par(
    oma=c(3,0,3,0),# Room for the title and legend
    mfrow=c(2,3)
)

matplot(t(apply(RMSEP(fit2cvA.red, estimate='adjCV', intercept=F)$val, 2, function(x) x)), main= 'A red')
matplot(t(apply(RMSEP(fit2cvB.red, estimate='adjCV', intercept=F)$val, 2, function(x) x)), main= 'B red')
matplot(t(apply(RMSEP(fit2cvA.grn, estimate='adjCV', intercept=F)$val, 2, function(x) x)))
matplot(t(apply(RMSEP(fit2cvB.grn, estimate='adjCV', intercept=F)$val, 2, function(x) x)))
matplot(t(apply(RMSEP(fit2cvA.II, estimate='adjCV', intercept=F)$val, 2, function(x) x)))
matplot(t(apply(RMSEP(fit2cvB.II, estimate='adjCV', intercept=F)$val, 2, function(x) x)))
legend(2.8,-1,c("group A", "group B"), pch = c(1,2), lty = c(1,2))

par(op)
mtext("Root mean square error of prediction ", line=2, font=2, cex=1.2)
#op <- par(usr=c(0,1,0,1), # Reset the coordinates
#          xpd=NA)         # Allow plotting outside the plot region
legend(-.1,1, # Find suitable coordinates by trial and error
       c("one", "two"), lty=1, lwd=3, col=c("navy", "orange"), box.col=NA)


pdf(file = 'graph.pdf')
par(mar=c(5.1, 4.1, 5.1, 2.1))
m <- matrix(c(1,1,1,2,3,4,5,6,7,8,8,8), nrow = 4, ncol = 3, byrow = TRUE)
#mtext("Root mean square error of prediction ", line=2, font=2, cex=1.2)
layout(mat = m, heights = c(0.2,0.35,0.35,0.2))

par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.4,"Root mean square error of prediction",cex=2,font=1.5)
#mtext("Root mean square error of prediction",cex=2,font=2)
#plot.new()
#mtext("Root mean square error of prediction ", line=2, font=2, cex=1.2)
par(mar=c(3, 3, 3, 3), mgp = c(2.0, 0.5, 0))
matplot(t(apply(RMSEP(fit2cvA.red, estimate='adjCV', intercept=F)$val, 2, function(x) x)), ylab="Error", main= 'A red', type = "l", col=1:ncmp.large, lty=1)
matplot(t(apply(RMSEP(fit2cvB.red, estimate='adjCV', intercept=F)$val, 2, function(x) x)), ylab="", main= 'B red', type = "l" , col=1:ncmp.large, lty=1)
matplot(t(apply(RMSEP(fit2cvA.grn, estimate='adjCV', intercept=F)$val, 2, function(x) x)), ylab="", main= 'A grn', type = "l" , col=1:ncmp.large, lty=1)
matplot(t(apply(RMSEP(fit2cvB.grn, estimate='adjCV', intercept=F)$val, 2, function(x) x)), xlab="Quantiles", ylab="Error", main= 'B grn', type = "l" , col=1:ncmp.large, lty=1)
matplot(t(apply(RMSEP(fit2cvA.II, estimate='adjCV', intercept=F)$val, 2, function(x) x)), xlab="Quantiles", ylab="", main= 'A II', type = "l" , col=1:ncmp.large, lty=1)
matplot(t(apply(RMSEP(fit2cvB.II, estimate='adjCV', intercept=F)$val, 2, function(x) x)), xlab="Quantiles", ylab="", main= 'B II', type = "l" , col=1:ncmp.large, lty=1)

par(mar=c(0, 0, 1, 0)) 
# c(bottom, left, top, right)
plot.new()
#plot.new()
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#plot_colors <- c("blue","black", "green", "orange", "pink")
legend(title='Number of PLS components:', x = "top",inset = 0,
       legend = 1:ncmp.large, 
       col=1:ncmp.large, lty=1, cex=1, horiz = TRUE)
#title("My 'Title' in a strange place", side = 3, line = -21, outer = TRUE)
dev.off()
