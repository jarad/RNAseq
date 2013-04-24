library(rjags)

load("jd-res.RData")

vp = function(str,max) {
  paste(str,"[",1:max,"]",sep="")
}
plot(res[,vp("c",dat$J)], ask=T)


#plot(res[,vp("alpha",dat$G)])

gr = gelman.diag(res, mult=F)
max(gr$psrf[,2])

res = window(res,start=60100)


gelman.diag(res[,hyper_parms])


alpha_summary = summary(res[,vp("alpha",dat$G)])
delta_summary = summary(res[,vp("delta",dat$G)])
deltanz_summary = summary(res[,vp("deltaNotZero",dat$G)])
z_summary = summary(res[,vp("Z",dat$G)])
c_summary     = summary(res[,vp("c",dat$J)])



summary(res[,"delta.pi"]) 

# Summary of non-differentially expressed genes
delta = deltanz_summary$quantiles[,3]
zp = z_summary$statistics[,1]
plot(delta, 1-zp, ylab="Probability of signal", xlab="Estimated non-zero difference", pch=19, cex=0.5,
     ylim=c(0,1), xlim=range(delta))
segments(deltanz_summary$quantiles[,2],
         1-zp,
         deltanz_summary$quantiles[,4],
         1-zp)





delta_summary$quantiles[,3]
ordr = order(absdelta)

par(mfrow=c(1,2))
sigs = which(absdelta>0.01)

plot(delta[ordr],1-zp[ordr], ylab="Probability of signal", xlab="Estimated difference", pch=19, cex=0.5)
text(delta[sigs], 1-zp[sigs], sigs, cex=0.5, pos=4)

plot(alpha_summary$quantiles[ordr,3], delta[ordr], pch=19, cex=1-zp[ordr], 
     xlab="Estimated signal", ylab="Estimated difference")
text(alpha_summary$quantiles[sigs,3], delta[sigs], sigs, cex=0.5, pos=4)


