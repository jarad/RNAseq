library(plyr)
library(rjags)

set.seed(1)

prior = list(delta = list(pi=.8, mean=0, sd=sqrt(1.93)), alpha = list(mean=0.011, sd=sqrt(.058))

# Genes
G = 1000
gene = data.frame(alpha = rnorm(G, prior$alpha$mean, prior$alpha$sd),
                  delta = (1-rbinom(G, 1, prior$delta$pi))*rnorm(G, prior$delta$mean, prior$delta$sd),
                  sigma = 1/sqrt(rgamma(G, 2.84, 2.84*.023)))

# Observations
sim_f = function(c) {
  mn = exp(c(c[1]+gene$alpha-gene$delta, c[2]+gene$alpha+gene$delta)+rnorm(2*G, 0, gene$sigma))
  data.frame(y    = rpois(2*G,mn),
             trt  = rep(1:2, each=G),
             gene = 1:G)
}

n_reps = 4
c = matrix(rnorm(n_reps*2,0,.1),n_reps)
d = adply(c, 1, sim_f)
names(d)[1] = "rep"
d$rep = as.numeric(d$rep)
d$sampleID = d$rep+n_reps*(d$trt-1)

model = "
model {
  for (i in 1:N) {
    y[i] ~ dpois( exp(mu[gene[i],treatment[i]]+e[i]+c[id[i]]) )
    e[i] ~ dnorm( 0 , tau[gene[i]])
  }

  for (g in 1:G) {
    tau[g] ~ dgamma(d*tau0/2,d/2)

    mu[g,1] <- alpha[g]+delta[g]
    mu[g,2] <- alpha[g]-delta[g]

    alpha[g] ~ dnorm(alpha.theta, alpha.tau)
    
    delta[g] <- (1-Z[g])*deltaNotZero[g]
    Z[g] ~ dbern(delta.pi)
    deltaNotZero[g] ~ dnorm(delta.theta, delta.tau)
  }

  for (j in 1:J) {
    c[j] ~ dnorm(0,c.tau)
  }

  d ~ dunif(0,10000)
  tau0 ~ dgamma(1,1)

  alpha.theta ~  dnorm(0, 1e-2)
  alpha.tau   <- 1/alpha.sigma^2
  alpha.sigma ~  dunif(0,100)

  delta.pi    ~  dunif(0,1)
  delta.theta ~  dnorm(0, 1)
  delta.tau   <- 1/delta.sigma^2
  delta.sigma ~  dunif(0,100)

  c.tau <- 1/c.sigma^2
  c.sigma ~ dunif(0,1)
}
"

dat = list(y = d$y, 
           treatment = d$trt,
           gene = d$gene,
           id = d$sampleID,
           N = nrow(d),
           G = max(d$gene),
           J = max(d$sampleID))  
mu = matrix(log(aggregate(y~gene+treatment,dat,mean)$y), dat$G,2)

# Initial values
alpha = rowMeans(mu)
delta = 2*(mu[,1]-mu[,2])

inits = list(c = numeric(dat$J), 
             alpha = alpha,
             delta = delta)
             
# MCMC
m = jags.model(textConnection(model), dat, n.chains=3, n.adapt=1e3)
hyper_parms = c("alpha.theta", "alpha.sigma", "delta.pi", "delta.theta", "delta.sigma", "d", "tau0", "c.sigma")
parms = c(hyper_parms,"Z","alpha","delta","c")
res = coda.samples(m, parms, 1e3, thin=1)

vp = function(str,max) {
  paste(str,"[",1:max,"]",sep="")
}
plot(res[,vp("c",dat$J)], ask=T)


#plot(res[,vp("alpha",dat$G)])

gelman.diag(res[,hyper_parms])


alpha_summary = summary(res[,vp("alpha",dat$G)])
delta_summary = summary(res[,vp("delta",dat$G)])
c_summary     = summary(res[,vp("c",n_reps*2)])


c = as.numeric(c)
ordr = order(c)
plot( c[ordr],
      c_summary$quantiles[ordr,3], ylim=range(c_summary$quantiles),
      pch=19, cex=0.5)
segments(c[ordr],
         c_summary$quantiles[ordr,1],
         c[ordr],
         c_summary$quantiles[ordr,5])
abline(0,1, col="gray")


ordr = order(gene$alpha)
plot( gene$alpha[ordr],
      alpha_summary$quantiles[ordr,3], ylim=range(alpha_summary$quantiles),
      pch=19, cex=0.5)
segments(gene$alpha[ordr],
         alpha_summary$quantiles[ordr,1],
         gene$alpha[ordr],
         alpha_summary$quantiles[ordr,5])
abline(0,1, col="gray")

ordr = order(-gene$delta)
plot( -gene$delta[ordr],
      delta_summary$quantiles[ordr,3], ylim=range(delta_summary$quantiles),
      pch=19, cex=0.5)
segments(-gene$delta[ordr],
         delta_summary$quantiles[ordr,1],
         -gene$delta[ordr],
         delta_summary$quantiles[ordr,5])
abline(0,1, col="gray")


summary(res[,"delta.pi"]) 
length(which(gene$delta==0))/G


# Summary of non-differentially expressed genes
z_summary = summary(res[,vp("Z",dat$G)])
zp = z_summary$statistics[,1]
absdelta = abs(delta_summary$quantiles[,3])
ordr = order(absdelta)
plot(absdelta[ordr],1-zp[ordr], ylab="Probability of signal", xlab="Estimated size of signal")

plot(abs(gene$delta[ordr]),1-zp[ordr], pch=19, cex= exp(gene$alpha[ordr]/2),
     ylab="Probability of signal", xlab="Size of signal")



par(mfrow=c(1,2))
plot(gene$alpha[ordr], abs(gene$delta[ordr]), pch=19, cex=1-zp[ordr], 
     xlab="True signal", ylab="Absolute value of signal difference", xlim=c(-3,3), ylim=c(0,3))

plot(alpha_summary$quantiles[ordr,3], absdelta[ordr], pch=19, cex=1-zp[ordr], 
     xlab="Estimated signal", ylab="Absolute value of estimated signal difference",
     xlim=c(-3,3), ylim=c(0,3))
