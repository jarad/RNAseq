
library(reshape2)
library(rjags)
library(stringr)
library(plyr)


d = read.table("../data/jd.txt", header=T)
d = d[which(apply(d,1,median)>0),]
set.seed(1)
d = d[sample(nrow(d), 200),]

d$gene = 1:nrow(d)
use = c(grep("T32", names(d)),grep("T16", names(d)),grep("gene", names(d)))

dm = melt(d[,use],id.vars=length(use)) 

dm$subject = factor(as.numeric(dm$variable))
t32 = paste("T32",c("",".1",".2",".3"),sep="")
dm$treatment = factor(ifelse(is.na(pmatch(dm$variable, t32, dup=TRUE)), "T16","T32"))



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



dat = list(y         = dm$value, 
           N         = nrow(dm), 
           G         = max(dm$gene),
           J         = max(as.numeric(dm$subject)),
           id        = as.numeric(dm$subject), 
           gene      = dm$gene, 
           treatment = as.numeric(dm$treatment))
mu = matrix(log(aggregate(y~gene+treatment,dat,mean)$y+.1), dat$G,2)

# Initial values
inits = list(alpha = rowMeans(mu), deltaNotZero = 2*(mu[,1]-mu[,2]), Z=rep(0,dat$G))

m = jags.model(textConnection(model), dat, inits, n.chains=3, n.adapt=1e5)

hyper_parms = c("alpha.theta", "alpha.sigma", "delta.pi", "delta.theta", "delta.sigma", "d", "tau0"
, "c.sigma")
parms = c(hyper_parms,"Z","alpha","delta","c","deltaNotZero")
res = coda.samples(m, parms, 1e5, thin=100)
save.image("jd-res.RData")


