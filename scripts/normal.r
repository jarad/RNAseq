library(plyr)
library(rjags)

set.seed(1)

prior = list(delta = list(pi=.8, mean=0, sd=1), alpha = list(mean=0, sd=1))

# Genes
G = 100
gene = data.frame(alpha = rnorm(G, prior$alpha$mean, prior$alpha$sd),
                  delta = (1-rbinom(G, 1, prior$delta$pi))*rnorm(G, prior$delta$mean, prior$delta$sd),
                  sigma = 1/sqrt(rgamma(G, 1, 1)))

# Observations
sim_f = function() {
  data.frame(y    = c(rnorm(G,gene$alpha-gene$delta, gene$sigma), rnorm(G,gene$alpha+gene$delta, gene$sigma)),
             trt  = rep(1:2, each=G),
             gene = 1:G)
}

d = rdply(4, sim_f)
names(d)[1] = "rep"

model = "
model {
  for (i in 1:N) {
    y[i] ~ dnorm( mu[gene[i],treatment[i]], tau[gene[i]] )
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

  d ~ dunif(0,10000)
  tau0 ~ dgamma(1,1)

  alpha.theta ~  dnorm(0, 1e-2)
  alpha.tau   <- 1/alpha.sigma^2
  alpha.sigma ~  dunif(0,100)

  delta.pi    ~  dunif(0,1)
  delta.theta ~  dnorm(0, 1)
  delta.tau   <- 1/delta.sigma^2
  delta.sigma ~  dunif(0,100)
}
"

dat = list(y = d$y, 
           treatment = d$trt,
           gene = d$gene,
           N = nrow(d),
           G = max(d$gene))  
m = jags.model(textConnection(model), dat, n.chains=3, n.adapt=1e4)
hyper_parms = c("alpha.theta", "alpha.sigma", "delta.pi", "delta.theta", "delta.sigma", "d", "tau0")
parms = c(hyper_parms,"Z","alpha","delta")
res = coda.samples(m, parms, 1e4, thin=10)


