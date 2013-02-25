

# Stan
library(rstan)
mod.stan = stan_model("model.stan")
res.stan = sampling(mod.stan, data=dat$dat); res.stan


# JAGS
library(rjags)
mod.jags = jags.model("model.jags", data=dat$dat, n.adapt=1, n.chains=4)
res.jags = coda.samples(mod.jags, c("mu","sigma"), n.iter=2000); summary(res.jags)

