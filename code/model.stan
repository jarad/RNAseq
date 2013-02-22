data {
  int<lower=0> G;
  int<lower=0> n;
  int<lower=0> y[G,3,n];
  int<lower=0> c[3,n];
}

parameters {  
  real epsilon[G,3,n];
  real phi[G];
  real alpha[G];
  real delta[G];

  real theta_alpha;
  real theta_delta;

  real<lower=0, upper=1> pa;
  real<lower=0, upper=1> pd;

  real<lower=0> sigma[G];
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_delta;
}

transformed parameters {  
  matrix[G,3] mu;
  for (g in 1:G) {
    mu[g,1] <- phi[g]-alpha[g];
    mu[g,2] <- phi[g]+delta[g];
    mu[g,3] <- phi[g]+alpha[g];
  }
}

model {
  real palpha[G];
  real pdelta[G];
  for (g in 1:G) {
    for (i in 1:3) {
      for (j in 1:n) {
        y[g,i,j] ~ poisson(c[i,j]*exp(mu[g,i]+epsilon[g,i,j]));
        epsilon[g,i,j] ~ normal(0,sigma[g]); 
      }
    }

    palpha[g] <- log(pa)   + normal_log(alpha[g], 0,           1e-6) +
                 log(1-pa) + normal_log(alpha[g], theta_alpha, sigma_alpha);
    pdelta[g] <- log(pd)   + normal_log(delta[g], 0,           1e-6) +
                 log(1-pd) + normal_log(delta[g], theta_delta, sigma_delta);

    sigma[g] ~ uniform(0,10);
  } 

  theta_alpha ~ normal(0.0,1.0);
  theta_delta ~ normal(0.0,1.0);

  sigma_alpha ~ uniform(0,10);
  sigma_delta ~ uniform(0,10);

  lp__ <- lp__ + log_sum_exp(palpha) + log_sum_exp(pdelta);
}

