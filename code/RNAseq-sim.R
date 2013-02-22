library(plyr)

sim_gene_parameters = function(G,hyper) 
{
  # Gene specific parameters
  delta = (1-rbinom(G, 1, hyper$delta$pi))*rnorm(G, hyper$delta$theta, hyper$delta$sigma)
  alpha = (1-rbinom(G, 1, hyper$alpha$pi))*rnorm(G, hyper$alpha$theta, hyper$alpha$sigma)
  phi   =                                  rnorm(G, hyper$phi  $theta, hyper$phi  $sigma)
  sigma = hyper$d*hyper$sigma^2/rchisq(G,hyper$d)
  mu = cbind(phi-alpha,phi+delta,phi+alpha)

  return(data.frame(gene=1:G, delta=delta, alpha=alpha, phi=phi, sigma=sigma, mu=mu))
}



sim = function(gene, c) {
  df = rdply(nrow(gene), c)
  names(df)[1] = "gene"

  mu = as.matrix(gene[,grep("mu", names(gene))])
  df$mu =mu[df$gene+(df$genotype-1)*nrow(mu)]

  df$epsilon = rnorm(nrow(df), 0, gene$sigma[df$gene])
  df$lambda = exp(df$mu+df$epsilon)

  df$y = rpois(nrow(df), df$c*df$lambda)

  return(df)
}

create_c = function(n) 
{
  df_list = list()
  for (i in 1:length(n))   
  {
    df_list[[i]] = data.frame(gene=i,genotype=rep(1:3, n[i]), rep=rep(1:n[i], each=3))
  }

  df = rbind.fill(df_list)
  df$c = 1
  return(df)
}

# Simulate gene-specific parameters from proposal
hyper = list(delta=list(pi=.616,theta=-.004,sigma=sqrt(.193)),
             alpha=list(pi=.889,theta=0.011,sigma=sqrt(.058)),
             phi  =list(      theta=0,sigma=.1),
             sigma=.023, d=2.84)

G=1000
df_gene = sim_gene_parameters(G,hyper)

# Simulate data
n = rep(3,G)
c = create_c(n)
d = sim(df_gene,c)

write.csv(d, file="data.csv", row.names=F)


