
library(reshape2)
library(rjags)
library(stringr)





d = read.table("../data/jd.txt", header=T)[1:100,]

d$gene = 1:nrow(d)
use = c(grep("T32", names(d)),grep("T16", names(d)),grep("gene", names(d)))

dm = melt(d[,use],id.vars=length(use)) 

t32 = paste("T32",c("",".1",".2",".3"),sep="")
dm$treatment = factor(ifelse(is.na(pmatch(dm$variable, t32, dup=TRUE)), "T16","T32"))

dat = list(y=dm$value, N=nrow(dm), G = max(dm$gene), gene=dm$gene, treatment=as.numeric(dm$treatment))

alpha.init = log(aggregate(value~gene,dm,mean)$value+.1)
tmp = aggregate(value~gene+treatment,dm,mean)
delta.init = 

inits = list(alpha=alpha.init)
m = jags.model("../models/treatments.txt", dat, inits, n.chains=3, n.adapt=1e4)
res = coda.samples(m, c("alpha","delta","Z","deltaNotZero"), 1000)
Z = res[[1]][,!is.na(str_match(colnames(res[[1]]), "Z\\[.\\]"))]
pZ = colMeans(Z)
pZ[wm <- which.min(pZ)]

plot(res,ask=T)



res = coda.samples(m, c("d","s0","alpha.theta","alpha.sigma","delta.pi","delta.theta","delta.sigma"), 1e4, thin=10)
plot(res,ask=T)
