source("https://bioconductor.org/biocLite.R")
# biocLite() 
# biocLite("topGO")
biocLite("minet") # 
biocLite("glasso")

library(topGO)
vignette("topGO")
vignette("goseq")
biocLite("GO.db")
library("GO.db")

setwd("~/Dropbox/Teaching/")
sachscow <- as.matrix(read.table())

############ 
library(glasso)
library(igraph)
library(bnlearn)
library(pcalg)
library(spacejam)

setwd("/Users/ruijuanli/Desktop/2016_summer/summer_institute/network")

sachs <- as.matrix(read.table("data.txt"))
dim(sachs)

##
## Conditional independence graphs
##

## glasso
est.1 <- glasso(s=sachs, rho=5, approx=FALSE, penalize.diagonal=FALSE)
est.1

## neighborhood selection
est.2 <- glasso(s=sachs, rho=5, approx=TRUE, penalize.diagonal=FALSE)

## plotting the graphs
A1 <- abs(est.1$wi) > 1E-8; diag(A1) <- 0
A2 <- abs(est.2$wi) > 1E-8; diag(A2) <- 0

g1 <- graph.adjacency(A1, mode="undirected")
g2 <- graph.adjacency(A2, mode="undirected")

par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(g1, layout=layout.circle(g1), main='glasso')
plot(g2, layout=layout.circle(g2), main='NS')

##
## Bayeisan networks
##
dat <- read.table('sachs.data')
p <- ncol(dat)
n <- nrow(dat)

#inf <- read.table('sachs.info')
ps <- c("praf","pmek","plcg","PIP2","PIP3","P44","pakts","PKA","PKC","P38","pjnk")
colnames(dat) <- ps

X2 <- as.data.frame(scale(dat))

## Grow-Shrink
A3 <- gs(x=X2, alpha=0.01)

## Hill climbing
A4 <- hc(X2)

compare(A3, A4)

## plot the graphs
par(mfrow=c(1,2))
plot(A3, main='Grow-Shrink')
plot(A4, main='Hill Climbing')


## pcalg
indepTest <- gaussCItest

## define sufficient statistics
suffStat <- list(C=cor(dat), n=n)

## estimate CPDAG
pc.fit <- pc(suffStat, indepTest, p, alpha=0.1, verbose=FALSE)
#plot(pc.fit, main='PC Algorithm')

A5 <- pc.fit@graph		#get a graphNEL obj
A5 <- as(A5, "matrix")	#get the adjacency matrix for graphNEL obj
g5 <- graph.adjacency(A5, mode='directed')	#create an igraph obj
V(g5)$name <- ps		#assign vertex names

##get edge modes
##NOTE: the adjmat from pcalg is in "col"(biology) format, but
##igraph reads them in "row"(cs) format, need to transpose!!
getemode4CPDAG <- function(amat, est){    
  if(est == "pcalg") amat = t(amat)
  
  emode <- 2*amat
  emode[(amat==t(amat)) & (amat!=0)] <- -2
  emode <- emode[amat!=0]
  emode[emode==-2] <- 0
  
  return(emode)
}

par(mar=c(0,0,1,0))
plot(g5, layout=layout.circle(g5), vertex.size=25, vertex.color=NA,
     main='PC Algorithm', edge.arrow.mode=getemode4CPDAG(A5,"pcalg"))

##
##drawing all estimates together
##
A3 <- amat(A3)	#GS
g3 <- graph.adjacency(A3, mode='directed')	#create an igraph obj
V(g3)$name <- ps		#assign vertex names

A4 <- amat(A4)	#HC (scaled data)
g4 <- graph.adjacency(A4, mode='directed')	#create an igraph obj
V(g4)$name <- ps		#assign vertex names

par(mfrow=c(1,3), mar=c(0,0,1,0))
plot(g5, layout=layout.circle(g5), vertex.size=25, vertex.color=NA,
     main='PC Algorithm', edge.arrow.mode=getemode4CPDAG(A5,"pcalg"))
plot(g3, layout=layout.circle(g3), vertex.size=25, vertex.color=NA,
     main='Grow-Shrink', edge.arrow.mode=getemode4CPDAG(A3,"pcalg"))
plot(g4, layout=layout.circle(g4), vertex.size=25, vertex.color=NA,
     main='Hill Climbing', edge.arrow.mode=getemode4CPDAG(A4,"pcalg"))
