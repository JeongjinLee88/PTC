##  Load packages
library(ismev)
library(evd)
library(igraph)

##  Set directory
setwd("/home/leej40/Documents/PTC/Code")
source("TransformedOperations.R") # transformed operations

##  Load Data
FileList <- list.files(path="Data")
for (i in 1:length(FileList)){
  load(paste("Data",FileList[i],sep="/"))
}
head(ComTSs)
dates <- ComTSs[,1]  ## dates
summary(dates)

##  Extract gauging stations 1-31
OriDat <- ComTSs[,-1]
d <- dim(OriDat)[2]  # d = # of variables
n <- dim(OriDat)[1]  # n = # of observations

##  Define an empirical CDF
Uhat <- apply(OriDat, 2, function(i) rank(i)/(n+1))

##  Simulate from a shifted Pareto distribution
shift=0.9352074
Xp <- apply(Uhat, 2, function(x) (1-x)^(-1/2)-shift)
frac=0.1
k=ceiling(frac*n)
##  Radial/angular
Rad=sqrt(apply(Xp^2,1,sum))
W=Xp/Rad
Keep=Rad > quantile(Rad, 1-k/n)
W_keep=W[Keep,]
##  TPDM estimate
TPDM.h=(t(W_keep)%*%W_keep / k) * d
TPDM.h_KK=TPDM.h[1:2,1:2]
TPDM.h_KL=TPDM.h[1:2,-(1:2)]
TPDM.h_LL=TPDM.h[-(1:2),-(1:2)]
Ahat=rbind(diag(2),t(-TPDM.h_KL%*%solve(TPDM.h_LL)))
##  Tail conditional covariance matrix estimate
TPDM.h_K.L=t(Ahat)%*%TPDM.h%*%Ahat
##  PTC estimate
rho.h=TPDM.h_K.L[1,2]/sqrt(TPDM.h_K.L[1,1]*TPDM.h_K.L[2,2])
##  Variance estimate
wa1 <- W_keep%*%Ahat[,1]
wa2 <- W_keep%*%Ahat[,2]
term_sq_sum <- sum((wa1*wa2)^2)
tau_sq=(d^2/k)*term_sq_sum

##  Create null matrices to store
TPDM.h_K.L_off=matrix(rep(NA,d^2),ncol=d)
Tau_sq_off=matrix(rep(NA,d^2),ncol=d)

##  Hypothesis testing for upper off-diagonals
for(s in 1:(d-1)){
  for(l in (s+1):d){
    # k=1        k=2
    # l=2,...,d/ l=3,...,d
    
    ##  Reorder variables' indices
    X_kl=cbind(Xp[,c(s,l)],Xp[,-c(s,l)])
    ##  Reorder TPDM matrices
    IndexNew=c(c(s,l),c(1:d)[-c(s,l)])
    TPDMhat=TPDM.h[IndexNew,IndexNew]
    
    TPDM.h_KK=TPDMhat[1:2,1:2]
    TPDM.h_KL=TPDMhat[1:2,-(1:2)]
    TPDM.h_LL=TPDMhat[-(1:2),-(1:2)]
    #lambda=1e-3
    #I_LL=diag(nrow(TPDM.h_LL))
    #inv_LL=solve(TPDM.h_LL+lambda*I_LL)
    #Ahat=rbind(diag(2),t(-TPDM.h_KL%*%inv_LL))
    Ahat=rbind(diag(2),t(-TPDM.h_KL%*%solve(TPDM.h_LL)))
    #Ahat=rbind(diag(2),t(-TPDM.h_KL%*%MASS::ginv(TPDM.h_LL)))
    TPDM.h_K.L=t(Ahat)%*%TPDMhat%*%Ahat
    rho.h=TPDM.h_K.L[1,2]/sqrt(TPDM.h_K.L[1,1]*TPDM.h_K.L[2,2])
    
    W_keep_kl=cbind(W_keep[,c(s,l)],W_keep[,-c(s,l)])
    TotalM=dim(TPDM.h)[1]
    ThresExc=ceiling(frac*n)
    wa1 <- W_keep_kl%*%Ahat[,1]
    wa2 <- W_keep_kl%*%Ahat[,2]
    term_sq_sum <- sum((wa1*wa2)^2)
    tau_sq=(TotalM^2/ThresExc)*term_sq_sum
    
    ##  Save a point estimate and variance estimate
    TPDM.h_K.L_off[s,l]=TPDM.h_K.L[1,2]
    Tau_sq_off[s,l]=tau_sq
    print(c(s,l))
  }
}

##  Save output
save(TPDM.h_K.L_off,file = "TPDM.h_K.L_off.Rdata")
save(Tau_sq_off,file = "Tau_sq_off.Rdata")
#load(file = "TPDM.h_K.L_off.Rdata")
#load(file = "Tau_sq_off.Rdata")

##  Matrix of test-statistics
ThresExc=ceiling(frac*n)  # number of threshold exceedances
Ngroup=d                  # number of groups or variables (31)
Ncomp=choose(Ngroup,2)    # number of pairwise comparisons (465)
zMtx=TPDM.h_K.L_off/sqrt(Tau_sq_off/ThresExc) # matrix of test statistics
z_crit=qnorm(0.05/Ncomp/2,lower.tail = F) # z-critical value

##  Create a graph induced by PTC
##  Fail to reject <=> disconnected edges
AdjMtx=zMtx
AdjMtx[abs(AdjMtx)<=z_crit]=0
AdjMtx[abs(AdjMtx)>z_crit]=1
AdjMtx[is.na(AdjMtx)]=0
AdjMtx=AdjMtx+t(AdjMtx)
rownames(AdjMtx)=colnames(AdjMtx)=as.character(1:d)
AdjMtx

##  Count connected edges
sum(AdjMtx[upper.tri(AdjMtx)])
# 25 significant edges: frac=0.1
 
##  Danube river network
load("/home/leej40/Documents/PTC/graphicalExtremes/data/danube.rda")
source("/home/leej40/Documents/PTC/Code/functions_graph.R")

# Plot the physical flow
danube_flow <- getDanubeFlowGraph()
plotDanubeIGraph(graph = danube_flow,vertexColors = "white",edgeColors = "black")

g_PTC=graph_from_adjacency_matrix(
  adjmatrix = AdjMtx,
  mode = c("undirected"),
  weighted = NULL,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
PTC_flowedges=as_edgelist(g_PTC, names = FALSE)
wMtx=zMtx
edge_weights_vec <- wMtx[PTC_flowedges]
edge_weights_vec=abs(edge_weights_vec)

# Plot only the flow graph
g_PTC=igraph::graph_from_edgelist(PTC_flowedges, FALSE)
plotDanubeIGraph(graph = g_PTC,vertexColors = "white",edgeColors = "black",edge.width = edge_weights_vec)


##  Extremal tree graph
# Fit a tree structure to the data
r_List <- list.files(path="/home/leej40/Documents/PTC/graphicalExtremes/R",pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_List,source))

Uhat <- apply(OriDat, 2, function(i) rank(i)/(n+1))
Y <- apply(Uhat, 2, function(x) (1-x)^(-1/2))
emst_fit <- emst(data = Y, method = "vario")

# Compute likelihood/ICs, and plot fitted graph, parameters
loglik_emst <- loglik_HR(
  data = Y,
  Gamma = emst_fit$Gamma,
  graph = emst_fit$graph
)
load("/home/leej40/Documents/PTC/graphicalExtremes/data/danube.rda")
source("/home/leej40/Documents/PTC/Code/functions_graph.R")
plotDanubeIGraph(graph = emst_fit$graph,vertexColors = "white",edgeColors = "black")
#plot_fitted_params(chi_hat, Gamma2chi(emst_fit$Gamma), xlab = 'Empirical')
