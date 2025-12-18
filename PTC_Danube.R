source("/home/leej40/Documents/RiverNetwork/Codes/Functions.R")
library(igraph)

StsChos <- c(1:31)
NoSt <- length(StsChos)

### Declustering and obtaining events ###
Years <- unique(as.numeric(substr(ComTSs[,1],1,4)))
U <- 0.1
Lag <- 4
Plotting <- 0
StNames <- as.character(StsChos)
YearsWithEvent <- numeric()
AllEvMat <- matrix(0,length(StNames),1)
AllRes <- list()
for (i in 1:length(Years)){
  Res <- ObtainMultVarEvents(TSs=ComTSs,U=U,Lag=Lag,Year=Years[i],StNames=StNames,Plotting=Plotting,mfrow=c(NoSt,1))  ## For Censored Likelihood
  Events <- Res$AllEvMat
  Locations <- Res$Location
  if (length(Events)>0) {
    YearsWithEvent <- c(YearsWithEvent,rep(Years[i],dim(Events)[2]))
    AllEvMat <- cbind(AllEvMat,Events)
  }
  AllRes[[i]] <- Res
}

DataEvents <- t(AllEvMat[,-1])
rownames(DataEvents) <- YearsWithEvent 


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
#OriDat <- DataEvents # independent events -> bad.
d <- dim(OriDat)[2]  # d = # of variables
n <- dim(OriDat)[1]  # n = # of observations

##  Define an empirical CDF
#Uhat <- apply(OriDat, 2, function(i) (n-rank(i)+0.5)/n)
#Uhat <- apply(OriDat, 2, function(i) (rank(i)-0.5)/n)
Uhat <- apply(OriDat, 2, function(i) rank(i)/(n+1))

##  Transform to (shifted) Pareto distribution
##  Xp: Pareto samples and Yp: Pre-image of XPa
shift=0.9352074
Xp <- apply(Uhat, 2, function(x) (1-x)^(-1/2)-shift)
#Xp <- apply(Uhat, 2, function(x) (1-x)^(-1/2))
#Yp <- InvT(x = XPa)
frac=0.01
frac=0.05  # 0.1, 0.07, 0.04, 0.01
k=ceiling(frac*n)
Rad=sqrt(apply(Xp^2,1,sum))
W=Xp/Rad
Keep=Rad > quantile(Rad, 1-k/n)
W_keep=W[Keep,]
TotalM=d
TPDM.h=(t(W_keep)%*%W_keep / k) * TotalM
TPDM.h_KK=TPDM.h[1:2,1:2]
TPDM.h_KL=TPDM.h[1:2,-(1:2)]
TPDM.h_LL=TPDM.h[-(1:2),-(1:2)]
Ahat=rbind(diag(2),t(-TPDM.h_KL%*%solve(TPDM.h_LL)))
TPDM.h_K.L=t(Ahat)%*%TPDM.h%*%Ahat
rho.h=TPDM.h_K.L[1,2]/sqrt(TPDM.h_K.L[1,1]*TPDM.h_K.L[2,2])
wa1 <- W_keep%*%Ahat[,1]
wa2 <- W_keep%*%Ahat[,2]
term_sq_sum <- sum((wa1*wa2)^2)
tau_sq=(TotalM^2/k)*term_sq_sum

point_est=TPDM.h_K.L[1,2]
var_est=tau_sq
LB <- point_est - 1.96 * sqrt(var_est) / sqrt(k)
UB <- point_est + 1.96 * sqrt(var_est) / sqrt(k)
point_true=0
LB <= point_true & point_true <= UB


##  Create null matrices to store pairwise estimates and its sample variance
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
#load(file = "GammaMtx.Rdata")
#load(file = "TauMtx.Rdata")

##  Matrix of test-statistics
#frac=0.05
ThresExc=ceiling(frac*n)  # number of threshold exceedances
Ngroup=d                  # number of groups or variables (31)
Ncomp=choose(Ngroup,2)    # number of pairwise comparisons (465)
tMtx=TPDM.h_K.L_off/sqrt(Tau_sq_off/ThresExc) # matrix of test statistics
z_crit=qnorm(0.05/Ncomp/2,lower.tail = F) # z-critical value
z_crit=Tukey_critical

##  Create extremal graphs induced by PTC
##  FTR <-> PTC is zero <-> no connections
##  Reject -> the edge is connected
AdjMtx=tMtx
AdjMtx[abs(AdjMtx)<=z_crit]=0
AdjMtx[abs(AdjMtx)>z_crit]=1
AdjMtx[is.na(AdjMtx)]=0
AdjMtx=AdjMtx+t(AdjMtx)
rownames(AdjMtx)=colnames(AdjMtx)=as.character(1:d)
AdjMtx

##  Count connected edges
sum(AdjMtx[upper.tri(AdjMtx)])
# 17 significant edges out of 465 : frac=0.05
# 25 edges: frac=0.1
# 19 edges: frac=0.07
# 8 edges: frac=0.04
# 0 edges: frac=0.01
 
##  Danube river network
load("/home/leej40/Documents/PTC/graphicalExtremes/data/danube.rda")
source("/home/leej40/Documents/PTC/Code/functions_graph.R") # transformed operations

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
wMtx=tMtx
edge_weights_vec <- wMtx[PTC_flowedges]
edge_weights_vec=abs(edge_weights_vec)

# Plot only the flow graph
g_PTC=igraph::graph_from_edgelist(PTC_flowedges, FALSE)
plotDanubeIGraph(graph = g_PTC,vertexColors = "white",edgeColors = "black",edge.width = edge_weights_vec)


##  Extremal tree graph
# Fit tree grpah to the data
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
source("/home/leej40/Documents/PTC/Code/functions_graph.R") # transformed operations
plotDanubeIGraph(graph = emst_fit$graph,vertexColors = "white",edgeColors = "black")
#plot_fitted_params(chi_hat, Gamma2chi(emst_fit$Gamma), xlab = 'Empirical')














































################################################################################
################################################################################
g=graph_from_adjacency_matrix(
  adjmatrix = AdjMtx,
  mode = c("undirected"),
  weighted = NULL,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
wMtx=tMtx
edge_indices <- as_edgelist(g, names = FALSE)
edge_weights_vec <- wMtx[edge_indices]
edge_weights_vec=abs(edge_weights_vec)
id=tkplot(g,
       canvas.width = 800,   # Width in pixels
       canvas.height = 800,  # Height in pixels
       vertex.color = NA,           # Node color
       vertex.frame.color = "black",            # Node border color
       #vertex.shape=c("circle","square"),             # One of ??none??, ??circle??, ??square??, ??csquare??, ??rectangle?? ??crectangle??, ??vrectangle??, ??pie??, ??raster??, or ??sphere??
       vertex.size=22,                          # Size of the node (default is 15)
       #vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle))
       #vertex.label=LETTERS[1:10],                    # Character vector used to label the nodes
       vertex.label.color="black",
       #vertex.label.family="Times",                   # Font family of the label (e.g.??Times??, ??Helvetica??)
       #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
       vertex.label.cex=2,                 # Font size (multiplication factor, device-dependent)
       #vertex.label.dist=0,                           # Distance between the label and the vertex
       #vertex.label.degree=0 ,                        # The position of the label in relation to the vertex (use pi)
       edge.color="black",           # Edge color
       edge.width=edge_weights_vec
       #edge.width=abs(tMtx[abs(tMtx)>z_crit]),                        # Edge width, defaults to 1
       #edge.arrow.size=1,                           # Arrow size, defaults to 1
       #edge.arrow.width=1,                          # Arrow width, defaults to 1
       #edge.lty=c("solid")                           # Line type, could be 0 or ??blank??, 1 or ??solid??, 2 or ??dashed??, 3 or ??dotted??, 4 or ??dotdash??, 5 or ??longdash??, 6 or ??twodash??
       #edge.curved=c(rep(0,5), rep(1,5))            # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
)
final_layout <- tk_coords(id)
pdf("Network_Graph.pdf", width = 6, height = 6) # Size in inches
dev.new()
par(lwd = 20)
plot(g,
     layout = final_layout,       # Use the coordinates you manually adjusted
     vertex.shape="circle",
     vertex.color = "white",
     vertex.frame.color = "red",   # Explicitly set the border color
     vertex.size = 30,
     vertex.label.color = "black",
     vertex.label.cex = 2,
     edge.color = "black",
     edge.width = edge_weights_vec,
     rescale = FALSE,              # Important: keeps your tkplot aspect ratio
     xlim = range(final_layout[,1]), 
     ylim = range(final_layout[,2]))

dev.off() # Close the file

# 1. Normalize your custom layout to the standard unit square (-1 to 1)
# This preserves the exact relative positions but fixes the rendering scale.
layout_norm <- norm_coords(final_layout, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

pdf("Network_Graph.pdf", width = 6, height = 6)
plot(g,
     layout = layout_norm,
     rescale = TRUE,                # MUST be TRUE for correct border rendering
     vertex.shape = "circle",
     vertex.color = "white",
     vertex.frame.color = "black",  # Borders will now appear crisp
     vertex.size = 25,              # Adjusted size (standard scale 0-200)
     vertex.label.color = "black",
     vertex.label.cex = 1.5,        # Adjusted for legibility
     edge.color = "black",
     edge.width = edge_weights_vec,
     xlim = c(-1, 1),               # Standard limits
     ylim = c(-1, 1))

dev.off()


dev.new()
par(mar=c(5.1,5.1,2,2))
plot(g,
     vertex.color = "grey92",           # Node color
     #vertex.frame.color = "Forestgreen",            # Node border color
     #vertex.shape=c("circle","square"),             # One of ??none??, ??circle??, ??square??, ??csquare??, ??rectangle?? ??crectangle??, ??vrectangle??, ??pie??, ??raster??, or ??sphere??
     vertex.size=22,                          # Size of the node (default is 15)
     #vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle))
     #vertex.label=LETTERS[1:10],                    # Character vector used to label the nodes
     vertex.label.color="black",
     #vertex.label.family="Times",                   # Font family of the label (e.g.??Times??, ??Helvetica??)
     #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=2,                 # Font size (multiplication factor, device-dependent)
     #vertex.label.dist=0,                           # Distance between the label and the vertex
     #vertex.label.degree=0 ,                        # The position of the label in relation to the vertex (use pi)
     edge.color="black",           # Edge color
     #edge.width=abs(AdjMtx[abs(AdjMtx)>7.1892]),                        # Edge width, defaults to 1
     #edge.arrow.size=1,                           # Arrow size, defaults to 1
     #edge.arrow.width=1,                          # Arrow width, defaults to 1
     #edge.lty=c("solid")                           # Line type, could be 0 or ??blank??, 1 or ??solid??, 2 or ??dashed??, 3 or ??dotted??, 4 or ??dotdash??, 5 or ??longdash??, 6 or ??twodash??
     #edge.curved=c(rep(0,5), rep(1,5))            # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
)



##  Return a Latex code
library(xtable)
tMtx[is.na(tMtx)]=0
xtable(x = round(tMtx+t(tMtx),2))
PvalueMtx=round(pt(q = abs(tMtx),df=ThresExc*Ngroup-Ngroup,lower.tail = FALSE),2)
save(PvalueMtx,file = "PvalueMtx.Rdata")
LB=GammaMtx-t_critical*sqrt(TauMtx)
UB=GammaMtx+t_critical*sqrt(TauMtx)
CI=round(cbind(LB,UB),3)
Count=LB<0 & UB>0

##  Tukey's procedure for pairwise comparisons
library(agricolae)
##  Find the q-value
Tukey_critical=qtukey(p = 0.95,nmeans = Ngroup,df = ThresExc*Ngroup-Ngroup)/sqrt(2)
# Tukey's critical value = 3.768605
