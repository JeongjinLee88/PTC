setwd("/home/leej40/Documents/PTC/Code")
## Load necessary source codes
source("TransformedOperations.R") # transformed operations
source("functions.R") # return relevant estimates

##  Case 1. A pxq matrix C from a uniform dist
set.seed(1234)
Nrow=7; Ncol=30; min=0; max=5
B <- matrix(runif(Nrow*Ncol,min = min, max = max), nrow = Nrow, ncol = Ncol)
B_norm=sqrt(apply(B^2,1,sum))
C <- B/B_norm
##  True TPDM
TPDM_X=C%*%t(C)

##  Case 2. A matrix based on an AR(1) structure
phi=0.7
Low_ar1=ar1_cov(n = 15, rho=phi)
Low_ar1[upper.tri(Low_ar1)] <-0
##  Reorder indices of the matrix
seq=1:dim(Low_ar1)[1]
ind_new=c(c(2,4),seq[-c(2,4)])
C=Low_ar1[ind_new,ind_new]
##  True TPDM
TPDM_X=C%*%t(C)
sum(diag(TPDM_X))

##  Estimation procedure
##  Set different combinations of (n,k)
n_samp <- c(1500, 2500, 4000)
frac <- c(0.1,0.05,0.02,0.01)
#frac <- c(0.1,0.05,0.02,0.01,0.007) # for AR(1)

##  Create a data frame of all 3x4 = 12 combinations
param_grid <- expand.grid(n = n_samp, frac = frac)
##  Calculate the corresponding 'k' for each 'n'
param_grid$k <- param_grid$n * param_grid$frac
##  Return all results for different (n,k)
Output <- mapply(
  estimates,
  n = param_grid$n,
  k = param_grid$k,
  MoreArgs = list(
    C = C,
    TPDM_X = TPDM_X,
    ite = 1000,
    gInv=F),SIMPLIFY = FALSE)

##  Create descriptive names for each list element
names(Output) <- paste0("n_", param_grid$n, "_frac_", param_grid$frac)
save(Output,file = "Output_unif.RData")
save(Output,file = "Output_ar1.RData")


