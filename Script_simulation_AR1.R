setwd("/home/leej40/Documents/PTC/Code")
source("TransformedOperations.R") # transformed operations
##  (NTD) Different transform
source("TPDM_estimate.R")

##  1. a matrix C_pxq from a uniform dist
set.seed(1234)
Nrow=7; Ncol=30; min=0; max=5
B <- matrix(runif(Nrow*Ncol,min = min, max = max), nrow = Nrow, ncol = Ncol)
B_norm=sqrt(apply(B^2,1,sum))
C <- B/B_norm
##  True TPDM
TPDM_X=C%*%t(C)

##  2. a AR(1) correlation
ar1_cov <- function(n, rho){
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

phi=0.7
Low_ar1=ar1_cov(n = 15, rho=phi)
Low_ar1[upper.tri(Low_ar1)] <-0

##  True TPDM, inverse TPDM
#C=Low_ar1
d=dim(Low_ar1)[1]
seq=1:d
ind_new=c(c(2,4),seq[-c(2,4)])
C=Low_ar1[ind_new,ind_new]
TPDM_X=C%*%t(C)
TPDM_inv=solve(TPDM_X)
eigen(TPDM_inv)$values

sum(diag(TPDM_X))

consistency_check <- function(n, k, C, TPDM_X, ite=200,Thres_rad=0.95, Thres_total=0.99){
  
  d=dim(C)[1]; Ncol=dim(C)[2]
  TPDM.h=array(NA, dim=c(d,d,ite))
  TPDM.t_K.L=array(NA, dim=c(2,2,ite))
  TPDM.h_K.L=array(NA, dim=c(2,2,ite))
  rho.h=rep(NA,ite)
  tau_sq=rep(NA,ite)
  rho_var=rep(NA,ite)
  
  shift <- 0.9352074  # make the mean of InvT(Z) centered.
  TPDM_KK=TPDM_X[1:2,1:2]
  TPDM_KL=TPDM_X[1:2,-(1:2)] # 2x(p-2)
  TPDM_LL=TPDM_X[-(1:2),-(1:2)] # (p-2)x(p-2)
  A=rbind(diag(2),t(-TPDM_KL%*%solve(TPDM_LL))) # px2
  
  TPDM_K.L=t(A)%*%TPDM_X%*%A
  TPDM_K.L=TPDM_KK-TPDM_KL%*%solve(TPDM_LL)%*%t(TPDM_KL) # identical
  rho=TPDM_K.L[1,2]/sqrt(TPDM_K.L[1,1]*TPDM_K.L[2,2])
  
  TotalM=sum(diag(TPDM_X))
  
  set.seed(1234)
  random_seed=sample(1:10000,size = ite,replace = F)
  
  for(i in 1:ite){
    
    set.seed(random_seed[i])
    
    ##  Simulate a shifted Pareto dist
    U <- runif(n*Ncol)
    Z <- matrix(1/sqrt(1-U)-shift,nrow=n,ncol=Ncol) # Necessary
    ##  Generate a vector X by matrix multiplication
    Xp <- t(Amul(C, t(Z)))
    
    ##  Finite AR1 case
    #X_1 <- Z[,1]
    #X_2 <- vSum(v1 = Z[,2],v2 = Cmul(c = phi,x = X_1))
    #X_3 <- vSum(v1 = Z[,3],v2 = Cmul(c = phi,x = X_2))
    #X_4 <- vSum(v1 = Z[,4],v2 = Cmul(c = phi,x = X_3))
    #Xp <- cbind(X_2,X_4,X_1,X_3)
    
    ##  Consistency as (n,k) varies
    ##  Proposition 5.1: TPDM_hat -> TPDM_X
    TPDM.h_tmp=TPDM_estimate(X = Xp, Thres = 1-k/n,Thres_total = 0.99, TotalMass=TotalM)
    TPDM.h[,,i]=TPDM.h_tmp
    ##  (NTD) Different Norm
    ##  (NTD) Estimation of the total mass via different norm
    
    ##  Proposition 5.3: conditional TPDM_tilde -> conditional TPDM
    TPDM.t_K.L[,,i]=t(A)%*%TPDM.h_tmp%*%A
    
    TPDM.h_KK=TPDM.h_tmp[1:2,1:2]
    TPDM.h_KL=TPDM.h_tmp[1:2,-(1:2)]
    TPDM.h_LL=TPDM.h_tmp[-(1:2),-(1:2)]
    #Ahat=rbind(diag(2),t(-TPDM.h_KL%*%solve(TPDM.h_LL)))
    Ahat=rbind(diag(2),t(-TPDM.h_KL%*%MASS::ginv(TPDM.h_LL)))
    
    inv_LL_part <- MASS::ginv(TPDM.h_LL)
    TPDM.h_K.L_tmp=TPDM.h_KK - TPDM.h_KL %*% inv_LL_part %*% t(TPDM.h_KL)
    TPDM.h_K.L[,,i]=TPDM.h_K.L_tmp
    
    ##  Theorem 5.4: conditional TPDM_hat -> conditional TPDM
    #TPDM.h_K.L_tmp=t(Ahat)%*%TPDM.h_tmp%*%Ahat
    #TPDM.h_K.L[,,i]=TPDM.h_K.L_tmp
    
    ##  Corollary 5.5: rho.hat -> rho
    rho.h[i]=TPDM.h_K.L_tmp[1,2]/sqrt(TPDM.h_K.L_tmp[1,1]*TPDM.h_K.L_tmp[2,2])
    
    ##  TPDM estimation
    Rad=sqrt(apply(Xp^2,1,sum))
    W=Xp/Rad
    Keep=Rad > quantile(Rad, Thres_rad)
    W_keep=W[Keep,]
    
    ##  Estimate a total mass
    #Rad_thres=quantile(Rad, Thres_total)
    #Keep_t=Rad > Rad_thres
    #k=sum(Keep_t)
    #TotalMass=(k/n)*(Rad_thres^2)
    tau_eq_tmp=(TotalM^2/k)*(t(Ahat)[1,]%*%t(W_keep)%*%W_keep%*%Ahat[,2])^2
    tau_sq[i]=tau_eq_tmp
    rho_var[i]=tau_eq_tmp/(TPDM.h_K.L_tmp[1,1]*TPDM.h_K.L_tmp[2,2])
  }
  
  return(list("TPDM_K.L"=TPDM_K.L, "rho"=rho, "TPDM.h"=TPDM.h, "TPDM.t_K.L"=TPDM.t_K.L, "TPDM.h_K.L"=TPDM.h_K.L, "rho.h"=rho.h, "tau_sq"=tau_sq, "rho_var"=rho_var))
}
frac=0.005
frac=0.01
frac=0.02
frac=0.05
frac=0.1
n_samp=1000
n_samp=2000
n_samp=4000
Out1=consistency_check(n = n_samp,k = n_samp*frac,C = C,TPDM_X = TPDM_X,ite = 1000,Thres_rad = 0.95,Thres_total = 0.99)

##  1. Box-plot for conditional TPDM estimates (off-diagonal)
# Combine your vectors into a named list
data_to_plot <- list(
  "Model A" = Out1$TPDM.t_K.L[1,2,],
  "Model B" = Out1$TPDM.h_K.L[1,2,]
)
dev.new()
boxplot(data_to_plot,
        main = "Comparison of (1,2) Elements",
        ylab = "Value",
        col = c("lightblue", "lightgreen"),
        las = 1 # Rotates y-axis labels to be horizontal
)
abline(h = Out1$TPDM_K.L[1,2], col = "red", lty = 2, lwd = 2)

##  2. Box-plot for rho estimates
dev.new()
boxplot(Out1$rho.h,
        main = "rho_hat",
        ylab = "Value",
        col = c("lightblue"),
        las = 1 # Rotates y-axis labels to be horizontal
)
abline(h = Out1$rho, col = "red", lty = 2, lwd = 2)

##  3. Box-plots for TPDM estimates (off-diagonals)
# 1. Extract Upper Triangle Data
upper_triangle_list <- list()
for (i in 1:dim(C)[1]) {
  for (j in 1:dim(C)[1]) {
    # Only capture elements where i < j
    if (i < j) {
      label_ij <- paste0("(", i, ",", j, ")")
      elements_ij <- Out1$TPDM.h[i,j,]
      upper_triangle_list[[label_ij]] <- elements_ij
    }
  }
}
dev.new()
par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(upper_triangle_list,
        main = "Unique Off-Diagonal Elements (Upper Triangle)",
        ylab = "Value",
        las = 1, # Labels can be horizontal
        col = "lightgreen"
)
# Create the Corresponding Vector of True Values
# We must loop in the *exact same order* (i < j).
true_values_vector <- c()
for (i in 1:dim(C)[1]) {
  for (j in 1:dim(C)[1]) {
    if (i < j) {
      # Append the true (i,j) value from TPDM_X
      true_values_vector <- c(true_values_vector, TPDM_X[i, j])
    }
  }
}
# Now 'upper_triangle_list' and 'true_values_vector' are aligned, element-for-element.
# Second, add the horizontal segments
num_boxes <- length(upper_triangle_list)
x_positions <- 1:num_boxes # Boxplots are at x = 1, 2, 3, ...
box_width <- 0.4           # Half-width of the segment (0.8 is full box)
segments(
  x0 = x_positions - box_width, # Vector of start x-coords
  y0 = true_values_vector,      # Vector of y-coords
  x1 = x_positions + box_width, # Vector of end x-coords
  y1 = true_values_vector,      # Vector of y-coords (same as y0)
  
  # Styling
  col = "red",
  lwd = 2
)

##  Asymptotic normality with (n,k)
##  AR(1)
##  Theorem 5.7: sqrt(k)[TPDM.h_K.L]_12 ~ N(0, tau^2)
k = n_samp*frac
Standard=sqrt(k)*(Out1$TPDM.h_K.L[1,2,]-Out1$TPDM_K.L[1,2])/sqrt(Out1$tau_sq)
Standard_ext=Standard[-3<=Standard&Standard<=3]
# 1. Plot the histogram, but scaled to density (freq=FALSE)
dev.new()
hist(Standard_ext,
     xlim = c(-3, 3),
     breaks = 200,
     freq = FALSE,  # This is the key: changes y-axis to density
     main = "Histogram with Density Curve",
     xlab = "Value",
     col = "lightblue"
)
# 2. Add the density line on top
density_estimate <- density(Standard_ext, na.rm = TRUE,bw = 0.01)
lines(density_estimate,
      lwd = 2,
      col = "blue")

# Add a normal curve for comparison
curve(dnorm(x, mean=0, sd=1), col="red", lty=2, lwd=2, add=TRUE)

##  Rho
k = n_samp*frac
Standard=sqrt(k)*(Out1$rho.h-Out1$rho)/sqrt(Out1$rho_var)
dev.new()
hist(x = Standard,xlim=c(-3,3),breaks = 50)
Standard

