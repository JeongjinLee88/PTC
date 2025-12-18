estimates <- function(n, k, k_m, C, TPDM_X, ite=200, true_totalmass, gInv=FALSE){
  
  d=dim(C)[1]; Ncol=dim(C)[2]
  TPDM.h=array(NA, dim=c(d,d,ite)) # TPDM estimate
  TPDM.t_K.L=array(NA, dim=c(2,2,ite)) # conditional TPDM with true A
  TPDM.h_K.L=array(NA, dim=c(2,2,ite)) # conditional TPDM with A hat
  rho.h=rep(NA,ite) # PTC estimate
  tau_sq=rep(NA,ite) # variance estimate
  rho_var=rep(NA,ite) # PTC variance estimate
  TotalM_est=rep(NA,ite)
  
  shift <- 0.9352074  # make the mean of InvT(Z) centered.
  TPDM_KK=TPDM_X[1:2,1:2]
  TPDM_KL=TPDM_X[1:2,-(1:2)] # 2x(p-2)
  TPDM_LL=TPDM_X[-(1:2),-(1:2)] # (p-2)x(p-2)
  A=rbind(diag(2),t(-TPDM_KL%*%solve(TPDM_LL))) # true A_px2
  
  TPDM_K.L=t(A)%*%TPDM_X%*%A # true conditional TPDM
  TPDM_K.L=TPDM_KK-TPDM_KL%*%solve(TPDM_LL)%*%t(TPDM_KL) # identical
  rho=TPDM_K.L[1,2]/sqrt(TPDM_K.L[1,1]*TPDM_K.L[2,2]) # true PTC
  
  TotalM=sum(diag(TPDM_X)) # total mass
  
  #tau_sq_true=(TotalM^2)*((t(A)[1,]%*%TPDM_X%*%A[,2])^2-TPDM_K.L[1,2]^2)
  
  set.seed(1234)
  random_seed=sample(1:10000,size = ite,replace = F)
  
  for(i in 1:ite){
    
    set.seed(random_seed[i])
    
    ##  Simulate a shifted Pareto dist
    U <- runif(n*Ncol)
    Z <- matrix(1/sqrt(1-U)-shift,nrow=n,ncol=Ncol) # Necessary
    ##  Generate a vector X by matrix multiplication
    Xp <- t(Amul(C, t(Z)))
    ##  Radial/angular components
    Rad=sqrt(apply(Xp^2,1,sum))
    W=Xp/Rad
    Keep=Rad > quantile(Rad, 1-k/n)
    W_keep=W[Keep,]
    ##  Estimate a total mass
    Keep_t=Rad > sqrt(n/k_m)
    k_t=sum(Keep_t)
    TotalM_est_tmp=(k_t/k_m)
    TotalM_est[i]=TotalM_est_tmp
    #Rad_thres=quantile(Rad, Thres_t)
    #Keep_t=Rad > Rad_thres
    #k_t=sum(Keep_t)
    #TotalM_est=(k_t/n)*(Rad_thres^2)
    
    ##  Consistency as (n,k) varies
    ##  Proposition 5.1: TPDM_hat -> TPDM_X
    if(true_totalmass){
      TPDM.h_tmp=(t(W_keep)%*%W_keep / k) * TotalM  
      TPDM.h[,,i]=TPDM.h_tmp
    }else{
      TPDM.h_tmp=(t(W_keep)%*%W_keep / k) * TotalM_est_tmp  
      TPDM.h[,,i]=TPDM.h_tmp
    }
    ##  (NTD) Different Norm
    ##  (NTD) Estimation of the total mass via different norm
    
    ##  Proposition 5.3: conditional TPDM_tilde -> conditional TPDM
    TPDM.t_K.L[,,i]=t(A)%*%TPDM.h_tmp%*%A
    
    TPDM.h_KK=TPDM.h_tmp[1:2,1:2]
    TPDM.h_KL=TPDM.h_tmp[1:2,-(1:2)]
    TPDM.h_LL=TPDM.h_tmp[-(1:2),-(1:2)]
    if(gInv){
      Ahat=rbind(diag(2),t(-TPDM.h_KL%*%MASS::ginv(TPDM.h_LL)))
    }else{
      Ahat=rbind(diag(2),t(-TPDM.h_KL%*%solve(TPDM.h_LL)))
    }
    
    ##  Theorem 5.4: conditional TPDM_hat -> conditional TPDM
    if(gInv){
      ##  when ginv used
      inv_LL_part <- MASS::ginv(TPDM.h_LL)
      TPDM.h_K.L_tmp=TPDM.h_KK - TPDM.h_KL %*% inv_LL_part %*% t(TPDM.h_KL)
      TPDM.h_K.L[,,i]=TPDM.h_K.L_tmp
    }else{
      TPDM.h_K.L_tmp=t(Ahat)%*%TPDM.h_tmp%*%Ahat
      TPDM.h_K.L[,,i]=TPDM.h_K.L_tmp
    }
    
    ##  Corollary 5.5: rho.hat -> rho
    rho.h[i]=TPDM.h_K.L_tmp[1,2]/sqrt(TPDM.h_K.L_tmp[1,1]*TPDM.h_K.L_tmp[2,2])
    
    ##  Estimate variance
    if(true_totalmass){
      #Keep=Rad > quantile(Rad, 0.95)
      #W_keep=W[Keep,]
      #k_w=sum(W_keep)
      #tau_eq_tmp=(TotalM^2/k_w)*(t(Ahat)[1,]%*%t(W_keep)%*%W_keep%*%Ahat[,2])^2
      #Q1=(TotalM^2/k)*(t(Ahat)[1,]%*%t(W_keep)%*%W_keep%*%Ahat[,2])^2
      #Q2=(TotalM^2)*((t(Ahat)[1,]%*%t(W_keep)%*%W_keep%*%Ahat[,2])/k)^2
      #tau_eq_tmp=Q1-Q2
      wa1 <- W_keep%*%Ahat[,1]
      wa2 <- W_keep%*%Ahat[,2]
      term_sq_sum <- sum((wa1*wa2)^2)
      tau_eq_tmp=(TotalM^2/k)*term_sq_sum
    }else{
      #Q1=(TotalM_est_tmp^2/k)*(t(Ahat)[1,]%*%t(W_keep)%*%W_keep%*%Ahat[,2])^2
      #Q2=(TotalM_est_tmp^2)*((t(Ahat)[1,]%*%t(W_keep)%*%W_keep%*%Ahat[,2])/k)^2
      #tau_eq_tmp=Q1-Q2
      wa1 <- W_keep%*%Ahat[,1]
      wa2 <- W_keep%*%Ahat[,2]
      term_sq_sum <- sum((wa1*wa2)^2)
      tau_eq_tmp=(TotalM_est_tmp^2/k)*term_sq_sum
    }
    tau_sq[i]=tau_eq_tmp
    rho_var[i]=tau_eq_tmp/(TPDM.h_K.L_tmp[1,1]*TPDM.h_K.L_tmp[2,2])
  }
  
  return(list("TPDM_K.L"=TPDM_K.L, "rho"=rho, "TPDM.h"=TPDM.h, "TPDM.t_K.L"=TPDM.t_K.L, "TPDM.h_K.L"=TPDM.h_K.L, "rho.h"=rho.h, "tau_sq"=tau_sq, "rho_var"=rho_var, "TotalM_est"=TotalM_est))
}
