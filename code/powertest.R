setwd("/home/leej40/Documents/PTC/Code")
## Load necessary library
library(ggplot2)
library(reshape2) # For reshaping the matrix
## Load necessary source codes
source("TransformedOperations.R") # transformed operations
source("functions.R") # return relevant estimates

##  Power test for the pair (1,2) for an AR(1)
phi_seq=seq(0,0.9,by=0.1)

power_rate1=matrix(NA,length(phi_seq),3)
power_rate2=matrix(NA,length(phi_seq),3)
for(s in 1:length(phi_seq)){
  
  ##  True TPDM
  phi=phi_seq[s]
  Low_ar1=ar1_cov(n = 15, rho=phi)
  Low_ar1[upper.tri(Low_ar1)] <-0
  C=Low_ar1
  TPDM_X=C%*%t(C)
  
  ##  Estimation
  ##  Different combinations of (n,k)
  n_samp <- c(2000, 3000, 4000)
  frac <- 0.05
  k_values=n_samp*frac
  # Create a data frame of all 3x4 = 12 combinations
  param_grid <- expand.grid(n = n_samp, frac = frac)
  # Calculate the corresponding 'k' for each 'n'
  param_grid$k <- param_grid$n * param_grid$frac
  
  Output_AR1 <- mapply(
    estimates,
    n = param_grid$n,
    k = param_grid$k,
    MoreArgs = list(
      C = C,
      TPDM_X = TPDM_X,
      ite = 1000,
      gInv=F),SIMPLIFY = FALSE)
  
  # Create descriptive names for each list element
  names(Output_AR1) <- paste0("n_", param_grid$n, "_frac_", param_grid$frac)
  save(Output_AR1,file = paste0("Output_AR1_",s,".RData"))
  
  ##  Rejection rates under the alternative hypothesis under 1000 iterations
  ##  1. use the variance equation (33)
  calculate_coverage_var <- function(n_val, k, frac) {
    
    list_name <- paste0("n_", n_val, "_frac_", frac)
    
    # Extract the data
    point_est <- Output_AR1[[list_name]]$TPDM.h_K.L[1,2,]
    var_est   <- Output_AR1[[list_name]]$tau_sq
    
    # Calculate confidence interval
    LB <- point_est - 1.96 * sqrt(var_est) / sqrt(k)
    UB <- point_est + 1.96 * sqrt(var_est) / sqrt(k)
    
    # Return the mean
    #mean(LB <= point_true & point_true <= UB)
    mean(LB <= 0 & 0 <= UB)
  }
  ##  2. use empirical SE
  calculate_coverage_se <- function(n_val, k, frac) {
    
    list_name <- paste0("n_", n_val, "_frac_", frac)
    
    # Extract the data
    point_est <- Output_AR1[[list_name]]$TPDM.h_K.L[1,2,]
    se_est   <- sd(point_est)
    
    # Calculate confidence interval
    LB <- point_est - 1.96 * se_est
    UB <- point_est + 1.96 * se_est
    
    # Coverage rates
    #mean(LB <= point_true & point_true <= UB)
    mean(LB <= 0 & 0 <= UB)
  }
  coverage_out1 <- mapply(
    calculate_coverage_var,
    n_val = n_samp,
    k = k_values,
    frac = frac
  )
  coverage_out2 <- mapply(
    calculate_coverage_se,
    n_val = n_samp,
    k = k_values,
    frac = frac
  )
  power_rate1[s,]=1-coverage_out1
  power_rate2[s,]=1-coverage_out2
}
save(power_rate1,file = "power_rate1.RData")
save(power_rate2,file = "power_rate2.RData")

##  Power plots
n_seq   <- c(2000, 3000, 4000)
df <- as.data.frame(power_rate1)
colnames(df) <- paste0("n = ", n_seq)
df$phi <- phi_seq
df_long <- melt(df, id.vars = "phi", variable.name = "Sample_Size", value.name = "Power")

##  Power plots
pdf("/home/leej40/Documents/PTC/Code/power1_0.03(456)_noseq.pdf",6,6)
ggplot(df_long, aes(x = phi, y = Power, shape = Sample_Size, linetype = Sample_Size)) +
  geom_line(size = 0.8) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  scale_x_continuous(breaks = phi_seq) +
  scale_linetype_manual(values = c("dotted","longdash","solid")) +
  labs(x = "", 
       y = "",
       shape = "Sample Size", 
       linetype = "Sample Size") +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.7, 0.45), # Inset legend to save space
    legend.background = element_rect(color = "black", size = 0.1),
    legend.key.size = unit(1.8, "cm"),      # Adjusts the size of the lines/points keys
    legend.key.height = unit(0.6, "cm"),    # Specifically adjusts height of keys (reduces vertical space)
    legend.spacing.y = unit(0.2, "cm"),     # Adjusts vertical spacing between legend items
    legend.margin = margin(5, 5, 5, 5),     # Adjusts the padding inside the legend box
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 14), # Title ("Sample Size")
    legend.text = element_text(size = 14),                  # Labels ("n = 1500", etc.)
    plot.margin = margin(t = 1, r = 1, b = 0.5, l = 0.5, unit = "cm")
  )
dev.off()









