setwd("/home/leej40/Documents/PTC/Code")
source("TransformedOperations.R") # transformed operations
source("estimates.R")

##  2. Setting: a matrix from a AR(1) structure
ar1_cov <- function(n, rho){
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
phi_seq=seq(0,0.9,by=0.1)

power_coverage1=matrix(NA,length(phi_seq),3)
power_coverage2=matrix(NA,length(phi_seq),3)
for(s in 1:length(phi_seq)){
  
  phi=phi_seq[s]
  Low_ar1=ar1_cov(n = 15, rho=phi)
  Low_ar1[upper.tri(Low_ar1)] <-0
  C=Low_ar1
  TPDM_X=C%*%t(C)
  
  ##  Estimation
  ##  Different combinations of (n,k)
  n_samp <- c(1500, 2500, 4000)
  frac <- 0.1
  # 1. Create a data frame of all 3x4 = 12 combinations
  param_grid <- expand.grid(n = n_samp, frac = frac)
  # 2. Calculate the corresponding 'k' for each 'n'
  param_grid$k <- param_grid$n * param_grid$frac
  param_grid$km <- param_grid$n * 0.01
  
  Out_trueTotal <- mapply(
    estimates,
    n = param_grid$n,
    k = param_grid$k,
    k_m = param_grid$km,
    MoreArgs = list(
      C = C,
      TPDM_X = TPDM_X,
      ite = 1000,
      true_totalmass = T,
      gInv=F
    ),SIMPLIFY = FALSE)
  
  # Create descriptive names for each list element
  names(Out_trueTotal) <- paste0("n_", param_grid$n, "_frac_", param_grid$frac)
  save(Out_trueTotal,file = paste0("AR_power_",s,".RData"))
  
  ##  Asymptotic normality with (n,k)
  ##  Coverage rates under 1000 iterations
  ##  1. use the formula
  #point_true <- Out_trueTotal$n_1500_frac_0.1$TPDM_K.L[1,2]
  calculate_coverage_var <- function(n_val, k, frac) {
    
    list_name <- paste0("n_", n_val, "_frac_", frac)
    
    # Extract the data
    point_est <- Out_trueTotal[[list_name]]$TPDM.h_K.L[1,2,]
    var_est   <- Out_trueTotal[[list_name]]$tau_sq
    
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
    point_est <- Out_trueTotal[[list_name]]$TPDM.h_K.L[1,2,]
    se_est   <- sd(point_est)
    
    # Calculate confidence interval
    LB <- point_est - 1.96 * se_est
    UB <- point_est + 1.96 * se_est
    
    # Return the mean
    #mean(LB <= point_true & point_true <= UB)
    mean(LB <= 0 & 0 <= UB)
  }
  
  n=c(1500,2500,4000)
  k_values=c(150,250,400)
  frac_values=0.1
  
  coverage_out1 <- mapply(
    calculate_coverage_var,
    n_val = n,
    k = k_values,
    frac = frac_values
  )
  coverage_out2 <- mapply(
    calculate_coverage_se,
    n_val = n,
    k = k_values,
    frac = frac_values
  )
  power_coverage1[s,]=coverage_out1
  power_coverage2[s,]=coverage_out2
}
save(power_coverage1,file = "power_coverage1.RData")
save(power_coverage2,file = "power_coverage2.RData")

empiricalpower1=1-power_coverage1
empiricalpower2=1-power_coverage2

# Load necessary library
library(ggplot2)
library(reshape2) # For reshaping the matrix

phi_seq <- seq(0, 0.9, by = 0.1)
phi_seq <- seq(0.1, 0.9, by = 0.1)
n_seq   <- c(1500, 2500, 4000)
power_matrix <- empiricalpower1
power_matrix <- empiricalpower2
df <- as.data.frame(power_matrix)
colnames(df) <- paste0("n = ", n_seq)
df$phi <- phi_seq
# Reshape from wide to long format
df_long <- melt(df, id.vars = "phi", variable.name = "Sample_Size", value.name = "Power")

# 4. Generate Professional Plot
pdf("/home/leej40/Documents/PTC/Code/power1.pdf",6,6)
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
  theme_bw() + # Clean, professional theme
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









