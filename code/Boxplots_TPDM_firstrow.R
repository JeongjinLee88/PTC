load(file = "Output_unif.RData")
load(file = "Output_ar1.RData")
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
C=Low_ar1
##  True TPDM
TPDM_X=C%*%t(C)

##  Box-plots for TPDM estimates (the first row off-diagonals)
# Extract upper triangle data
d_set=dim(TPDM_X)[1]
TPDMest_out=Output$n_1500_frac_0.1$TPDM.h
upper_triangle_list <- list()
for (i in 1) {
  for (j in 1:d_set) {
    # Only capture elements where i < j
    if (i < j) {
      label_ij <- paste0("(", i, ",", j, ")")
      elements_ij <- TPDMest_out[i,j,]
      upper_triangle_list[[label_ij]] <- elements_ij
    }
  }
}

#pdf("/home/leej40/Documents/PTC/Code/Figures/TPDM_n_1500_frac_0.1.pdf",6,6)
#dev.new()
par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(upper_triangle_list,
        main = "",
        ylab = "",
        ylim = c(0,3),
        las = 1, # Labels can be horizontal
        cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.3
)
# Create the corresponding vector of true values
# Loop in the exact same order (i < j).
true_values_vector <- c()
for (i in 1) {
  for (j in 1:d_set) {
    if (i < j) {
      # Append the true (i,j) value from TPDM_X
      true_values_vector <- c(true_values_vector, TPDM_X[i, j])
    }
  }
}
# 'upper_triangle_list' and 'true_values_vector' are aligned
# Add the horizontal segments
num_boxes <- length(upper_triangle_list)
x_positions <- 1:num_boxes # Boxplots are at x = 1, 2, 3, ...
box_width <- 0.4           # Half-width of the segment (0.8 is full box)
segments(
  x0 = x_positions - box_width, # Vector of start x-coords
  y0 = true_values_vector,      # Vector of y-coords
  x1 = x_positions + box_width, # Vector of end x-coords
  y1 = true_values_vector,      # Vector of y-coords (same as y0)
  col = "red",
  lwd = 2
)
dev.off()