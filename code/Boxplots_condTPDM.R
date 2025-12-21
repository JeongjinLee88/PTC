load(file = "Output_unif.RData")
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

##  Box-plot for conditional TPDM estimates (off-diagonal)
##  1) n=1500
# Combine your vectors into a named list
#pdf("/home/leej40/Documents/PTC/Code/condTPDM_n_1500.pdf",6,6)
par(mar=c(5.1,5.1,2,2))
data_to_plot <- list(
  "A" = Output$n_1500_frac_0.1$TPDM.h_K.L[1,2,],
  "B" = Output$n_1500_frac_0.05$TPDM.h_K.L[1,2,],
  "C" = Output$n_1500_frac_0.02$TPDM.h_K.L[1,2,],
  "D" = Output$n_1500_frac_0.01$TPDM.h_K.L[1,2,]
)
#dev.new()
boxplot(data_to_plot,
        main = " ",
        ylab = " ",
        ylim = c(-0.23,0.23),
        las = 1, # Rotates y-axis labels to be horizontal,
        cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.3,
        xaxt="n")
mtext(side = 1, # 1 = bottom axis
      at = 1,   # Position for the first boxplot
      text = expression(k == 150),
      line = 1.5, # adjust this value to move down
      cex = 1.5) 
mtext(side = 1, 
      at = 2,   # Position for the second boxplot
      text = expression(k == 75),
      line = 1.5, 
      cex = 1.5)
mtext(side = 1, 
      at = 3,   # Position for the third boxplot
      text = expression(k == 30),
      line = 1.5, 
      cex = 1.5)
mtext(side = 1, 
      at = 4,   # Position for the last boxplot
      text = expression(k == 15),
      line = 1.5, 
      cex = 1.5)
abline(h = Output$n_1500_frac_0.1$TPDM_K.L[1,2], col = "red", lty = 2, lwd = 2)
dev.off()

##  2) n=2000
#pdf("/home/leej40/Documents/PTC/Code/condTPDM_n_2500.pdf",6,6)
par(mar=c(5.1,5.1,2,2))
data_to_plot <- list(
  "A" = Output$n_2500_frac_0.1$TPDM.h_K.L[1,2,],
  "B" = Output$n_2500_frac_0.05$TPDM.h_K.L[1,2,],
  "C" = Output$n_2500_frac_0.02$TPDM.h_K.L[1,2,],
  "D" = Output$n_2500_frac_0.01$TPDM.h_K.L[1,2,]
)
#dev.new()
boxplot(data_to_plot,
        main = " ",
        ylab = " ",
        #ylim = c(-0.2,0.15),
        ylim = c(-0.23,0.23),
        las = 1, # Rotates y-axis labels to be horizontal,
        cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.3,
        xaxt="n")
mtext(side = 1, # 1 = bottom axis
      at = 1,   # Position for the first boxplot
      text = expression(k == 250),
      line = 1.5, # adjust this value to move down
      cex = 1.5) 
mtext(side = 1, 
      at = 2,   # Position for the second boxplot
      text = expression(k == 125),
      line = 1.5, 
      cex = 1.5)
mtext(side = 1, 
      at = 3,   # Position for the third boxplot
      text = expression(k == 50),
      line = 1.5, 
      cex = 1.5)
mtext(side = 1, 
      at = 4,   # Position for the last boxplot
      text = expression(k == 25),
      line = 1.5, 
      cex = 1.5)
abline(h = Output$n_2500_frac_0.1$TPDM_K.L[1,2], col = "red", lty = 2, lwd = 2)
dev.off()

##  3) n=4000
#pdf("/home/leej40/Documents/PTC/Code/condTPDM_n_4000.pdf",6,6)
par(mar=c(5.1,5.1,2,2))
data_to_plot <- list(
  "A" = Output$n_4000_frac_0.1$TPDM.h_K.L[1,2,],
  "B" = Output$n_4000_frac_0.05$TPDM.h_K.L[1,2,],
  "C" = Output$n_4000_frac_0.02$TPDM.h_K.L[1,2,],
  "D" = Output$n_4000_frac_0.01$TPDM.h_K.L[1,2,]
)
#dev.new()
boxplot(data_to_plot,
        main = " ",
        ylab = " ",
        ylim = c(-0.23,0.23),
        las = 1, # Rotates y-axis labels to be horizontal,
        cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1.3,
        xaxt="n")
mtext(side = 1, # 1 = bottom axis
      at = 1,   # Position for the first boxplot
      text = expression(k == 400),
      line = 1.5, # adjust this value to move down
      cex = 1.5) 
mtext(side = 1, 
      at = 2,   # Position for the second boxplot
      text = expression(k == 200),
      line = 1.5, 
      cex = 1.5)
mtext(side = 1, 
      at = 3,   # Position for the third boxplot
      text = expression(k == 80),
      line = 1.5, 
      cex = 1.5)
mtext(side = 1, 
      at = 4,   # Position for the last boxplot
      text = expression(k == 40),
      line = 1.5, 
      cex = 1.5)
abline(h = Output$n_1500_frac_0.1$TPDM_K.L[1,2], col = "red", lty = 2, lwd = 2)
dev.off()
