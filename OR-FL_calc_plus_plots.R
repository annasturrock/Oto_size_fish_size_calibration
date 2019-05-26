library(segmented)

## 1) FL RECONSTRUCTIONS
#--------------------------------------------------------------------------------------------------------------------------
# Brokenstick FL estimation using "segmented"
#--------------------------------------------------------------------------------------------------------------------------
calib = read.csv("OR_FL_FINALforR.csv")
FL <- calib$FL; OR <- calib$OR

res.lm <- lm(FL ~ OR)
res.bs <- segmented(res.lm, seg.Z = ~ OR)

# With y intercept fixed at 30 mm so that size at first feeding (~32mm) fitted the data and reflected literature (Titus et al 2004). 
# Without this action the first segment had an negative slope (an artefact of the patchy calibration data between 25 and 35mm FL).
forced.intercept <- 30
res.lm <- lm(I(FL-forced.intercept) ~ 0 + OR)
res.bs <- segmented(res.lm, seg.Z = ~ 0 + OR)
res.lm.large <- lm(FL ~ OR, data = subset(calib, OR >= res.bs$psi[2]))

FL.bs <- predict(res.bs) + forced.intercept

#returnsOR$FL <- predict(res.bs, newdata = returnsOR) + forced.intercept

#--------------------------------------------------------------------------------------------------------------------------
# Residuals estimation
#--------------------------------------------------------------------------------------------------------------------------
resid.segm <- FL - FL.bs
resid_small <- resid.segm[OR < res.bs$psi[2]]
resid_large <- resid.segm[OR >= res.bs$psi[2]]
SD1 <- round(sd(resid.segm[OR < res.bs$psi[2]]),2)
SD2 <- round(sd(resid.segm[OR >= res.bs$psi[2]]),2)

#--------------------------------------------------------------------------------------------------------------------------
# Otoliths calibration plots
#--------------------------------------------------------------------------------------------------------------------------

tiff(filename="OR-FL calibration.tif", width = 2000, height = 1000, res=220, compression = "lzw")
par(mfrow = c(1,2), mar = c(5,5,1,1))
plot(FL ~ OR, data = calib, pch = 16, col = alpha("black", 0.7),type = "p", ylab = "Fork length (mm)", 
     xlab = "Otolith radius (µm)", las = 1, cex.lab = 1.3) #main = "Fork length prediction"
points(FL.bs[order(OR)] ~ OR[order(OR)], col = "red", pch = 16, type = "l", lwd = 2.5)
# abline(v = res.bs$psi[2], lwd = 1, col = "grey")
# abline(h = forced.intercept, lwd = 1, col = "grey")
text(198,128,"A", cex = 2)

boxplot(resid.segm[OR < res.bs$psi[2]], resid.segm[OR >= res.bs$psi[2]], names = c("Before breakpoint", "After breakpoint"), col = "light grey", 
        las = 1, ylim = c(-25, 27), xlab = "Residual error", outpch=16, outcol=alpha("black", 0.7), ylab = "Residual value", cex.lab = 1.3) #, main = "Broken stick models residuals"
text(1, -25, bquote(~ sigma[r] == .(SD1)))
text(2, -25, bquote(~ sigma[r] == .(SD2)))
text(0.58,25,"B", cex = 2)

dev.off()
