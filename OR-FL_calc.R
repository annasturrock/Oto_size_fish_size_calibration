##-------------------------------------------------------------------
## FL RECONSTRUCTIONS (go to "Selection" Project for full description/workings)
#--------------------------------------------------------------------------------------------------------------------------
# Brokenstick FL estimation using "segmented"
#--------------------------------------------------------------------------------------------------------------------------
calib = read.csv("OR_FL_FINALforR.csv")
calib = subset(calib, select=c("Sample_ID","OR","FL"))
FL <- calib$FL; OR <- calib$OR
forced.intercept <- 30
res.lm <- lm(I(FL-forced.intercept) ~ 0 + OR)
res.bs <- segmented(res.lm, seg.Z = ~ 0 + OR)
FL.bs <- predict(res.bs) + forced.intercept
plot(FL.bs,FL)

## BACKCALCULATE FL AT NATAL AND FW EXIT BASED ON OTOLITH RADIUS
FWExit_SUMMARY2 = subset(FWExit_SUMMARY, select=c("ASN","OR")) ##It kept running into problems if it had diff no. cols to the calibration
FWExit_SUMMARY2$FL <- predict(res.bs, newdata = FWExit_SUMMARY2) + forced.intercept
plot(FWExit_SUMMARY2$OR, FWExit_SUMMARY2$FL)
