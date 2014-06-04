library(survival) # for survival analysis
library(glmnet) # for LASSO regularization
#library(ROCR) # for ROC analysis
#library(gpclib) # for plotting confidence intervals for x and y by calculating intersection of two polygons 
library(Hmisc) # c-index
library(plyr) # data manipulation
output.folder = "./outtemp/" # name of the folder where plots and tables will be saved (use "" for current folder, or use e.g., "./foldername/")
source("R_myfunctions.R") # functions

############
# DATA
############
# log ratios (stim to costim) of relative frequencies
#dat = read.csv("CMVdata_44subjects_log2ratio.csv", check.names = FALSE)
dat = read.csv("CMVdata_44subjects_log2ratio_UPDATED.csv", check.names = FALSE) # data updates: the 32 patients originally with CMVstutus=0 in the 44 original cohort had been updated to new censor date or death. 
rownames(dat) = dat[,1]
dat = as.matrix(dat[,-1])

# same for relative frequencies: will be used only for descriptive analysis
dat.RF = read.csv("CMVdata_44subjects_relfreq_UPDATED.csv", check.names = FALSE) # data updates: the 32 patients originally with CMVstutus=0 in the 44 original cohort had been updated to new censor date or death. 
rownames(dat.RF) = dat.RF[,1]
dat.RF = as.matrix(dat.RF[,-1])

# validation cohort (previously on prophy but got off prophy after updates)
datvalid = read.csv("CMVdata_validation18_log2ratio.csv", check.names = F, head = T)
head(datvalid[,1:5])

################## 
######## Descriptive analysis for log ratios of relative frequencies
######## survival outcome matrix Y
colnames(dat)[1:2] # [1] "offprophyCMVfreedays" "CMVstatus" 
Y = Surv(dat[,1], dat[,2]) 
## Median (off-prophylaxis) follow-up time among censored
summary(Y[Y[,2]==0,1])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 107.0   350.8   534.0   597.3   781.5  1761.0 
## Median (off-prophylaxis) follow-up time 
plot(survfit(Surv(Y[,1],1-Y[,2]) ~ 1))
survfit(Surv(Y[,1],1-Y[,2]) ~ 1) # reverse Kaplan-Meier estimate
# same as summary(survfit(Surv(Y[,1],1-Y[,2]) ~ 1))$table[5]
# records   n.max n.start  events  median 0.95LCL 0.95UCL 
# 44      44      44      32     539     502     777 
######### repeat the same thing for the validation cohort of 18 patients
## Median (off-prophylaxis) follow-up time among censored n=15
summary((datvalid$cmv_freedays - datvalid$Total_prophy_days)[datvalid$CMVstatus==0])
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 1.0    85.5   154.0   156.0   207.5   373.0 
## Median (off-prophylaxis) follow-up time 
plot(survfit(Surv(datvalid$cmv_freedays - datvalid$Total_prophy_days,1-datvalid$CMVstatus) ~ 1))
survfit(Surv(datvalid$cmv_freedays - datvalid$Total_prophy_days,1-datvalid$CMVstatus) ~ 1) # reverse Kaplan-Meier estimate
# same as summary(survfit(Surv(datvalid$cmv_freedays - datvalid$Total_prophy_days,1-datvalid$CMVstatus) ~ 1))$table[5]
# records   n.max n.start  events  median 0.95LCL 0.95UCL 
     # 18      18      18      15     154     151     274 

######## univariate analysis - concordance and score test
(conc.pv.44 <- t(sapply(data.frame(dat[, -c(1:2)], check.names = F), function(x) {s = summary(coxph(Y ~ x)); c(s$concordance, sctest = s$sctest[c("df", "pvalue")])})))[order(conc.pv.44[,4]),]
as.matrix((conc.pv2.44 <- sapply(data.frame(dat[, -c(1:2)], check.names = F), function(x) wilcox.test(x~Y[,2])$p.))[order(conc.pv2.44)])
conc.pv.44[order(conc.pv.44[,4]),]
# write.csv(conc.pv.44[order(conc.pv.44[,4]),], file = "table_univariate.csv", row.names = T)

################################## 
######## Descriptive analysis using relative frequencies (not log ratios of relative frequencies)
RF.basic = dat.RF[, 5:68] # 64 cell subsets 
RF.basic.ie1 = RF.basic[, 1:32] # 32 cell subsets
RF.basic.pp65 = RF.basic[, 33:64] # 32 cell subsets

######### bar plots for CMV-/CMV+ ratio (within 6 month) 
is.case = Y[,"status"]==1 & Y[,"time"]<180
cbind(Y)[order(Y[,"time"]),]
#: four of the 32 CMV- patients have follow-up time < 180 days (107, 124, 125, 129 days)
#: one of the 12 CMV+ patient has time-to-event = 494
which.case = which(is.case==1)
which.control = which(is.case==0)
RF.basic.colMeans = ddply(as.data.frame(RF.basic), .(is.case), colMeans) #column 1: is.case (levels 1=FALSE, 2=TRUE)
RF.basic.diff = log(RF.basic.colMeans[1, -1])-log(RF.basic.colMeans[2, -1]) 
diff.CD4.ie1 = RF.basic.diff[, grep("IE-1/CD4/", colnames(RF.basic.diff))]
diff.CD4.pp65 = RF.basic.diff[, grep("pp65/CD4/", colnames(RF.basic.diff))]
diff.CD8.ie1 = RF.basic.diff[, grep("IE-1/CD8/", colnames(RF.basic.diff))]
diff.CD8.pp65 = RF.basic.diff[, grep("pp65/CD8/", colnames(RF.basic.diff))]
colnames(diff.CD4.ie1) = gsub(".+/CD(4|8)/", "", colnames(diff.CD4.ie1))
colnames(diff.CD4.pp65) = gsub(".+/CD(4|8)/", "", colnames(diff.CD4.pp65))
colnames(diff.CD8.ie1) = gsub(".+/CD(4|8)/", "", colnames(diff.CD8.ie1))
colnames(diff.CD8.pp65) = gsub(".+/CD(4|8)/", "", colnames(diff.CD8.pp65))

combinations.ordered = c("C+I-2-T-", "C-I+2-T-", "C-I-2+T-", "C-I-2-T+", "C+I+2-T-", "C+I-2+T-", "C+I-2-T+", "C-I+2+T-", "C-I+2-T+", "C-I-2+T+", "C+I+2+T-", "C+I+2-T+", "C+I-2+T+", "C-I+2+T+", "C+I+2+T+")
combi_names_long = function(nm) { #switch from short names C+I+2+T+ to long names CD107+INFgamma+IL2+TNFalpha+
  nm = gsub("C\\+", "CD107\\+", nm)
  nm = gsub("C\\-", "CD107\\-", nm)  
  nm = gsub("I", "INF*gamma", nm)
  nm = gsub("2", "IL2", nm)
  nm = gsub("T", "TNF*alpha", nm)
  nm = gsub("$", "phantom()", nm)
  nm
}
combinations.ordered.long = combi_names_long(combinations.ordered)
##### bar plots
col1="darkgreen"; col2="blue"; col3 = "orange"; col4 = "red"; col.ordered = rep(c(col1, col2, col3, col4), c(4, 6, 4, 1))
NAMES.ARG = parse(text = combinations.ordered.long)
NAMES.ARG.space = parse(text = gsub("phantom\\(\\)", "phantom(00)", combinations.ordered.long))
lab.ratio = "Ratio of mean relative frequencies\n (CMV within 6 months to no CMV)"
barplot_ratio.3 = function(x, YAXT = T, XLAB = lab.ratio, CEX = .9) { # horizontal
  x = rev(x)
  NAMES.ARG = rev(NAMES.ARG.space)
  col.ordered = rev(col.ordered)
  bar = barplot(x, names.arg = NAMES.ARG, las = 2, xaxt = "n", yaxt = "n", col = col.ordered, space = .3, xlab = XLAB, horiz = T) 
  abline(v = 1, lty = 2)
  axis(1, tck = 0.02)
  if (YAXT) text(0, bar, labels = NAMES.ARG, cex = CEX, xpd = NA, adj = 1)
    #mtext(side = 2, line = 1, xpd = NA, at = bar)
}
pdf(file = paste(output.folder, "FigureX_ratio_of_meanRFs_pp65ANDie1_horizontal.pdf", sep = ""), width = 8, height = 8) # horizontal
par(mfrow = c(2, 2), mgp = c(1.5,0.2,0), mar = c(5,4,4,2)-c(2,4,1.5,1.5)+.1, oma = c(0,10.5,0,0))
barplot_ratio(c(t(exp(diff.CD4.pp65[, combinations.ordered]))), XLAB = "")
title(main = "CD4+ pp65 stimulation", xlab = "CMV-/CMV+ ratio")
barplot_ratio(c(t(exp(diff.CD8.pp65[, combinations.ordered]))), XLAB = "", YAXT = F)
title(main = "CD8+ pp65 stimulation", xlab = "CMV-/CMV+ ratio")
barplot_ratio(c(t(exp(diff.CD4.ie1[, combinations.ordered]))), XLAB = "")
title(main = "CD4+ IE-1 stimulation", xlab = "CMV-/CMV+ ratio")
barplot_ratio(c(t(exp(diff.CD8.ie1[, combinations.ordered]))), XLAB = "", YAXT = F)
title(main = "CD8+ IE-1 stimulation", xlab = "CMV-/CMV+ ratio")
dev.off()
#########

######################
######## log ratio variables for main analysis
######## split into groups (basic, basic.ie1, ..., maturational, ...)
# log (base 2) ratios for CD8/IFNg
LR.CD8IFNg = dat[, 3:4] # 2 cell subsets (ie1 and pp65)
# log (base 2) ratios for basic cell subsets
LR.basic = dat[, 5:68] # 64 cell subsets 
LR.basic.ie1 = LR.basic[, 1:32] # 32 cell subsets
LR.basic.pp65 = LR.basic[, 33:64] # 32 cell subsets
# log (base 2) ratios for maturational cell subsets
LR.matu = dat[, 69:388] # 320 cell subsets 
LR.matu.ie1 = LR.matu[, 1:160] # 160 cell subsets
LR.matu.pp65 = LR.matu[, 161:320] # 160 cell subsets

################
# MAIN ANALYSIS
################

######## fit the Cox model with adaptive LASSO 
# fold id for leave-one-out cross-validation for tuning regularization parameters:
foldid.loo = seq(nrow(LR.basic)) # or 1:44
# fit the model, print coefficients, and save log risk score plot (as "plot_logriskscore_....pdf"):
family = "cox"
fit.CD8IFNg = fit.finalmodel(x = LR.CD8IFNg, y = Y, family = family, plot.name = "CD8IFNg", nopenalty = TRUE, foldid.list = list(foldid.loo))
fit.CD8IFNg.ie1 = fit.finalmodel(x = LR.CD8IFNg[,"IE-1/CD8/IFNg", drop = F], y = Y, family = family, plot.name = "CD8IFNg.ie1", nopenalty = TRUE, foldid.list = list(foldid.loo))
fit.basic = fit.finalmodel(x = LR.basic, y = Y, family = family, plot.name = "basic", foldid.list = list(foldid.loo))
fit.basic.ie1 = fit.finalmodel(x = LR.basic.ie1, y = Y, family = family, plot.name = "basic_ie1", foldid.list = list(foldid.loo))
fit.basic.pp65 = fit.finalmodel(x = LR.basic.pp65, y = Y, family = family, plot.name = "basic_pp65_vline", foldid.list = list(foldid.loo), vline = -1.126087)
fit.matu = fit.finalmodel(x = LR.matu, y = Y, family = family, plot.name = "maturational", foldid.list = list(foldid.loo)) 
fit.matu.ie1 = fit.finalmodel(x = LR.matu.ie1, y = Y, family = family, plot.name = "maturational_ie1", foldid.list = list(foldid.loo)) 
fit.matu.pp65 = fit.finalmodel(x = LR.matu.pp65, y = Y, family = family, plot.name = "maturational_pp65", foldid.list = list(foldid.loo)) 
 
######## find best cutoff for the basic pp65 model
cutoff_pp65 = find.cutoff(pred = fit.basic.pp65$fitted, label = Y[, "status"], time = Y[,"time"], type.measure.cutoff = "concordance", best.only = FALSE)
# : best cutoff log risk = -1.192406 (c-index=0.8388626 with SE=0.07128725)# with unupdated data: -1.126087 
cutoff_reliableHR_pp65 = sort(fit.basic.pp65$fitted)[match(1, Y[order(fit.basic.pp65$fitted), "status"])] 
# : the smallest cutoff for which both (high- and low- risk) groups have at least one event
# (otherwise, HR estimate is unreliable)
# = -1.177128 # with unupdated data: -1.116892 
cutoff_pp65_best = max(cutoff_pp65$Best, cutoff_reliableHR_pp65)
# : # = -1.177128 # with unupdated data: -1.116892

cutoff_pp65$All[cutoff_pp65$All[,2] > cutoff_pp65$All[25,2] - cutoff_pp65$All[25,3],]
# with updated data: ll=-1.37818283, ul=-0.09193318
# with unupdated data: ll=-1.3480484, ul=-0.1073164

### vertical line(s) to be placed in plots 
#after data updates: 
vline_pp65 = mean(c(-1.19240577, -1.17712796))
vrange_pp65 = c(mean(c(-1.37818283, -1.35583526)), mean(c(-0.09193318, 0.01816602)))
# vline_pp65 = -1.192406
# vrange_pp65 = c(-1.37818283, -0.09193318)
# #before data updates: 
# vline_pp65 = mean(c(-1.126087088, -1.116891907))
# vrange_pp65 = c(mean(c(-1.348048426, -1.313782273)), mean(c(-0.107316396, -0.009552492)))
# # vline_pp65 = -1.126087088
# # vrange_pp65 = c(-1.348048426, -0.107316396)


########## coefficients and relative importance for fit.basic.pp65 (Table X)
finalcoef = fit.basic.pp65$coefficients
finalcoef.scale = scale(LR.basic.pp65[, names(finalcoef)])
finalcoef.adj = finalcoef*attr(finalcoef.scale, "scaled:scale")
finalcoef.adj.perc = abs(finalcoef.adj)/(max(abs(finalcoef.adj))) * 100
cbind(coef = as.vector(finalcoef[order(-abs(finalcoef*attr(finalcoef.scale, "scaled:scale")))]), rel.imp = finalcoef.adj.perc[order(-finalcoef.adj.perc)]) 

########### Validation data of 18 patients (medium risk, no history CMV (same characteristics as original 44 but independent of original 44)
datval = datvalid
Xb.val = as.matrix(datval[, names(finalcoef)]) %*% finalcoef
pdf(paste(output.folder, "Rplot_validation18.pdf", sep = ""))
op = par(mar = par("mar")-c(0,0,3,0))
plot.concordance(y = Surv((datval$cmv_freedays-datval$Total_prophy_days), datval$CMVstatus), fitt = Xb.val, pch = c(pch1, pch2), col = c(col1.heavy, col2.heavy), legend = c("CMV infection", "Censoring"), log.y=F)
par(op)
dev.off()
# c-index 0.88 (with original 17) 0.9230769 (with 18)
funconc(time = (datval$cmv_freedays-datval$Total_prophy_days), status = datval$CMVstatus, score = Xb.val, more = T)

######## perform 10 x stratified 5-fold cross-validation
family = "cox"
set.seed(100);cv.CD8IFNg = run.cv(x = LR.CD8IFNg, y = Y, family = family, nopenalty = TRUE, nrepeat = 10, nfolds = 5) 
set.seed(108);cv.CD8IFNg.ie1 = run.cv(x = LR.CD8IFNg[,1,drop = F], y = Y, family = family, nopenalty = TRUE, nrepeat = 10, nfolds = 5)
set.seed(101);cv.basic = run.cv(x = LR.basic, y = Y, family = family, nrepeat = 10, nfolds = 5) 
set.seed(102);cv.basic.ie1 = run.cv(x = LR.basic.ie1, y = Y, family = family, nrepeat = 10, nfolds = 5) 
set.seed(103);cv.basic.pp65 = run.cv(x = LR.basic.pp65, y = Y, family = family, nrepeat = 10, nfolds = 5) 
set.seed(104);cv.matu = run.cv(x = LR.matu, y = Y, family = family, nrepeat = 10, nfolds = 5) 
set.seed(105);cv.matu.ie1 = run.cv(x = LR.matu.ie1, y = Y, family = family, nrepeat = 10, nfolds = 5) 
set.seed(106);cv.matu.pp65 = run.cv(x = LR.matu.pp65, y = Y, family = family, nrepeat = 10, nfolds = 5) 

# resubstitution c-index etc (values saved as "table_resubstitution....txt")
my.perf(fit.object.list = list(fit.CD8IFNg.ie1, fit.CD8IFNg, fit.basic, fit.basic.ie1, fit.basic.pp65, fit.matu, fit.matu.ie1, fit.matu.pp65), 
        fit.name.list = list("CD8 IFNg IE-1", "CD8 IFNg", "basic", "basic IE-1", "basic pp65", "maturational", "maturational IE-1", "maturational pp65"), 
        prefix = "resubstitution_updated_", type.response = "time-to-event", label = Y[,"status"], timelabel = Y[, "time"], is.cv = F, plot.se = F)

# average cross-validation c-index (values saved as "table_cv_....txt", a plot saved as "plot_ROC_...pdf")
set.seed(100);my.perf(fit.object.list = list(cv.CD8IFNg, cv.basic, cv.basic.ie1, cv.basic.pp65, cv.matu, cv.matu.ie1, cv.matu.pp65), 
        fit.name.list = list("CD8 IFNg", "basic", "basic IE-1", "basic pp65", "maturational", "maturational IE-1", "maturational pp65"), 
        prefix = "cv_updated_", type.response = "time-to-event", label = Y[,"status"], timelabel = Y[, "time"], is.cv = T, plot.se = F)
set.seed(200);my.perf(fit.object.list = list(cv.CD8IFNg.ie1), 
        fit.name.list = list("CD8 IFNg IE-1"), 
        prefix = "cv_updated_", type.response = "time-to-event", label = Y[,"status"], timelabel = Y[, "time"], is.cv = T, plot.se = F)

######## Summary (c-index and plots) for paper submission (basic_pp65 for original data, mock-Quantiferon for original data, and basic_pp65 for validation data)
col.final = c("red", "blue") # c(col1.heavy, col2.heavy)
pdf(paste(output.folder, "Figure_basic_pp65_NOLINE.pdf", sep = ""), width = 5, height = 5)
op = par(mar = par("mar")-c(1,0,3,0))
plot.concordance(y = Y, fitt = fit.basic.pp65$fitted, pch = c(pch1, pch2), col = col.final, legend = c("CMV infection", "Censoring"), log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.basic.pp65$fitted, more = T) #0.8791469 (0.0870618) # before updates: 0.8746594 (0.09053199)
par(op)
dev.off()
#
pdf(paste(output.folder, "Figure_basic_pp65_CONSERVATIVE.pdf", sep = ""), width = 5, height = 5)
op = par(mar = par("mar")-c(1,0,3,0))
plot.concordance(y = Y, fitt = fit.basic.pp65$fitted, pch = c(pch1, pch2), col = col.final, legend = c("CMV infection", "Censoring"), vline = vrange_pp65[1], log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.basic.pp65$fitted, more = T) #0.8791469 (0.0870618) # before updates: 0.8746594 (0.09053199)
par(op)
dev.off()
#
pdf(paste(output.folder, "Figure_basic_pp65_BEST.pdf", sep = ""), width = 5, height = 5)
op = par(mar = par("mar")-c(1,0,3,0))
plot.concordance(y = Y, fitt = fit.basic.pp65$fitted, pch = c(pch1, pch2), col = col.final, legend = c("CMV infection", "Censoring"), vline = vline_pp65, log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.basic.pp65$fitted, more = T) #0.8791469 (0.0870618) # before updates: 0.8746594
par(op)
dev.off()
#
pdf(paste(output.folder, "Figure_basic_pp65_RANGE.pdf", sep = ""), width = 5, height = 5)
op = par(mar = par("mar")-c(1,0,3,0))
plot.concordance(y = Y, fitt = fit.basic.pp65$fitted, pch = c(pch1, pch2), col = col.final, legend = c("CMV infection", "Censoring"), vline = vline_pp65, vrange = vrange_pp65, log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.basic.pp65$fitted, more = T) #0.8791469 (0.0870618) # before updates: 0.8746594
par(op)
dev.off()
# 
pdf(paste(output.folder, "Figure_CD8_INFg.pdf", sep = ""), width = 5, height = 5)
op = par(mar = par("mar")-c(1,0,3,0))
plot.concordance(y = Y, fitt = fit.CD8IFNg$fitted, pch = c(pch1, pch2), col = col.final, legend = c("CMV infection", "Censoring"), log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.CD8IFNg$fitted, more = T) #0.5829384 (0.0870618) # before updates: 0.5940054 (0.09053199)
par(op)
dev.off()
#
pdf(paste(output.folder, "Figure_CD8_INFg_IE1.pdf", sep = ""), width = 5, height = 5)
op = par(mar = par("mar")-c(1,0,3,0))
plot.concordance(y = Y, fitt = fit.CD8IFNg.ie1$fitted, pch = c(pch1, pch2), col = col.final, legend = c("CMV infection", "Censoring"), log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.CD8IFNg.ie1$fitted, more = T) #0.5924171 (0.0870618)
par(op)
dev.off()
######
pdf(paste(output.folder, "FigureXX_Comparison.pdf", sep = ""), width = 7, height = 4)
par(mfrow = c(1, 2), mgp = c(2,0.5,0), mar = c(5,4,4,2)-c(2,4,3.5,1.5)+.1, oma = c(0,3,0,0), xpd = NA)
#par(mfrow = c(1,2), mar = c(5,4,4,2)-c(1,0,3,0)+.1)#, mar = c(5,4,4,2)-c(1,0,3,0)+.1)
plot.concordance(y = Y, fitt = fit.basic.pp65$fitted, pch = c(pch1, pch2), col = col.final, cex.axis = .8, legend = c("CMV infection", "Censoring"), log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.basic.pp65$fitted, more = T) #0.8791469 (0.0870618) # before updates: 0.8746594 (0.09053199)
plot.concordance(y = Y, fitt = fit.CD8IFNg$fitted, ylab = "", pch = c(pch1, pch2), col = col.final, yaxt = F, cex.axis = .8, legend = c("CMV infection", "Censoring"), log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.CD8IFNg$fitted, more = T) #0.5829384 (0.0870618) # before updates: 0.5940054 (0.09053199)
dev.off()
#
pdf(paste(output.folder, "FigureXXX_Cutoff.pdf", sep = ""), width = 3.8, height = 4)
par(mgp = c(2,0.5,0), mar = c(5,4,4,2)-c(2,1.1,3.5,1.5)+.1)
plot.concordance(y = Y, fitt = fit.basic.pp65$fitted, pch = c(pch1, pch2), col = col.final, cex.axis = .8, legend = c("CMV infection", "Censoring"), vline = vline_pp65, vrange = vrange_pp65, log.y=T)
funconc(time = Y[,1], status = Y[,2], score = fit.basic.pp65$fitted, more = T) #0.8791469 (0.0870618) # before updates: 0.8746594 (0.09053199)
dev.off()
#
pdf(paste(output.folder, "FigureXXXXX_Validation.pdf", sep = ""), width = 3.8, height = 4)
par(mgp = c(2,0.5,0), mar = c(5,4,4,2)-c(2,1.1,3.5,1.5)+.1)
plot.concordance(y = Surv((datval$cmv_freedays-datval$Total_prophy_days), datval$CMVstatus), fitt = Xb.val, pch = c(pch1, pch2), col = col.final, cex.axis = .8, legend = c("CMV infection", "Censoring"), vline = vline_pp65, vrange = vrange_pp65, log.y=F)
funconc(time = (datval$cmv_freedays-datval$Total_prophy_days), status = datval$CMVstatus, score = Xb.val) #0.9230769 for both before and after data updates
dev.off()
#######

######## perform bootstrap analysis for the basic.pp65 model
set.seed(1001); boot.basic.pp65 = run.boot(x = LR.basic.pp65, y = Y, B = 500, maxit = 1000000)
# This is the final figure (will be manually renamed as FIGURE XXXX)
my.dendro(boot.basic.pp65, freq.th = 50, varnames = combi_names_long(gsub("pp65/", "", colnames(LR.basic.pp65))), names.finalcoef = combi_names_long(gsub("pp65/", "", names(finalcoef))), plot.name = "basic_pp65_fullname", horiz = T, longwidth = 6, shortwidth = 4, horizmar = c(0,0,0,17)+.1, cex = 1, B = 500, grid = F, height = F, col.pos = "red", plotmath = T)

##### Power analysis and sample size determination for CTOT proposal
# estimated hazard ratio with best cutoff point
cox_pp65 = coxph(Y ~ I(fit.basic.pp65$fitted > cutoff_pp65_best))
hr_pp65 = exp(coef(cox_pp65)) # estimated hazard ratio
#after data updates: 28.36322 #before data updates: 30.30421
ci_pp65 = exp(confint(cox_pp65))
# after data updates:
   # 2.5 %   97.5 %
# 3.522205 228.4002
# before data updates:
# 2.5 %   97.5 %
# 3.849437 238.5661

cutoff_CD8IFNg = find.cutoff(pred = fit.CD8IFNg$fitted, label = Y[, "status"], time = Y[,"time"], type.measure.cutoff = "concordance", best.only = FALSE)
# best cutoff log risk = 0.3577959 # before data updates: 0.3625796 
cutoff_reliableHR_CD8IFNg = sort(fit.CD8IFNg$fitted)[match(1, Y[order(fit.CD8IFNg$fitted), "status"])] 
cutoff_CD8IFNg_best = max(cutoff_CD8IFNg$Best, cutoff_reliableHR_CD8IFNg)

cox_CD8IFNg = coxph(Y ~ I(fit.CD8IFNg$fitted > cutoff_CD8IFNg_best))
hr_CD8IFNg = exp(coef(cox_CD8IFNg)) # estimated hazard ratio
# after data updates: 4.135228 # before data updates: 2.589347
ci_CD8IFNg = exp(confint(cox_CD8IFNg))
# after data updates:
# 2.5 %   97.5 %
# 1.295074 13.20397
# before data updates:
# 2.5 %   97.5 %
#0.777368 8.624896
library(xtable)

#6-month mortality
table(upper.group = fit.basic.pp65$fitted > cutoff_pp65_best, CMV.in6mon = Y[, "status"] == 1 & Y[, "time"] < 180)
           # CMV.in6mon
# upper.group FALSE TRUE
      # FALSE    25    1
      # TRUE      8   10
#table(upper.group = fit.basic.pp65$fitted > cutoff_pp65_best, CMV.status = Y[, "status"] == 1)
cens.obs = which( Y[, "status"] == 0 & Y[, "time"] < 180 )
table(upper.group = fit.basic.pp65$fitted[-cens.obs] > cutoff_pp65_best, CMV.in6mon = Y[-cens.obs, "status"] == 1 & Y[-cens.obs, "time"] < 180)
           # CMV.in6mon
# upper.group FALSE TRUE
      # FALSE    23    1
      # TRUE      6   10
#table(upper.group = fit.basic.pp65$fitted[-cens.obs] > cutoff_pp65_best, CMV.status = Y[-cens.obs, "status"] == 1)
#after data updates:
#10/18 [1] 0.5555556
#10/16 [1] 0.625
# 18/44 [1] 0.4090909
# 16/40 [1] 0.4
#before data updates: 
#10/18 [1] 0.5555556
#10/13 [1] 0.7692308
##### upper.group: control group
##### lower.group: intervention group (reduced mortality)
# compute % reduction in mc (= tref-year mortality in ctl) given mc, and hr (of trt to ctl)
compute.r = function(hr = 1/5, mc = .6) ((1-mc)^hr - (1-mc))/mc * 100
morts = 0.6 # observed 6-month mortality in upper.group = 0.56 ~ 0.63 # before data updates: 0.56 ~ 0.77
hratios = c(5, 10, 30, 60) # hr of upper to lower = 28.36322 with a 95% CI (3.522205, 228.4002) #before data updates: 30.3 with a 95% CI (3.8, 238.6)
mylty = c(2,3,1,4); mycol = c("black", "black", "red", "black")
reds = compute.r(1/hratios) # % reduction in mc by intervention 
n = seq(20, 70, .1)
nc.props = seq(20, 50, 10)/100 # observed upper group 18/44 = 0.4090909
mycex = 1
for (nc.prop in nc.props) {
pdf(paste(output.folder, "power_highriskgroup_prop", round(nc.prop*100), "_mortality", round(morts*100),".pdf", sep = ""), width = 4, height = 4)
#op = par(mfrow = c(2, 2), oma = c(3,0,1,0), mar = c(4,3,3,1)) # number of plots per page = length of morts
op = par(mar = par("mar")-c(1,0,3,1))
for (mort in morts) {
	plot(0, 0, xlim=range(n), ylim=c(.5,1),
           xlab="Total sample size",
           ylab="Power", type="n", cex = mycex)
#    title(paste("6-month CMV in high-risk group ", round(mort*100), "%", sep = ""))
    for (i in seq_along(reds)) {
    		power = sapply(n, function(n) cpower(tref = .5, mc = mort, r = reds[i], accrual = .5, tmin = .5, nc = nc.prop*n, ni = (1-nc.prop)*n, pr = FALSE))
    		lines(n, power, lty = mylty[i], col = mycol[i])
	}
    #abline(h=c(0.8,0.9,0.95), col = "grey")
    #points(c(27.57131, 36.91024, 45.64754), c(.8,.9,.95))
	legend("bottomright", legend = paste("hazard ratio =", hratios), text.col = mycol, lty = mylty, col = mycol, cex = mycex)
}
#mtitle(paste("High-risk group sample size ", round(nc.prop*100), "%", sep = ""), ll=paste(" alpha=.05, 2-tailed"), cex.l=1, cex = 1.5)
par(op)
dev.off()
}
red = reds[3] # from hratio=30
nc.prop = .4
mort = .6
cbind(power=c(.8,.9,.95), samplesize=sapply(c(.8,.9,.95), function(power) uniroot(function(x) cpower(tref=.5, mc=mort, r=red, accrual=.5, tmin=.5, nc = nc.prop*x, ni = (1-nc.prop)*x, pr=FALSE) - power, c(1,40000))$root))
     # power samplesize
# [1,]  0.80   27.57131
# [2,]  0.90   36.91024
# [3,]  0.95   45.64754
red = reds[1] # from hratio=5
nc.prop = .4
mort = .6
cbind(power=c(.8,.9,.95), samplesize=sapply(c(.8,.9,.95), function(power) uniroot(function(x) cpower(tref=.5, mc=mort, r=red, accrual=.5, tmin=.5, nc = nc.prop*x, ni = (1-nc.prop)*x, pr=FALSE) - power, c(1,40000))$root))
     # power samplesize
# [1,]  0.80   31.36941
# [2,]  0.90   41.99483
# [3,]  0.95   51.93573





