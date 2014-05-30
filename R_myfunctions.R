al = 0.2
c1 = c(1,.5,0)
c2 = c(0,.5,1)
c3 = (c1+c2)/2
col1 = rgb(red = c1[1], green = c1[2], blue = c1[3], alpha = al)
col2 = rgb(red = c2[1], green = c2[2], blue = c2[3], alpha = al)
col3 = rgb(red = c3[1], green = c3[2], blue = c3[3], alpha = al)
col1.heavy = rgb(red = c1[1], green = c1[2], blue = c1[3], alpha = 1)
col2.heavy = rgb(red = c2[1], green = c2[2], blue = c2[3], alpha = 1)
pch1 = 4 #=x #15=filled square # 16=filled circle
pch2 = 3 #=plus # 2=open triangle # 17 filled triangle
# 'fitted' means fitted linear predictors or log risk scores when family == "cox" and fitted probabilities when family == "binomial"
######### fit the model, print coefficients, and save log risk score plot: 
fit.finalmodel = function(fit.object = NULL, x, y, family = c("cox", "binomial"), plot.name, standardize = FALSE, alpha.wt = 1, use.glmnet = FALSE, nopenalty = FALSE, foldid.list = NULL, nfolds = 10, nrepeat = 5, maxit = 100000, cex.lab = 1, cex.axis = 1, cex.legend = 1, legend = c("CMV infection", "free of CMV"), vline = NULL, vrange = NULL) {
	# Note that y is a Surv object (nx2 matrix)
	family = match.arg(family)
	grouped = family == "cox"
	if (family == "binomial" & is.Surv(y)) y = y[, "status"]
	if (!nopenalty) {
		if (is.null(fit.object)) {
			fit = my.adalasso.1dim(x = x, y = y, family = family, alpha.wt = alpha.wt, foldid.list = foldid.list, nfolds = nfolds, nrepeat = nrepeat, adaptive = TRUE, standardize = standardize, maxit = maxit, grouped = grouped, stratify.folds = TRUE)
		}
		print(fit$coefficients)
		fitt.lasso = fit$fitted.lasso
		fitt = fit$fitted
	} else {
		if (is.null(fit.object)) {
			if (use.glmnet) {
				fit = my.cvm.nopenalty(x = x, y = y, family = family, foldid.list = foldid.list, nrepeat = nrepeat, nfolds = nfolds, standardize = standardize, alpha.wt = 1, stratify.folds = TRUE, maxit = maxit)
			} else {
				if (family == "cox") {
					fit = coxph(y~x)
					fit$fitted = predict(fit, type="lp") # type=c("lp", "risk", "expected", "terms")
				}
				if (family == "binomial") {
					fit = glm(y~x, family = family)
					fit$fitted = predict(fit, type="response") # type = c("link", "response", "terms")					
				}				
			}
		}
		fitt = fit$fitted # fitted relative-risk for "cox" and fitted probabilities for "binomial"		
		print(fit$coef)
	}
		
	if (!nopenalty) Method = c("lasso", "adalasso") else Method = "nopenalty"

	for (method in Method[length(Method)]) {
		fitt = switch(method, lasso = fitt.lasso, adalasso = fitt, nopenalty = fitt)
		if (family == "cox") {
# 			my.col = rep(NA, nrow(y))
# 			my.col[y[,"status"]==1] = col1.heavy
# 			my.col[y[,"status"]==0] = col2.heavy
# 			my.pch = rep(NA, nrow(y))
# 			my.pch[y[,"status"]==1] = pch1
# 			my.pch[y[,"status"]==0] = pch2		
#       time.obs = log(y[,"time"]) 
# 			time.seq = log(floor(exp(seq(0, max(time.obs), length = 15))/20)*20)
# 			plot(fitt, time.obs, type = "n", ylab = "Number of off-prophylaxis CMV-free days", xlab = "log risk score", yaxt = "n")
#       axis(2, at = time.seq, labels = exp(time.seq))
#       points(fitt, time.obs, pch = my.pch, cex = pt.cex, col = my.col)
# 			#text(fitt, time.obs, label = rank(y[,"time"]), cex = text.cex)
# 			legend("topright", pch = c(pch1, pch2), pt.cex = pt.cex, legend = c("CMV infection", "Censoring"), col = c(col1.heavy, col2.heavy))
		  pdf(paste(output.folder, "plot_cox_logriskscore_", plot.name, "_", method,  ".pdf", sep = ""))
		  op = par(mar = par("mar")+c(0,1,-3,0))		  
      plot.concordance(y = y, fitt = fitt, pch = c(pch1, pch2), col = c(col1.heavy, col2.heavy), legend = legend, cex.lab = cex.lab, cex.axis = cex.axis, cex.legend = cex.legend, vline = vline, vrange = vrange)
      par(op)
			dev.off()
		}
		if (family == "binomial") {
			my.col = rep(NA, length(y))
			my.col[y==1] = col1
			my.col[y==0] = col2
			# pdf(paste(output.folder, "plot_logistic_probability_", plot.name, "_", method,  ".pdf", sep = ""), height = 4, width = 6)
			# op = par(mar = par("mar")-c(0,0,3,0))
			# plot(seq_along(fitt), fitt, type = "n", ylab = "Probability of CMV infection", ylim = c(0,1), xlab = "Subject")
			# abline(h = .5, lty = 3)
			# points(seq_along(fit), fitt, , pch = 19, cex = pt.cex, col = my.col)
			# legend("topleft", pch = 19, pt.cex = pt.cex, legend = c("CMV+", "CMV-"), col = c(col1, col2))
			# par(op)
			# dev.off()
		}
	}
	if (is.null(fit.object)) return(fit)
}

plot.concordance = function(y, fitt, pch = c(16, 17), col = c("orange","blue"), legend = c("Event", "Censoring"), log.y = T, yaxt = T, ylab = "Number of off-prophylaxis CMV-free days", time.seq.len = 15, time.seq.round = 20, pt.cex = 1.5, text.cex = .75, pos.legend = "topright", cex.lab = 1, cex.axis = 1, cex.legend = 1, vline = NULL, vrange = NULL) {
  my.col = rep(NA, nrow(y))
  my.col[y[,"status"]==1] = col[1]
  my.col[y[,"status"]==0] = col[2]
  my.pch = rep(NA, nrow(y))
  my.pch[y[,"status"]==1] = pch[1]
  my.pch[y[,"status"]==0] = pch[2] 
  if (log.y) {
#     time.obs = log(y[,"time"]) 
#     time.seq = log(floor(exp(seq(0, max(time.obs), length = time.seq.len))/time.seq.round)*time.seq.round)
#     plot(fitt, time.obs, type = "n", ylab = "Number of off-prophylaxis CMV-free days", xlab = "Log risk score", yaxt = "n", cex.lab = cex.lab, cex.axis = cex.axis)
#     axis(2, at = time.seq, labels = exp(time.seq))
#     : same as
    time.obs = y[,"time"]; hrange = c(1e-15, max(time.obs, na.rm=T)*2)
    plot(fitt, time.obs, type = "n", ylab = ylab, xlab = "Log risk score", yaxt = "n", cex.lab = cex.lab, cex.axis = cex.axis, log = "y")
  } else {
    time.obs = y[,"time"]; hrange = c(-100, max(time.obs, na.rm=T)*2)
    plot(fitt, time.obs, type = "n", ylab = ylab, xlab = "Log risk score", yaxt = "n", cex.lab = cex.lab, cex.axis = cex.axis)    
  }
  if (yaxt) axis(2, cex.axis = cex.axis)
  if (!is.null(vrange)) rect(vrange[1], hrange[1], vrange[2], hrange[2], col = grey(.95), border = "transparent")
#  if (!is.null(vrange)) rect(vrange[1], yr[1]-abs(yr[1]), vrange[2], yr[2]+abs(yr[2]), col = grey(.95), border = "transparent")
  if (!is.null(vline)) abline(v=vline, col = grey(.6))
  box()
  points(fitt, time.obs, pch = my.pch, cex = pt.cex, col = my.col)
  legend(pos.legend, pch = pch, pt.cex = pt.cex, legend = legend, cex = cex.legend, col = col, bg = "white")
}

Intersect = function(z) Reduce('intersect', z) 

Rug = function(x, col = 1, side = 1) by(data.frame(x, col), col, function(df) rug(df[, 1], col = as.character(df[1, 2]), side = side))



######### perform cross-validation
run.cv = function(x, y, family = c("cox", "binomial"), use.glmnet = FALSE, nopenalty = FALSE, scale.within = FALSE, standardize = FALSE, nrepeat = 1, nfolds = 10, foldid.outer.list = NULL, maxit = 100000) {
	# Note that y is a Surv object (nx2 matrix)
	family = match.arg(family)
	grouped = family == "cox"
	if (family == "binomial" & is.Surv(y)) {
		y = y[, "status"]
	}
	if (is.null(foldid.outer.list)) {
		foldid.outer.list = vector("list", nrepeat)
		for (rr in seq_len(nrepeat)) {
			foldid.outer.list[[rr]] = my.cv.glmnet(x = x, y = y, family = family, foldid.only = TRUE, nfolds = nfolds, stratify.folds = (nfolds < nrow(x)))
		}
	} else {
		nrepeat = length(foldid.outer.list)
		nfolds = max(foldid.outer.list[[1]])
	}
	
	fit.list = as.list(seq_len(nfolds*nrepeat))
	pred.test.list = as.list(seq(nfolds*nrepeat))
	label.test.list = as.list(seq(nfolds*nrepeat))
	pred.train.list = as.list(seq(nfolds*nrepeat))
	label.train.list = as.list(seq(nfolds*nrepeat))
	nobs.test.list = as.list(seq(nfolds*nrepeat))
  if (family == "cox") {
	  time.train.list = as.list(seq(nfolds*nrepeat))
	  time.test.list = as.list(seq(nfolds*nrepeat))  
  }
  
	for (rr in seq(nrepeat)) {
		foldid.outer = foldid.outer.list[[rr]]
		for (i in seq(nfolds)) {
			cat(paste("repeat =", rr, "of", nrepeat, "\nfold =", i, "of", nfolds, "\n"))
			omit = which(foldid.outer == i)
			x.train = x[-omit, , drop = FALSE]
			if (family == "cox") y.train = y[-omit, , drop = FALSE] else y.train = y[-omit]
			x.test = x[omit, , drop = FALSE]
			if (family =="cox") y.test = y[omit, , drop = FALSE] else y.test = y[omit]
			if (scale.within) {
				x.train = scale(x.train)
				x.test = scale(x.test, center = attr(x.train, "scaled:center"), scale = attr(x.train, "scaled:scale"))
			}	
			if (!nopenalty) {
				fit = my.adalasso.1dim(x = x.train, y = y.train, family = family, foldid.list = list(seq(nrow(x.train))), standardize = standardize, grouped = grouped, stratify.folds = TRUE, maxit = maxit) # find best lambda using (inner) LOOcv
				int = fit$intercept	
				coef = fit$coefficients
			} else {
				if (use.glmnet) {
					fit = my.cvm.nopenalty(x = x.train, y = y.train, family = family, foldid.list = list(seq(nrow(x.train))), standardize = standardize, grouped = grouped, stratify.folds = TRUE, maxit = maxit) # find best lambda using (inner) LOOcv
				} else {
					if (family == "cox") {
						fit = coxph(y.train~x.train)
						fit$fitted = predict(fit, type="lp") # type=c("lp", "risk", "expected", "terms")
					}
					if (family == "binomial") {
						fit = glm(y.train~x.train, family = family); print(fit$coef)
						fit$fitted = predict(fit, type="response") # type = c("link", "response", "terms")				
					}				
				}
				if (family == "cox") {
					int = NULL
					coef = fit$coef
				} else {
					int = fit$coef[1]
					coef = fit$coef[-1]
				}
				names(coef) = colnames(x)
			}
			fit.list[[i + nfolds*(rr-1)]] = fit							
print(coef)
      
			if (!length(coef)==0) {
				lp = as.vector(x.test[, names(coef), drop=FALSE] %*% coef + ifelse(length(int), int, 0))
				if (family == "cox") pred.test.list[[i + nfolds*(rr-1)]] = lp # log relative-risk, log risk score, or linear predictor
				if (family == "binomial") {
					exp.lp = exp(lp)
					pred.test.list[[i + nfolds*(rr-1)]] = exp.lp/(1+exp.lp) # probability of 1, same as predict(..., type = "response")
				}
			} else {
				if (family == "cox") pred.test.list[[i + nfolds*(rr-1)]] = rep(0, length(omit))
				if (family == "binomial") {
					exp.int = exp(int)
					pred.test.list[[i + nfolds*(rr-1)]] = rep(exp.int/(exp.int+1), length(omit))		
				}
			}
			if (family == "cox") {
				label.train.list[[i + nfolds*(rr-1)]] = y.train[, "status"]
				label.test.list[[i + nfolds*(rr-1)]] = y.test[, "status"]
				time.train.list[[i + nfolds*(rr-1)]] = y.train[, "time"]
				time.test.list[[i + nfolds*(rr-1)]] = y.test[, "time"]
			} else {
				label.train.list[[i + nfolds*(rr-1)]] = y.train
				label.test.list[[i + nfolds*(rr-1)]] = y.test				
			}
			nobs.test.list[[i + nfolds*(rr-1)]] = length(omit)	
		} # end of i
	} # end of rr
	#if (!nopenalty) 
		fitted = lapply(fit.list, '[[', 'fitted') 
	#	else fitted = lapply(fit.list, '[[', 'fitted')
	if (family == "cox") {
		list(pred.test.list = pred.test.list, label.test.list = label.test.list, nobs.test.list = nobs.test.list, pred.train.list = fitted, label.train.list = label.train.list, fit.list = fit.list
		, time.train.list = time.train.list, time.test.list = time.test.list
		)
	} else {
		list(pred.test.list = pred.test.list, label.test.list = label.test.list, nobs.test.list = nobs.test.list, pred.train.list = fitted, label.train.list = label.train.list, fit.list = fit.list)		
	}
}

# for final model, print resubstitution AUC
# for cross-validation results, print average AUC and plot average ROC curve (also print and plot average performance measures when threshold optimization was used)
# summary.perf = function(fit.object.list, fit.name.list, label = NULL, is.cv = TRUE, type.th = c("0.5", "optim", "none"), th.optim = TRUE) {
	# nfit = length(fit.object.list)
	# if (nfit==1) {
# #		pdf(paste(output.folder, "plot_ROC_", fit.name.list[[1]], ".pdf", sep = ""))
# #		op = par(mar = par("mar")-c(0,0,3,0))
		# my.perf(fit.object = fit.object.list[[1]], fit.name = fit.name.list[[1]], label = label, is.cv = is.cv, type.th = match.arg(type.th), th.optim = th.optim, plot.add = FALSE, plot.se = TRUE, plot.col = "black", plot.legend = TRUE)
	# } else {
		# pdf(paste(output.folder, "plot_ROC_compare_", fit.name.list[[1]], ".pdf", sep = ""))
		# op = par(mar = par("mar")-c(0,0,3,0))
		# for (i in seq(nfit)) {
			# my.perf(fit.object.list[[i]], fit.name = fit.name.list[[i]], 
			# label = label, is.cv = is.cv, type.th = match.arg(type.th), th.optim = th.optim, plot.add = ifelse(i==1, FALSE, TRUE), plot.se = FALSE, plot.col = i, plot.legend = FALSE)
		# }
		# legend("bottomright", legend = unlist(fit.name.list), lty = seq(nfit), col = seq(nfit))
	# }
# #	par(op)
# #	dev.off()
# }

# > my.perf(fit.object.list = list(cv.matu.pp65, cv.INT.basic.pp65), fit.name.list = list("maturational pp65", "basic pp65 with interactions"), prefix = "cv_raw_", label = Y[,"status"], is.cv = T, plot.se = F)
# > my.perf(fit.object.list = list(cv.matu.ie1, cv.INT.basic.pp65), fit.name.list = list("maturational ie1", "basic pp65 with interactions"), prefix = "cv_raw_", label = Y[,"status"], is.cv = T, plot.se = F)
# Error in if (result == "NaN") result = 1 : 
  # missing value where TRUE/FALSE needed
# In addition: Warning message:
# In max(ind) : no non-missing arguments to max; returning -Inf

#### AUC is a measure for assessing discrimination that ignores time-to-event. For time-to-event outcomes, use Harell's c-index. In the case of binary outcomes, Harell's c-index is the same as AUC. 
my.perf = function(fit.object.list, fit.name.list, prefix, is.cv = TRUE, type.response = NULL, label = NULL, timelabel = NULL, type.th = c("0.5", "optim", "none"), th.optim = TRUE, plot.se = TRUE) {
	if (is.null(type.response)) {
		if ("time.test.list" %in% names(fit.object.list)) {
			type.response = "time-to-event"
		} else type.response = "binary"
	} else type.response = match.arg(type.response, c("binary", "time-to-event"))
	nfit = length(fit.object.list)
	if (type.response == "time-to-event") {
		for (i in seq(nfit)) {
			fit.object = fit.object.list[[i]]
			fit.name = fit.name.list[[i]]
			if (!is.cv) {
				pred = fit.object$fitted
		        perf = c(funconc(time = timelabel, status = label, score = pred, more = T), funlogrank(time = timelabel, status = label, score = pred, more = T))
				sink(file = paste(output.folder, "table_", prefix, fit.name, ".txt", sep = ""))
				cat(
				  "c-index (with se)", perf[1:2], 
				  "\nlog rank (with p-value)", perf[3:4], 
				  "\n\nNumber of nonzero coefficients", length(fit.object$coefficients)
				)
				sink()
		    } else {
				cv = fit.object
				pred = cv$pred.test.list
				label = cv$label.test.list	
				timelabel = cv$time.test.list
				concordance = mapply(funconc, time = timelabel, status = label, score = pred)
		        logrank = mapply(funlogrank, time = timelabel, status = label, score = pred)
				cvm.conc = weighted.mean(concordance, w = unlist(cv$nobs.test.list))
				cvsd.conc = sqrt(weighted.mean((concordance-cvm.conc)^2, w = unlist(cv$nobs.test.list))/(length(concordance)-1))
				cvm.logrank = weighted.mean(logrank, w = unlist(cv$nobs.test.list))
				cvsd.logrank = sqrt(weighted.mean((logrank-cvm.logrank)^2, w = unlist(cv$nobs.test.list))/(length(logrank)-1))
						
		        nzero = sapply(cv$fit.list, function(x) length(x$coefficients))
				# cvm.nzero = mean(nzero)
				# cvsd.nzero = sd(nzero)
				cvm.nzero = weighted.mean(nzero, w = unlist(cv$nobs.test.list))
				cvsd.nzero = sqrt(weighted.mean((nzero-cvm.nzero)^2, w = unlist(cv$nobs.test.list))/(length(nzero)-1))
				
				sink(file = paste(output.folder, "table_", prefix, fit.name, ".txt", sep = ""))
				cat(
				  "c-index", c(cvm.conc, cvsd.conc),
				  "\nlog rank test", c(cvm.logrank, cvsd.logrank),  
				  "\n\nNumber of nonzero coefficients", c(cvm.nzero, cvsd.nzero)
				)
				sink()			
			} # end of else, that is, if (is.cv)
		} # end of for i		
	} # end of if (type.response == "time-to-event")
	if (type.response == "binary") {
		type.th = match.arg(type.th)
		for (i in seq(nfit)) {
			fit.object = fit.object.list[[i]]
			fit.name = fit.name.list[[i]]
			if (!is.cv) {
				pred = fit.object$fitted
			} else {
				cv = fit.object
				pred = cv$pred.test.list
				pred = lapply(pred, function(x) {if (length(unique(x))==1) x+rnorm(length(x),0,.Machine$double.eps*10) else x}) # avoid error due to ties (when there is only one unique value)
				label = cv$label.test.list	
			}
			perf.auc = performance(prediction(pred, label), 'auc')
			auc = unlist(perf.auc@y.values)
			perf.ROC = performance(prediction(pred, label), measure = 'tpr', x.measure = 'fpr') 
			if (!is.cv) {
				sink(file = paste(output.folder, "table_", prefix, fit.name, ".txt", sep = ""))
				cat("AUC", auc)
				cat("\n\nNumber of nonzero coefficients", length(fit.object$coefficients))
				
				cat("\n\n0.5")
				cat("\nsens", with.given.cutoff(obj = fit.object, label = label, type.measure = 'sens', type.measure.cutoff = "0.5", is.cv = F))
				cat("\tspec", with.given.cutoff(obj = fit.object, label = label, type.measure = 'spec', type.measure.cutoff = "0.5", is.cv = F))
				cat("\tppv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'ppv', type.measure.cutoff = "0.5", is.cv = F))
				cat("\tnpv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'npv', type.measure.cutoff = "0.5", is.cv = F))
		
				cat("\n\nMedian")
				cat("\nsens", with.given.cutoff(obj = fit.object, label = label, type.measure = 'sens', type.measure.cutoff = "median", is.cv = F))
				cat("\tspec", with.given.cutoff(obj = fit.object, label = label, type.measure = 'spec', type.measure.cutoff = "median", is.cv = F))
				cat("\tppv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'ppv', type.measure.cutoff = "median", is.cv = F))
				cat("\tnpv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'npv', type.measure.cutoff = "median", is.cv = F))
				
				cat("\n\nBalanced Accuracy")
				cat("\nsens", with.given.cutoff(obj = fit.object, label = label, type.measure = 'sens', is.cv = F))
				cat("\tspec", with.given.cutoff(obj = fit.object, label = label, type.measure = 'spec', is.cv = F))
				cat("\tppv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'ppv', is.cv = F))
				cat("\tnpv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'npv', is.cv = F))
				
				cat("\n\nF-score")
				cat("\nsens", with.given.cutoff(obj = fit.object, label = label, type.measure = 'sens', type.measure.cutoff = "fscore", is.cv = F))
				cat("\tspec", with.given.cutoff(obj = fit.object, label = label, type.measure = 'spec', type.measure.cutoff = "fscore", is.cv = F))
				cat("\tppv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'ppv', type.measure.cutoff = "fscore", is.cv = F))
				cat("\tnpv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'npv', type.measure.cutoff = "fscore", is.cv = F))
				
				cat("\n\nMCC")
				cat("\nsens", with.given.cutoff(obj = fit.object, label = label, type.measure = 'sens', type.measure.cutoff = "mcc", is.cv = F))
				cat("\tspec", with.given.cutoff(obj = fit.object, label = label, type.measure = 'spec', type.measure.cutoff = "mcc", is.cv = F))
				cat("\tppv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'ppv', type.measure.cutoff = "mcc", is.cv = F))
				cat("\tnpv", with.given.cutoff(obj = fit.object, label = label, type.measure = 'npv', type.measure.cutoff = "mcc", is.cv = F))
				sink()
			} else { # if (is.cv)
				cvm.auc = weighted.mean(auc, w = unlist(cv$nobs.test.list))
				cvsd.auc = sqrt(weighted.mean((auc-cvm.auc)^2, w = unlist(cv$nobs.test.list))/(length(auc)-1))
				nzero = sapply(cv$fit.list, function(x) length(x$coefficients))
				# cvm.nzero = mean(nzero)
				# cvsd.nzero = sd(nzero)
				cvm.nzero = weighted.mean(nzero, w = unlist(cv$nobs.test.list))
				cvsd.nzero = sqrt(weighted.mean((nzero-cvm.nzero)^2, w = unlist(cv$nobs.test.list))/(length(nzero)-1))
		
				sink(file = paste(output.folder, "table_", prefix, fit.name, ".txt", sep = ""))
				cat("AUC", c(cvm.auc, cvsd.auc))
				cat("\n\nNumber of nonzero coefficients", c(cvm.nzero, cvsd.nzero))
		
				cat("\n\n0.5")
				cat("\nsens", unlist(with.given.cutoff(obj = cv, type.measure = 'sens', type.measure.cutoff = "0.5")[c("mean", "se")]))
				cat("\tspec", unlist(with.given.cutoff(obj = cv, type.measure = 'spec', type.measure.cutoff = "0.5")[c("mean", "se")]))
				cat("\tppv", unlist(with.given.cutoff(obj = cv, type.measure = 'ppv', type.measure.cutoff = "0.5")[c("mean", "se")]))
				cat("\tnpv", unlist(with.given.cutoff(obj = cv, type.measure = 'npv', type.measure.cutoff = "0.5")[c("mean", "se")]))
				point5.tpr = unlist(with.given.cutoff(obj = cv, type.measure = 'tpr', type.measure.cutoff = "0.5")[c("mean", "se")])
				point5.fpr = unlist(with.given.cutoff(obj = cv, type.measure = 'fpr', type.measure.cutoff = "0.5")[c("mean", "se")])
		
				cat("\n\nMedian")
				cat("\nsens", unlist(with.given.cutoff(obj = cv, type.measure = 'sens', type.measure.cutoff = "median")[c("mean", "se")]))
				cat("\tspec", unlist(with.given.cutoff(obj = cv, type.measure = 'spec', type.measure.cutoff = "median")[c("mean", "se")]))
				cat("\tppv", unlist(with.given.cutoff(obj = cv, type.measure = 'ppv', type.measure.cutoff = "median")[c("mean", "se")]))
				cat("\tnpv", unlist(with.given.cutoff(obj = cv, type.measure = 'npv', type.measure.cutoff = "median")[c("mean", "se")]))
				median.tpr = unlist(with.given.cutoff(obj = cv, type.measure = 'tpr', type.measure.cutoff = "median")[c("mean", "se")])
				median.fpr = unlist(with.given.cutoff(obj = cv, type.measure = 'fpr', type.measure.cutoff = "median")[c("mean", "se")])
	
				cat("\n\nMinimum Distance")
				cat("\nsens", unlist(with.given.cutoff(obj = cv, type.measure = 'sens', type.measure.cutoff = "mindist")[c("mean", "se")]))
				cat("\tspec", unlist(with.given.cutoff(obj = cv, type.measure = 'spec', type.measure.cutoff = "mindist")[c("mean", "se")]))
				cat("\tppv", unlist(with.given.cutoff(obj = cv, type.measure = 'ppv', type.measure.cutoff = "mindist")[c("mean", "se")]))
				cat("\tnpv", unlist(with.given.cutoff(obj = cv, type.measure = 'npv', type.measure.cutoff = "mindist")[c("mean", "se")]))
				mindist.tpr = unlist(with.given.cutoff(obj = cv, type.measure = 'tpr', type.measure.cutoff = "mindist")[c("mean", "se")])
				mindist.fpr = unlist(with.given.cutoff(obj = cv, type.measure = 'fpr', type.measure.cutoff = "mindist")[c("mean", "se")])
	
				cat("\n\nHarmonic Mean")
				cat("\nsens", unlist(with.given.cutoff(obj = cv, type.measure = 'sens', type.measure.cutoff = "harmonic")[c("mean", "se")]))
				cat("\tspec", unlist(with.given.cutoff(obj = cv, type.measure = 'spec', type.measure.cutoff = "harmonic")[c("mean", "se")]))
				cat("\tppv", unlist(with.given.cutoff(obj = cv, type.measure = 'ppv', type.measure.cutoff = "harmonic")[c("mean", "se")]))
				cat("\tnpv", unlist(with.given.cutoff(obj = cv, type.measure = 'npv', type.measure.cutoff = "harmonic")[c("mean", "se")]))
				harmonic.tpr = unlist(with.given.cutoff(obj = cv, type.measure = 'tpr', type.measure.cutoff = "harmonic")[c("mean", "se")])
				harmonic.fpr = unlist(with.given.cutoff(obj = cv, type.measure = 'fpr', type.measure.cutoff = "harmonic")[c("mean", "se")])
	
				cat("\n\nAnti-harmonic Mean")
				cat("\nsens", unlist(with.given.cutoff(obj = cv, type.measure = 'sens', type.measure.cutoff = "antiharmonic")[c("mean", "se")]))
				cat("\tspec", unlist(with.given.cutoff(obj = cv, type.measure = 'spec', type.measure.cutoff = "antiharmonic")[c("mean", "se")]))
				cat("\tppv", unlist(with.given.cutoff(obj = cv, type.measure = 'ppv', type.measure.cutoff = "antiharmonic")[c("mean", "se")]))
				cat("\tnpv", unlist(with.given.cutoff(obj = cv, type.measure = 'npv', type.measure.cutoff = "antiharmonic")[c("mean", "se")]))
				antiharmonic.tpr = unlist(with.given.cutoff(obj = cv, type.measure = 'tpr', type.measure.cutoff = "antiharmonic")[c("mean", "se")])
				antiharmonic.fpr = unlist(with.given.cutoff(obj = cv, type.measure = 'fpr', type.measure.cutoff = "antiharmonic")[c("mean", "se")])
				
				cat("\n\nBalanced Accuracy")
				cat("\nsens", unlist(with.given.cutoff(obj = cv, type.measure = 'sens')[c("mean", "se")]))
				cat("\tspec", unlist(with.given.cutoff(obj = cv, type.measure = 'spec')[c("mean", "se")]))
				cat("\tppv", unlist(with.given.cutoff(obj = cv, type.measure = 'ppv')[c("mean", "se")]))
				cat("\tnpv", unlist(with.given.cutoff(obj = cv, type.measure = 'npv')[c("mean", "se")]))
				balacc.tpr = unlist(with.given.cutoff(obj = cv, type.measure = 'tpr')[c("mean", "se")])
				balacc.fpr = unlist(with.given.cutoff(obj = cv, type.measure = 'fpr')[c("mean", "se")])
				
				cat("\n\nF-score")
				cat("\nsens", unlist(with.given.cutoff(obj = cv, type.measure = 'sens', type.measure.cutoff = "fscore")[c("mean", "se")]))
				cat("\tspec", unlist(with.given.cutoff(obj = cv, type.measure = 'spec', type.measure.cutoff = "fscore")[c("mean", "se")]))
				cat("\tppv", unlist(with.given.cutoff(obj = cv, type.measure = 'ppv', type.measure.cutoff = "fscore")[c("mean", "se")]))
				cat("\tnpv", unlist(with.given.cutoff(obj = cv, type.measure = 'npv', type.measure.cutoff = "fscore")[c("mean", "se")]))
				fscore.tpr = unlist(with.given.cutoff(obj = cv, type.measure = 'tpr', type.measure.cutoff = "fscore")[c("mean", "se")])
				fscore.fpr = unlist(with.given.cutoff(obj = cv, type.measure = 'fpr', type.measure.cutoff = "fscore")[c("mean", "se")])
				
				cat("\n\nMCC")
				cat("\nsens", unlist(with.given.cutoff(obj = cv, type.measure = 'sens', type.measure.cutoff = "mcc")[c("mean", "se")]))
				cat("\tspec", unlist(with.given.cutoff(obj = cv, type.measure = 'spec', type.measure.cutoff = "mcc")[c("mean", "se")]))
				cat("\tppv", unlist(with.given.cutoff(obj = cv, type.measure = 'ppv', type.measure.cutoff = "mcc")[c("mean", "se")]))
				cat("\tnpv", unlist(with.given.cutoff(obj = cv, type.measure = 'npv', type.measure.cutoff = "mcc")[c("mean", "se")]))
				mcc.tpr = unlist(with.given.cutoff(obj = cv, type.measure = 'tpr', type.measure.cutoff = "mcc")[c("mean", "se")])
				mcc.fpr = unlist(with.given.cutoff(obj = cv, type.measure = 'fpr', type.measure.cutoff = "mcc")[c("mean", "se")])
				sink()
				TH = sort(unique(unlist(perf.ROC@alpha.values)), decreasing = TRUE)
				avex = sapply(TH, function(th) {
					extract.perf = Map(function(alpha, perf) {
					ind = max(which(alpha>=th)) 
					perf[ind]
					}, perf.ROC@alpha.values, perf.ROC@x.values)
					perf.vec = unlist(extract.perf)
					cvm.perf = weighted.mean(perf.vec, w = unlist(cv$nobs.test.list))
					cvsd.perf = sqrt(weighted.mean((perf.vec-cvm.perf)^2, w = unlist(cv$nobs.test.list))/(length(perf.vec)-1))
					c(cvm.perf, cvsd.perf)
					}	
				)
				avey = sapply(TH, function(th) {
					extract.perf = Map(function(alpha, perf) {
					ind = max(which(alpha>=th)) 
					perf[ind]
					}, perf.ROC@alpha.values, perf.ROC@y.values)
					perf.vec = unlist(extract.perf)
					cvm.perf = weighted.mean(perf.vec, w = unlist(cv$nobs.test.list))
					cvsd.perf = sqrt(weighted.mean((perf.vec-cvm.perf)^2, w = unlist(cv$nobs.test.list))/(length(perf.vec)-1))
					c(cvm.perf, cvsd.perf)
					}	
				)
				# Note: avex = fpr (= 1-spec), avey = tpr (= sens)
				# cbind(cutoff = TH, spec = 1-avex, sens = avey, dist = (0-avex)^2 + (1-avey)^2)
				if (i == 1) {
					if (nfit > 1) {
						print(nfit)
						pdf(paste(output.folder, "plot_ROC_compare_", prefix, ".pdf", sep = "")) 
					} else {
						print(nfit)
						pdf(paste(output.folder, "plot_ROC_", prefix, fit.name, ".pdf", sep = "")) 
					}
					op = par(mar = par("mar")-c(0,0,3,0))
					plot(avex[1,], avey[1,], type = 'n', xlab = perf.ROC@x.name , ylab = perf.ROC@y.name)
				}
				if (plot.se) {	
					X.Vec1 <- c(avex[1,], tail(avex[1,], 1), rev(avex[1,]), avex[1,][1])
					Y.Vec1 <- c(avey[1,]-avey[2,], tail(avey[1,]+avey[2,], 1), rev(avey[1,]+avey[2,]), (avey[1,]-avey[2,])[1])
					Y.Vec2 <- c(avey[1,], tail(avey[1,], 1), rev(avey[1,]), avey[1,][1])
					X.Vec2 <- c(avex[1,]-avex[2,], tail(avex[1,]+avex[2,], 1), rev(avex[1,]+avex[2,]), (avex[1,]-avex[2,])[1])
					p1 = as(cbind(X.Vec1, Y.Vec1), "gpc.poly")
					p2 = as(cbind(X.Vec2, Y.Vec2), "gpc.poly")
					plot(append.poly(p1, p2), poly.args = list(border = F), add = T)
					plot(p1, poly.arg = list(col = col1, border = F), add = TRUE)
					plot(p2, poly.arg = list(col = col2, border = F), add = TRUE)
				}
				lines(avex[1,], avey[1,], col = i)
				abline(a=0, b=1, lty = 3)
				if (type.th == "optim") {
					points(mindist.fpr[1], mindist.tpr[1], pch = 19, col = "red")
					arrows(mindist.fpr[1]-mindist.fpr[2], mindist.tpr[1], mindist.fpr[1]+mindist.fpr[2], mindist.tpr[1], angle = 90, lty = 1, col = "red", length = 0.05, code = 3)
					arrows(mindist.fpr[1], mindist.tpr[1]-mindist.tpr[2], mindist.fpr[1], mindist.tpr[1]+mindist.tpr[2], angle = 90, lty = 1, col = "red", length = 0.05, code = 3)
					text(mindist.fpr[1], mindist.tpr[1], labels = "Minimum distance", pos = 2, offset = 1.5, col = "red", cex = .7)
	
					points(balacc.fpr[1], balacc.tpr[1], pch = 19, col = "red")
					arrows(balacc.fpr[1]-balacc.fpr[2], balacc.tpr[1], balacc.fpr[1]+balacc.fpr[2], balacc.tpr[1], angle = 90, lty = 1, col = "red", length = 0.05, code = 3)
					arrows(balacc.fpr[1], balacc.tpr[1]-balacc.tpr[2], balacc.fpr[1], balacc.tpr[1]+balacc.tpr[2], angle = 90, lty = 1, col = "red", length = 0.05, code = 3)
					text(balacc.fpr[1], balacc.tpr[1], labels = "Balanced accuracy", pos = 2, offset = 1.5, col = "red", cex = .7)
					
					points(fscore.fpr[1], fscore.tpr[1], pch = 19, col = "red")
					arrows(fscore.fpr[1]-fscore.fpr[2], fscore.tpr[1], fscore.fpr[1]+fscore.fpr[2], fscore.tpr[1], angle = 90, lty = 1, col = "red", length = 0.05, code = 3)
					arrows(fscore.fpr[1], fscore.tpr[1]-fscore.tpr[2], fscore.fpr[1], fscore.tpr[1]+fscore.tpr[2], angle = 90, lty = 1, col = "red", length = 0.05, code = 3)
					text(fscore.fpr[1], fscore.tpr[1], labels = "F-score", pos = 2, offset = 1.5, col = "red", cex = .7)
					
					points(mcc.fpr[1], mcc.tpr[1], pch = 19, col = "red")
					arrows(mcc.fpr[1]-mcc.fpr[2], mcc.tpr[1], mcc.fpr[1]+mcc.fpr[2], mcc.tpr[1], angle = 90, lty = 1, col = "red", length = 0.05, code = 3)
					arrows(mcc.fpr[1], mcc.tpr[1]-mcc.tpr[2], mcc.fpr[1], mcc.tpr[1]+mcc.tpr[2], angle = 90, lty = 1, col = "red", length = 0.05, code = 3)
					text(mcc.fpr[1], mcc.tpr[1], labels = "MCC", pos = 2, offset = 1.5, col = "red", cex = .7)
					if (nfit == 1) {
						legend("bottomright",legend = c(paste("(Average ", tolower(perf.ROC@x.name), ", Average ", tolower(perf.ROC@y.name),")", sep = ""), paste("Average ", c(perf.ROC@x.name, perf.ROC@y.name)," +/- 1SE", sep = ""), "With threshold optimization"), lty = c(1, 1, 1, NA), lwd = c(par("lwd"), 10,10, NA), col = c(par("col"),col2, col1, "red"), pch = c(NA, NA, NA, 19))
					} 	
				} else if (type.th == "0.5") {
					points(point5.fpr[1], point5.tpr[1], pch = 19, col = ifelse(nfit==1, "red", i))
					arrows(point5.fpr[1]-point5.fpr[2], point5.tpr[1], point5.fpr[1]+point5.fpr[2], point5.tpr[1], angle = 90, lty = 1, col = ifelse(nfit==1, "red", i), length = 0.05, code = 3)
					arrows(point5.fpr[1], point5.tpr[1]-point5.tpr[2], point5.fpr[1], point5.tpr[1]+point5.tpr[2], angle = 90, lty = 1, col = ifelse(nfit==1, "red", i), length = 0.05, code = 3)
					if (nfit == 1) {
						legend("bottomright",legend = c(paste("(Average ", tolower(perf.ROC@x.name), ", Average ", tolower(perf.ROC@y.name),")", sep = ""), paste("Average ", c(perf.ROC@x.name, perf.ROC@y.name)," +/- 1SE", sep = ""), "With cutoff = 0.5"), lty = c(1, 1, 1, NA), lwd = c(par("lwd"), 10,10, NA), col = c(par("col"),col2, col1, "red"), pch = c(NA, NA, NA, 19))
					}	
				} else { 
					if (nfit == 1) {
						legend("bottomright",legend = c(paste("(Average ", tolower(perf.ROC@x.name), ", Average ", tolower(perf.ROC@y.name),")", sep = ""), paste("Average ", c(perf.ROC@x.name, perf.ROC@y.name)," +/- 1SE", sep = "")), lty = c(1, 1, 1), lwd = c(par("lwd"), 10,10), col = c(par("col"),col2, col1))	
					}
				}
			}	
		} # end of for i
		if (is.cv) {
			if (nfit > 1) legend("bottomright", legend = unlist(fit.name.list), lty = seq(nfit), col = seq(nfit))
			par(op)
			dev.off()
		}
	} # end of if (type.response == "binary")
}

# functions for performance measures
funconc = function(time, status, score, more = F) {
  out = summary(coxph(Surv(time, status) ~ score))$concor # concordance
  #           survConcordance(Surv(time, status) ~ score)$concordance
  if (more) out else out[1] # without se
}
funlogrank = function(time, status, score, more = F) {
  out = summary(coxph(Surv(time, status) ~ score))$sctest # logrank test
  if (more) out[-2] else out[1] # without p-value			
}

find.cutoff = function(pred, label, time = NULL, type.measure.cutoff = c("mindist", "harmonic", "antiharmonic", "balanced", "fscore", "mcc", "median", "0.5", "logrank", "LR", "concordance"), best.only = TRUE) {
	type.measure.cutoff = match.arg(type.measure.cutoff)
	#if (length(unique(pred))==1) pred = pred+rnorm(length(pred),0,.Machine$double.eps*10) 
	if (type.measure.cutoff == "median") {
		cutoff = median(pred)
	} else if (type.measure.cutoff == "0.5") {	
		cutoff = 0.5	
	} else {
		if (length(unique(pred)) > 1) {
			if (type.measure.cutoff == "mindist") {
				perf1 = performance(prediction(pred, label), 'sens')
				perf2 = performance(prediction(pred, label), 'spec')
				th = perf1@x.values[[1]]
				meas = -sqrt((1-perf1@y.values[[1]])^2+(1-perf2@y.values[[1]])^2)
			} 	
			if (type.measure.cutoff == "harmonic") {
				perf1 = performance(prediction(pred, label), 'sens')
				perf2 = performance(prediction(pred, label), 'spec')
				th = perf1@x.values[[1]]
				meas = 2*perf1@y.values[[1]]*perf2@y.values[[1]]/(perf1@y.values[[1]]+perf2@y.values[[1]])
			} 
			if (type.measure.cutoff == "antiharmonic") {
				perf1 = performance(prediction(pred, label), 'sens')
				perf2 = performance(prediction(pred, label), 'spec')
				th = perf1@x.values[[1]]
				meas = perf1@y.values[[1]]*perf2@y.values[[1]]/(2-perf1@y.values[[1]]*perf2@y.values[[1]])
			} 	
			if (type.measure.cutoff == "balanced") {
				perf1 = performance(prediction(pred, label), 'sens')
				perf2 = performance(prediction(pred, label), 'spec')
				th = perf1@x.values[[1]]
				meas = (perf1@y.values[[1]]+perf2@y.values[[1]])/2
			} 	
			if (type.measure.cutoff == "fscore") {
				perf1 = performance(prediction(pred, label), 'f')
				th = perf1@x.values[[1]]
				meas = perf1@y.values[[1]]
			} 	
			if (type.measure.cutoff == "mcc") {
				perf1 = performance(prediction(pred, label), 'mat')
        th = perf1@x.values[[1]]
				meas = perf1@y.values[[1]]
			}	
      if (type.measure.cutoff == "logrank") {
        if (is.null(time)) stop('Need a time vector!')
        th = sort(pred)       
        meas = sapply(th, function(x) {summary(coxph(Surv(time, label) ~ I(pred > x)))$sctest["test"]}) 
      }
			if (type.measure.cutoff == "LR") {
			  if (is.null(time)) stop('Need a time vector!')
			  th = sort(pred)       
			  meas = sapply(th, function(x) {summary(coxph(Surv(time, label) ~ I(pred > x)))$logtest["test"]}) 
			}
			if (type.measure.cutoff == "concordance") {
			  if (is.null(time)) stop('Need a time vector!')
			  th = sort(pred)       
			  meas2 = t(sapply(th, function(x) {summary(coxph(Surv(time, label) ~ I(pred > x)))$concordance})) 
        colnames(meas2) = c("performance", "se")
        meas = meas2[, "performance"]
			}
			maxmeas = max(meas)
			maxmeas.ind = which.max(meas) # only the first
			
# 			if (maxmeas.ind < length(meas)) {
# 				incre = perf1@x.values[[1]][maxmeas.ind] - perf1@x.values[[1]][maxmeas.ind+1]
# 			} else {
# 				incre = perf1@x.values[[1]][maxmeas.ind-1] - perf1@x.values[[1]][maxmeas.ind]
# 			} 
 			cutoff = th[maxmeas.ind]
# 			cutoff = perf1@x.values[[1]][maxmeas.ind] - incre/2 # use midpoint for continuity (?) correction. 
# 			If not check the condition 'maxmeas.ind < length(meas)', have "Warning messages: In max(ind) : no non-missing arguments to max; returning -Inf" in with.give.cutoff(), because cutoff = NA
		} else cutoff = unique(pred) 
	}
	if (best.only || type.measure.cutoff %in% c("median", "0.5")) cutoff 
  else if(type.measure.cutoff == "concordance") list(Best = cutoff, All = cbind(cutoff = th, meas2))
  else list(Best = cutoff, All = cbind(cutoff = th, performance = meas))
	
	#list(cutoff = cutoff, measure = meas, maxmeasure = maxmeas, which = maxmeas.ind)
}
with.given.cutoff = function(obj, label = NULL, type.measure = 'sens', type.measure.cutoff = "balanced", is.cv = TRUE) {
	if (is.cv) {
		cutoff.list = Map(function(x, y) find.cutoff(x, y, type.measure.cutoff = type.measure.cutoff), obj$pred.train.list, obj$label.train.list)
		perf.list = Map(function(cutoff, pred, label) {
			perf = performance(prediction(pred, label), type.measure)
			ind = which(perf@x.values[[1]]>=cutoff) 
			result = perf@y.values[[1]][max(ind)] 
			if (result=="NaN") result = 1
			result
		}, cutoff.list, obj$pred.test.list, obj$label.test.list)
		
		perf.vec = unlist(perf.list)
		cvm = weighted.mean(perf.vec, w = unlist(obj$nobs.test.list))
		cvsd = sqrt(weighted.mean((perf.vec-cvm)^2, w = unlist(obj$nobs.test.list))/(length(perf.vec)-1))
	
		list(performance = perf.vec, measure = type.measure, mean = cvm, se = cvsd)
	} else {
		cutoff = find.cutoff(obj$fitted, label, type.measure.cutoff = type.measure.cutoff)				
		perf = performance(prediction(obj$fitted, label), type.measure)
		ind = which(perf@x.values[[1]]>=cutoff) 
		result = perf@y.values[[1]][max(ind)] 
		if (result=="NaN") result = 1

		result
	}
}
######## Perform bootstrap analysis
run.boot = function(x, y, family = "cox", standardize = FALSE, stratify = T, B=1, maxit = 100000) {
	epsilon = .Machine$double.eps*10
	colnames(x) = 1:ncol(x)
	lambda.lasso = rep(NA, B)
	lambda.adalasso = rep(NA, B)
	which.adalasso = vector("list", B)
	coef.adalasso = vector("list", B)
	for (b in seq_len(B)) {
	repeat{ 
    #proportional (stratified) sampling
    if (family == "binomial" && stratify)	{
      ind1.orig = which(y==unique(y)[1])
      ind2.orig = which(y==unique(y)[2])
      ind1 = ind1.orig[sample.int(length(ind1.orig), replace = T)]
      ind2 = ind2.orig[sample.int(length(ind2.orig), replace = T)]
      ind = c(ind1, ind2)
    } else ind = sample(1:nrow(x), replace = T)
		xb = x[ind, , drop = F]; if (is.vector(y)) yb = y[ind] else yb = y[ind, , drop = F]
		fit = try(my.adalasso.1dim(x = xb, y = yb, family = family, standardize = standardize, foldid.list = list(seq(nrow(xb))), maxit = maxit), silent = F)
    # for error checking     
    if (!inherits(fit, "try-error")) {
      break # Error in y[, c(2, 1)] : incorrect number of dimensions: this happened when proportional sampling was not used (e.g., with only one positive case, we have only negative cases in a training fold)
    }
	}
		cat(paste("b = ", b,"  of ", B, "\n", sep = ""))
		if (length(fit$lambda.lasso)) lambda.lasso[b] = fit$lambda.lasso
		if (length(fit$lambda.adalasso)) lambda.adalasso[b] = fit$lambda.adalasso
		if (length(fit$coefficients)) {
			which.adalasso[[b]] = as.numeric(names(fit$coefficients))
			coef.adalasso[[b]] = fit$coefficients
		}
	}
	list(lambda.lasso = lambda.lasso, lambda.adalasso = lambda.adalasso, which.adalasso = which.adalasso, coef.adalasso = coef.adalasso)	
}

######## Plot dendrogram

my.dendro = function(boot, freq.th = 50, varnames, names.finalcoef, horiz = TRUE, plot.name, longwidth = 6, shortwidth = 4.5, horizmar = c(4,1,1,14)+.1, cex = 1
, B = 500, grid = T, height = T, col.pos = "orange", col.neg = "blue", reverse = F, plotmath = F) {
	DISTmat = sort.freq.pair(freq.pair(boot$which.adalasso, varnames = varnames))
	hc = function.hclust(DISTmat, th = freq.th, hang = -3, xlab = NA, sub = NA)
	sgn.vec = signs(boot$which.adalasso, boot$coef.adalasso, varnames = varnames)$max.sgn
	sgn.vec = ifelse(sgn.vec==1,"+", "-")
	freq.vec = diag(DISTmat)
	colLab <- function(n) {
	  	if (is.leaf(n)) {
		    a <- attributes(n)
		    if (!plotmath) {
          attr(n, "label") <- paste(a$label, " (", freq.vec[names(freq.vec)==a$label], "/", B, ")", sep = "") # change the node label  
		    } else {
          attr(n, "label") <- parse(text = paste(a$label, "~(", freq.vec[names(freq.vec)==a$label], "/", B, ")", sep = "")) # change the node label    
		    }
        #if (a$label %in% names.finalcoef) my.lab.col = "black" else my.lab.col = gray(.6)
        #: same as 
        my.lab.col = c(grey(.6), "black")[(a$label %in% names.finalcoef)+1]        
		    #if (sgn.vec[names(sgn.vec)==a$label]=="+") my.col = col.pos else my.col = col.neg
        #: same as 
        my.col = c(col.neg, col.pos)[(sgn.vec[names(sgn.vec)==a$label]=="+")+1]
			attr(n, "nodePar") <- c(a$nodePar, list(lab.col = my.lab.col, col = my.col, pch = 19, cex = 1.2)) # change the node color	
		} 
		n
	}
	clusDendro = as.dendrogram(hc)
	clusDendro = dendrapply(clusDendro, colLab)
  if (reverse) clusDendro = rev(clusDendro)

	pdf(paste(output.folder, "plot_bootstrap_dendr", freq.th, "_horiz", ifelse(horiz, "T", "F"), "_", plot.name, ".pdf", sep = ""), height = ifelse(horiz, shortwidth, longwidth), width = ifelse(horiz, longwidth, shortwidth))
	vertmar = horizmar[c(4,1,2,3)]
	if (horiz) plotmar = horizmar else plotmar = vertmar
	op <- par(mar = plotmar, cex = cex) #par(mar = par("mar") + c(9,0,-3,-1)) 
  plot(clusDendro, horiz = horiz, cex = cex, axes = height)
	if (grid) {if (horiz) abline(v = seq(0,500, 20), lty = 3, col = "grey") else abline(h = seq(0,500, 20), lty = 3, col = "grey")}
	par(op)
	dev.off()
}

signs = function(wlist, clist, varnames = NULL, use.allvar = FALSE, nallvar = length(varnames)) {
	wvec = unlist(wlist)
	if (use.allvar) {
		wvec = factor(wvec, levels = seq(nallvar))
	}
	Sgn = table(factor(sign(unlist(clist)), levels = c(-1,1)), wvec)
	names(dimnames(Sgn))[2] = "" # remove name 'wvec'
	if (!is.null(varnames)) {
		newnames = varnames[as.numeric(colnames(Sgn))]
		colnames(Sgn) = newnames
	}
	max.sgn = (Sgn[1,]<Sgn[2,])*2-1
	max.count = pmax(Sgn[1,], Sgn[2,])
	list(table = Sgn, max.sgn = max.sgn, max.count = max.count)
}

# frequencies (pairwise)
freq.pair = function(wlist, varnames = NULL, use.allvar = FALSE, nallvar = length(varnames)) {
	if (use.allvar) allvar = seq(nallvar) else allvar = sort(unique(unlist(wlist))) 
	freqlist = lapply(wlist, function(vec) outer(allvar, allvar, FUN = function(x, y) x %in% vec & y %in% vec))
	freq.pair = Reduce("+", freqlist)
	if (!is.null(varnames)) {
		rownames(freq.pair) = varnames[allvar]
		colnames(freq.pair) = varnames[allvar]
	} 
	freq.pair
}

sort.freq.pair = function(mat) {
	mat[order(-diag(mat)), order(-diag(mat))]
}

function.hclust = function(mat, th=50, method = "complete", plot.it = FALSE, ...) {
	ind = diag(mat)>=th
	hc = hclust(d = as.dist(mat[ind, ind]), method = method)
	if (plot.it) plot(hc, ... )
	return(hc)
}

######## Fit adaptive LASSO
#### 'my.adalasso.1dim' fits adaptive lasso, with best lambda chosen by cross-validation using the deviance as performance measure (through 'cv.glmnet'). Before 'cv.glmnet' is called, foldid is generated by 'my.cv.glmnet' (option 'foldid.only = TRUE') and used as input for 'cv.glmnet'. 
my.adalasso.1dim = function (x, y, family = c("gaussian", "binomial", "cox"), foldid.list = NULL, nfolds = 10, nrepeat = 1, standardize = FALSE, adaptive = TRUE, alpha.wt = 1, grouped = TRUE, stratify.folds = TRUE, subsample.train = FALSE, maxit = 10000) 
{
	epsilon = .Machine$double.eps*10
	nobs = nrow(x)
	family = match.arg(family)
	if (family != "binomial" && family != "cox") {stratify.folds = FALSE; subsample.train = FALSE}
    ### find (representative) seq of lambdas
    lambdas = glmnet(x, y, family = family, standardize = standardize, alpha = alpha.wt, maxit = maxit)$lambda
    nlamb = length(lambdas)
    
    #### foldid.list
	if (is.null(foldid.list)) {
		foldid.list = vector("list", nrepeat)
		for (rr in seq_len(nrepeat)) {
			foldid.list[[rr]] = my.cv.glmnet(x = x, y = y, family =family, foldid.only = TRUE, nfolds = nfolds, stratify.folds = stratify.folds)
		}
	} else {
		nrepeat = length(foldid.list)
		nfolds = max(foldid.list[[1]])
	}
    
	cvmvec = rep(0, length = nlamb)
	for (rr in seq_len(nrepeat)) {
		foldid = foldid.list[[rr]]
		cvmvec = cvmvec + cv.glmnet(x=x, y=y, family = family, lambda = lambdas, standardize = standardize, alpha = alpha.wt, foldid = foldid, grouped = grouped, maxit = maxit)$cvm
	}
	cvm.lasso = cvmvec/nrepeat
	cvmin.lasso = min(cvm.lasso)
	lambda.lasso.min = lambdas[cvm.lasso == cvmin.lasso]
	# final lasso model
	fit = glmnet(x, y, family = family, standardize = standardize, alpha = alpha.wt, lambda = lambda.lasso.min[1], maxit= maxit)
	# type = c("link", "response", "coefficients", "nonzero", "class")
	if (family == "cox") {
		fitted.lasso = as.vector(predict(fit, newx = x, type = "link")) # linear predictor or log risk score, same as log(predict(fit, ..., type = "response")
	} else {
		fitted.lasso = as.vector(predict(fit, newx = x, type = "response")) # probabilities if family == "binomial" 
	}
	coef.lasso = fit$beta[,1]
	print(coef.lasso)	
	coef.lasso[abs(coef.lasso)<epsilon] = 0
	intercept.lasso = fit$a0
	cvm.ada = NULL
	lambdas.ada = NULL
    cvmin.ada = NULL
    lambda.ada.min = NULL
    weights = NULL
    intercept.ada = NULL 
    coef.ada = NULL
	if (adaptive && !all(coef.lasso==0)) {
		which0 = as.vector(which(coef.lasso==0))
		weights = 1/pmax(abs(coef.lasso), epsilon)
		weights[which0] = 0
		lambdas.ada = glmnet(x, y, family = family, standardize = standardize, penalty.factor = weights, exclude = which0, alpha = 1, maxit = maxit)$lambda
		cvmvec = rep(0, length = length(lambdas.ada))
		for (rr in seq_len(nrepeat)) {
			foldid = foldid.list[[rr]]
			cvmvec = cvmvec + cv.glmnet(x=x, y=y, family = family, lambda = lambdas.ada, standardize = standardize, alpha = 1, foldid = foldid, maxit = maxit, penalty.factor = weights, exclude = which0, grouped = grouped)$cvm
		}
		cvm.ada = cvmvec/nrepeat
		cvmin.ada = min(cvm.ada)
		lambda.ada.min = lambdas.ada[cvm.ada == cvmin.ada]
		#final adaptive lasso model
    	fit.ada = glmnet(x, y, standardize = standardize, family = family, alpha = 1, penalty.factor = weights, exclude = which0, lambda = lambda.ada.min[1], maxit = maxit)
    	coef.ada = fit.ada$beta[,1]
    	intercept.ada = fit.ada$a0
    	# type = c("link", "response", "coefficients", "nonzero", "class")
    	# if type = "response", fitted relative-risk for family == "cox", fitted probabilities for family == "binomial"
    	# use type = "link" for linear predictors
    	if (family == "cox") {
    		fitted.adalasso = as.vector(predict(fit.ada, newx = x, type = "link")) # linear predictor or log risk score, same as log(predict(fit, ..., type = "response")
    	} else {
    		fitted.adalasso = as.vector(predict(fit.ada, newx = x, type = "response")) # probabilities if family == "binomial"
    	}
	} else {
		# if (family == "cox") fitted.adalasso = rep(0, nrow(x))
		# if (family == "binomial") fitted.adalasso = rep(.5, nrow(x)) 
		intercept.ada = intercept.lasso 
		fitted.adalasso = fitted.lasso
	} 
#  print(fitted.adalasso)  
    return(list(xnames = colnames(x), cv = data.frame(lambda.ada = lambdas.ada, cvm.ada = cvm.ada),
	foldid.list = foldid.list,
    cvmin.lasso = cvmin.lasso, 
    lambda.lasso = lambda.lasso.min, 
    intercept.lasso = intercept.lasso, 
    coefficients.lasso = coef.lasso[coef.lasso != 0],
    fitted.lasso = fitted.lasso,
    cvmin.adalasso = cvmin.ada, 
    lambda.adalasso = lambda.ada.min, 
    weights = weights[weights != 0],
    intercept = intercept.ada, 
    coefficients = coef.ada[coef.ada != 0], 
    fitted = fitted.adalasso))
}
my.cv.glmnet = function(glmnet.object.list = NULL, x, y, family = c("gaussian", "binomial", "cox"), foldid, lambda = NULL, foldid.only = FALSE, nfolds = 10, standardize = TRUE, alpha = 1, stratify.folds = FALSE, subsample.train = FALSE, grouped = T, maxit = 1000)
{
	family = match.arg(family)
	lambdas = lambda
	nobs = nrow(x)
	if (!missing(foldid)) {
		k = max(foldid)
		all.folds = split(seq(nobs), foldid)
		counts.infold = sapply(all.folds, length)	
	} else {
		if (family != "binomial" && family != "cox") {
			stratified.folds = FALSE; subsample.train = FALSE
		}
		k = nfolds

		if (stratify.folds) {	
			if (family == "cox") {
				status = y[, "status"]	
			} else {
				status = y
			}
			class.counts = table(status)
			all.class = names(class.counts)
		 	smallest.class = all.class[class.counts == min(class.counts)]
		 	other.class = all.class[class.counts > min(class.counts)]
		 	mode(all.class) = mode(status)
		 	mode(smallest.class) = mode(status)
		 	mode(other.class) = mode(status)
			all.folds = vector("list", k)
			counts.infold = rep(0, k)
			for (cls in seq_along(all.class)) {
				repeat {
					permu = sample(rep(sample(1:k), length = class.counts[cls])) 
					# Tried `rep(sample(1:k), length = class.counts[cls])' 
					# but it produces a pattern, e.g., for k=10, 11,21,31,41,... always in the same fold (so, repeated CV not efficient)
					# -> avoid this by permuting the vector again
					if (class.counts[cls] < k) permu = factor(permu, levels = 1:k) # to avoid losing folds with zero count in table(permu)
					counts.infold.temp = counts.infold + table(permu)
					if (max(counts.infold.temp)-min(counts.infold.temp) <= 1) 
					{
						counts.infold = counts.infold.temp
						break
					}
				}
				all.folds = Map(c, all.folds, split(which(status == all.class[cls]), permu)) # list of length k
			}
		} else {
			all.folds = split(sample(seq(nobs)), rep(1:k, length = nobs)) # size of fold i >= size of fold j when i<j
			counts.infold = sapply(all.folds, length)
		}
		
		# optional: vector of fold id
		foldid = rep(1:k, counts.infold)[order(unlist(all.folds))]; if (any(is.na(foldid))) {print(all.folds); print(all.class); print(permu)}
	} # end of if (!missing(foldid)) ... else
		
	if (foldid.only) return(foldid)	else {
	evalmat = matrix(NA, k, length(lambdas))
	predmat = matrix(NA, nrow(x), length(lambdas))
	for (i in seq_len(k)) {
		trainfit = glmnet.object.list[[i]]		
		omit = all.folds[[i]]
		if (family == "cox") xtrain = x[-omit, , drop = FALSE]
		
		xtest = x[omit, , drop = FALSE]

		if (family == "cox") {	
			ytrain = y[-omit,]; ytest = y[omit,] 
		} else { ytest = y[omit] }

		if (family != "cox") 
			pred = predict(trainfit, newx = xtest, type = "response", s = lambdas)	
		if (family == "cox") {
		 	pred = NA
			coefmat = predict(trainfit, type = "coeff", s = lambdas)
		}
		if (length(omit) == 1) pred  = matrix(pred, nrow = 1)
		# compute deviance
		if (family == "gaussian") {
				evalmat[i,] = apply((ytest-pred)^2, 2, mean)
			}
		if (family == "binomial") {
				evalmat[i,] = apply(ytest*log(pred)+(1-ytest)*log(1-pred), 2, mean) * (-2)
			}
		if (family == "cox") {
			# use function 'coxnet.deviance' in package 'glmnet'
				plfull = coxnet.deviance(x=x, y=y, beta=coefmat)
      			plminusk= coxnet.deviance(x=xtrain, y=ytrain, beta=coefmat)
				evalmat[i,seq(along=plfull)] = plfull-plminusk
			}
		predmat[omit, ] = pred
	} 
	
	if (grouped) {
		cvraw = evalmat
		wts = counts.infold
		if (family == "cox")  { 
			status=y[,"status"]
  			wts=as.vector(tapply(y[,"status"],foldid,sum))
  			cvraw=cvraw/wts
   		}
	} else {
		cvraw = switch(family, "gaussian" = (y-predmat)^2, "binomial" = -2 * (y*log(predmat)+(1-y)*log(1-predmat)))
		wts = rep(1, nrow(cvraw))
	} 
	cvm = apply(cvraw, 2, weighted.mean, w = wts)#, na.rm = T) # same for grouped = TRUE or FALSE
	cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, w = wts)/(nrow(cvraw)-1)) 
	cvmin = min(cvm)
	getmin = getmin(lambdas,cvm,cvsd) # from glmnet
	lambda.min = getmin$lambda.min
	lambda.1se = getmin$lambda.1se
	# lambda.min = lambdas[cvm == cvm.min]
	# lambda.1se = max(lambdas[cvm <= cvm.min + cvsd])
	fit.best = glmnet(x, y, family = family, standardize = standardize, lambda = lambdas, maxit = maxit)
	nzero = sapply(predict(fit.best, type = "nonzero"), length)
	list(lambda = lambdas, cvm = cvm, cvsd = cvsd, cvmin = cvmin, nzero = nzero, lambda.min = lambda.min, lambda.1se = lambda.1se
	 , all.folds = all.folds, foldid = foldid#, dist = dist, dist.train = dist.train
	)
	} # end of if (foldid.only) ... else
}
# copied from glmnet getmin.R
getmin=function(lambda,cvm,cvsd){
  cvmin=min(cvm,na.rm=TRUE)
  idmin=cvm<=cvmin
  lambda.min=max(lambda[idmin],na.rm=TRUE)
  idmin=match(lambda.min,lambda)
  semin=(cvm+cvsd)[idmin]
  idmin=cvm<=semin
  lambda.1se=max(lambda[idmin],na.rm=TRUE)
  list(lambda.min=lambda.min,lambda.1se=lambda.1se)
}

######## Fit non-regularized (in fact minimally-regularized) Cox model ('no penalty')
my.cvm.nopenalty = function (x, y, family = c("gaussian", "binomial", "cox"), foldid.list = NULL, nfolds = 10, nrepeat = 1, standardize = FALSE, adaptive = TRUE, alpha.wt = 1, grouped = TRUE, stratify.folds = TRUE, subsample.train = FALSE, maxit = 10000) {
	epsilon = .Machine$double.eps*10
	nobs = nrow(x)
	family = match.arg(family)
	if (family != "binomial" && family != "cox") {stratify.folds = FALSE; subsample.train = FALSE}
    ### find (representative) seq of lambdas
    lambdas = c(glmnet(x, y, family = family, standardize = standardize, alpha = alpha.wt, maxit = maxit*10)$lambda,0)
    nlamb = length(lambdas)
    
    #### foldid.list
	if (is.null(foldid.list)) {
		foldid.list = vector("list", nrepeat)
		for (rr in seq_len(nrepeat)) {
			foldid.list[[rr]] = my.cv.glmnet(x = x, y = y, family =family, foldid.only = TRUE, nfolds = nfolds, stratify.folds = stratify.folds, maxit = maxit*10)
		}
	} else {
		nrepeat = length(foldid.list)
		nfolds = max(foldid.list[[1]])
	}
    
	cvmvec = rep(0, length = nlamb)
	for (rr in seq_len(nrepeat)) {print(rr)
		foldid = foldid.list[[rr]]; print(foldid)
		cvmvec = cvmvec + cv.glmnet(x=x, y=y, family = family, lambda = lambdas, standardize = standardize, alpha = alpha.wt, foldid = foldid, maxit = maxit)$cvm
	}
	cvm.lasso = cvmvec/nrepeat
	fit = glmnet(x=x, y=y, family = family, lambda = 0, standardize = standardize, alpha = alpha.wt, maxit = maxit*10)
	coef = fit$beta[,1]
	intercept = fit$a0
	# type = c("link", "response", "coefficients", "nonzero", "class")
	if (family == "cox") {
		fitted = predict(fit, newx = x, type = "link") # linear predictor or log risk score, same as log(predict(fit, ..., type = "response")
	} else {
		fitted = predict(fit, newx = x, type = "response") # probabilities if family == "binomial"
	}
    return(list(cvm = cvm.lasso, intercept = intercept, coefficients = coef[coef!=0], fitted = fitted))
}
