######################################################################
### Q-based statistics
######################################################################

qlss <- function(...) UseMethod("qlss")

qlss.numeric <- function(x, probs = 0.1, ...){

if(any(!probs > 0) | any(!probs < 0.5)) stop("Quantile index out of range: probs must be > 0 and < 0.5")
if(any(duplicated(probs))) probs <- unique(probs)
call <- match.call()

nq <- length(probs)

vec3 <- as.numeric(do.call(what = quantile, args = list(x = x, probs = 1:3/4, ...)))
Me <- vec3[2]
IQR <- vec3[3] - vec3[1]

IPR <- Ap <- Tp <- rep(NA, nq)
for(i in 1:nq){
	vecp <- as.numeric(do.call(what = quantile, args = list(x = x, probs = probs[i], ...)))
	vecq <- as.numeric(do.call(what = quantile, args = list(x = x, probs = 1 - probs[i], ...)))

	IPR[i] <- vecq - vecp
	Ap[i] <- (vecq - 2*Me + vecp)/IPR[i]
	Tp[i] <- IPR[i]/IQR
}
names(IPR) <- names(Ap) <- names(Tp) <- probs

val <- list(location = list(median = Me), scale = list(IQR = IQR, IPR = IPR), shape = list(skewness = Ap, shape = Tp))
val$probs <- probs
val$call <- call
class(val) <- "qlss"
return(val)
}

qlss.default <- function(fun = "qnorm", probs = 0.1, ...){

if(any(!probs > 0) | any(!probs < 0.5)) stop("Quantile index out of range: probs must be > 0 and < 0.5")
if(any(duplicated(probs))) probs <- unique(probs)
call <- match.call()

nq <- length(probs)

vec3 <- as.numeric(do.call(what = match.fun(fun), args = list(p = 1:3/4, ...)))
Me <- vec3[2]
IQR <- vec3[3] - vec3[1]

IPR <- Ap <- Tp <- rep(NA, nq)
for(i in 1:nq){
	vecp <- as.numeric(do.call(what = match.fun(fun), args = list(p = probs[i], ...)))
	vecq <- as.numeric(do.call(what = match.fun(fun), args = list(p = 1 - probs[i], ...)))

	IPR[i] <- vecq - vecp
	Ap[i] <- (vecq - 2*Me + vecp)/IPR[i]
	Tp[i] <- IPR[i]/IQR
}
names(IPR) <- names(Ap) <- names(Tp) <- probs

val <- list(location = list(median = Me), scale = list(IQR = IQR, IPR = IPR), shape = list(skewness = Ap, shape = Tp))
val$probs <- probs
val$call <- call
class(val) <- "qlss"
return(val)
}

qlss.formula <- function(formula, probs = 0.1, data = sys.frame(sys.parent()), subset, weights, na.action, contrasts = NULL, method = "fn", type = "rq", tsf = "mcjI", symm = TRUE, dbounded = FALSE, lambda = NULL, conditional = FALSE, ...){

if(any(!probs > 0) | any(!probs < 0.5)) stop("Quantile index out of range: probs must be > 0 and < 0.5")
if(any(duplicated(probs))) probs <- unique(probs)
if(type == "rqt" && tsf == "mcjII") stop("'mcjII' not available for qlss")
if (!inherits(formula, "formula") || length(formula) != 3) {
	stop("\nFixed-effects model must be a formula of the form \"y ~ x\"")
}

call <- match.call()
mc <- match.call(expand.dots = TRUE)
fitLs <- names(mc)

mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("get_all_vars")
mf <- eval(mf, parent.frame())
if (method == "model.frame") 
        return(mf)

taus <- c(1:3/4, probs, 1-probs)
nq <- length(probs)

val <- list()
if(type == "rq"){
	fitLs <- as.list(mc[fitLs %in% names(formals(rq))])
	fitLs$data <- mf
	fitLs$tau <- sort(taus) # rq automatically sorts tau, just enforcing here
	val$fit <- do.call(rq, args = fitLs)
}

if(type == "rqt"){
	if(conditional & length(lambda) < length(taus)) stop(paste0("Length of 'lambda' must be ", length(taus), ". See details in '?qlss'"))
	fitLs <- as.list(mc[fitLs %in% names(formals(tsrq))])
	fitLs$data <- mf
	fitLs$tau <- sort(taus) # rqt does not sort, but rq does, so...
	val$fit <- do.call(tsrq, args = fitLs)
}

# sort wrt taus
ii <- match(taus,fitLs$tau)
Fitted <- val$fit$fitted.values[,ii]

vec3 <- Fitted[,1:3]
vecp <- Fitted[ , -c(1:3), drop = FALSE][ , 1:nq, drop = FALSE]
vecq <- Fitted[ , -c(1:3), drop = FALSE][ , (nq+1):(2*nq), drop = FALSE]
Me <- drop(Fitted[,2])
IQR <- drop(vec3[,3] - vec3[,1])
IPR <- vecq - vecp
Ap <- (vecq - 2 * Me + vecp)/IPR
Tp <- IPR/IQR
colnames(IPR) <- colnames(Ap) <- colnames(Tp) <- probs

val$location <- list(median = matrix(Me))
val$scale <- list(IQR = matrix(IQR), IPR = IPR)
val$shape <- list(skewness = Ap, shape = Tp)
val$call <- call
val$probs <- probs
val$type <- type
val$method <- method

class(val) <- "qlss"
return(val)

}

predict.qlss <- function(object, newdata, interval = FALSE, level = 0.95, R = 200, na.action = na.pass, trim = 0.05, ...){

type <- object$type
if(!length(type)) return(object)

val <- switch(type,
	rq = qlssPredRq(object = object, newdata = newdata, interval = interval, level = level, R = R, na.action = na.action, trim = trim),
	rqt = qlssPredRqt(object = object, newdata = newdata, interval = interval, level = level, R = R, na.action = na.action, trim = trim)
)

val$probs <- object$probs

class(val) <- "qlss"
return(val)

}

qlssPredRq <- function(object, newdata, interval, level, R, na.action, trim){

probs <- object$probs
nq <- length(probs)

fit <- object$fit
taus <- c(1:3/4, probs, 1-probs)

if(missing(newdata))
	{x <- fit$x}	else {
	mt <- terms(fit)
	Terms <- delete.response(mt)
	mf <- model.frame(Terms, newdata, na.action = na.action, xlev = fit$xlevels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, mf)
	x <- model.matrix(Terms, mf, contrasts.arg = fit$contrasts)
}

# sort wrt taus
ii <- match(taus, fit$tau)
Fitted <- x%*%fit$coefficients
Fitted <- Fitted[,ii]

vecp <- Fitted[ , -c(1:3), drop = FALSE][ , 1:nq, drop = FALSE]
vecq <- Fitted[ , -c(1:3), drop = FALSE][ , (nq+1):(2*nq), drop = FALSE]
Me <- drop(Fitted[,2])
IQR <- drop(Fitted[,3] - Fitted[,1])
IPR <- vecq - vecp
Ap <- (vecq - 2 * Me + vecp)/IPR
Tp <- IPR/IQR
colnames(IPR) <- colnames(Ap) <- colnames(Tp) <- probs
val <- list(location = list(median = matrix(Me)), scale = list(IQR = matrix(IQR), IPR = IPR), shape = list(skewness = Ap, shape = Tp))

if(interval){
	fit.s <- summary(fit, se = "boot", R = R, covariance = TRUE)
	B <- lapply(fit.s, function(x) x$B)
	
	Fitted <- lapply(B, function(b,x) x%*%t(b), x = x)
	Fitted <- Fitted[ii] # sort wrt to taus
	
	sdtrim <- function(u, trim){
		sel1 <- u >= quantile(u, probs = trim/2, na.rm = TRUE)
		sel2 <- u <= quantile(u, probs = 1-trim/2, na.rm = TRUE)
		sd(u[sel1 & sel2])
	}
	
	level <- level + (1-level)/2

	vecp <- Fitted[-c(1:3)][1:nq]
	vecq <- Fitted[-c(1:3)][(nq+1):(2*nq)]

	tmp <- qt(level, R - 1) * apply(Fitted[[2]], 1, sdtrim, trim = trim)
	Me.ci <- cbind(Me - tmp, Me + tmp)
	tmp <- qt(level, R - 1) * apply(Fitted[[3]] - Fitted[[1]], 1, sdtrim, trim = trim)
	IQR.ci <- cbind(IQR - tmp, IQR + tmp)
	
	IPR.ci <- Ap.ci <- Tp.ci <- list()
	for(j in 1:nq){
		ipr <- vecq[[j]] - vecp[[j]]
		tmp <- qt(level, R - 1) * apply(ipr, 1, sdtrim, trim = trim)
		IPR.ci[[j]] <- cbind(IPR[,j] - tmp, IPR[,j] + tmp)
		ap <- (vecq[[j]] - 2 * Fitted[[2]] + vecp[[j]])/ipr
		tmp <- qt(level, R - 1) * apply(ap, 1, sdtrim, trim = trim)
		Ap.ci[[j]] <- cbind(Ap[,j] - tmp, Ap[,j] + tmp)
		tp <- ipr/(Fitted[[3]] - Fitted[[1]])
		tmp <- qt(level, R - 1) * apply(tp, 1, sdtrim, trim = trim)
		Tp.ci[[j]] <- cbind(Tp[,j] - tmp, Tp[,j] + tmp)
	}
	
	CI <- list(Me = Me.ci, IQR = IQR.ci, IPR = IPR.ci, Ap = Ap.ci, Tp = Tp.ci)
	names(CI$IPR) <- names(CI$Ap) <- names(CI$Tp) <- probs
	val$CI <- CI
}

return(val)
}

qlssPredRqt <- function(object, newdata, interval, level, R, na.action, trim){

probs <- object$probs
nq <- length(probs)

fit <- object$fit
taus <- c(1:3/4, probs, 1-probs)
tsf <- fit$tsf
symm <- attributes(tsf)$symm
dbounded <- attributes(tsf)$dbounded
isBounded <- attributes(tsf)$isBounded
#conditional <- attributes(tsf)$conditional
conditional <- TRUE

if(missing(newdata))
	{x <- fit$x}	else {
	mt <- terms(fit)
	Terms <- delete.response(mt)
	mf <- model.frame(Terms, newdata, na.action = na.action, xlev = fit$xlevels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, mf)
	x <- model.matrix(Terms, mf, contrasts.arg = fit$contrasts)
}

n <- nrow(x)

lambdahat <- fit$lambda
linpred <- x%*%fit$coefficients
Fitted <- matrix(NA, nrow(linpred), ncol(linpred))
for (j in 1:length(taus)) {
	Fitted[, j] <- switch(tsf, mcjI = invmcjI(linpred[, 
		j], lambdahat[j], symm, dbounded), bc = invbc(linpred[, 
		j], lambdahat[j]), ao = invao(linpred[, j], lambdahat[j], 
		symm))
}

if (isBounded) {
	Fitted <- apply(Fitted, 2, function(x, x.r) invmap(x, 
		x.r), x.r = range(fit$y))
}

# sort wrt taus
ii <- match(taus, fit$tau)
Fitted <- Fitted[,ii]

vecp <- Fitted[ , -c(1:3), drop = FALSE][ , 1:nq, drop = FALSE]
vecq <- Fitted[ , -c(1:3), drop = FALSE][ , (nq+1):(2*nq), drop = FALSE]
Me <- drop(Fitted[,2])
IQR <- drop(Fitted[,3] - Fitted[,1])
IPR <- vecq - vecp
Ap <- (vecq - 2 * Me + vecp)/IPR
Tp <- IPR/IQR
colnames(IPR) <- colnames(Ap) <- colnames(Tp) <- probs
val <- list(location = list(median = matrix(Me)), scale = list(IQR = matrix(IQR), IPR = IPR), shape = list(skewness = Ap, shape = Tp))

if(interval){
	B <- summary(fit, se = "boot", R = R, conditional = conditional, covariance = TRUE)$B

	linpred <- lapply(B, function(b,x) x%*%t(b), x = x)
	Fitted <- list()
	for (j in 1:length(taus)) {
		Fitted[[j]] <- switch(tsf,
		mcjI = apply(linpred[[j]], 2, function(x,l,s,d) invmcjI(x,l,s,d), l = lambdahat[j], s = symm, d = dbounded),
		bc = apply(linpred[[j]], 2, function(x,l,s,d) invbc(x,l), l = lambdahat[j]),
		ao = apply(linpred[[j]], 2, function(x,l,s) ao(x,l,s), l = lambdahat[j], s = symm)
		)
		if (isBounded) {
			Fitted[[j]] <- apply(Fitted[[j]], 2, function(x, x.r) invmap(x, 
				x.r), x.r = range(fit$y))
		}
	}
	Fitted <- Fitted[ii] # sort wrt to taus
}

if(interval){

	sdtrim <- function(u, trim){
		sel1 <- u >= quantile(u, probs = trim/2, na.rm = TRUE)
		sel2 <- u <= quantile(u, probs = 1-trim/2, na.rm = TRUE)
		sd(u[sel1 & sel2])
	}
	
	level <- level + (1-level)/2

	vecp <- Fitted[-c(1:3)][1:nq]
	vecq <- Fitted[-c(1:3)][(nq+1):(2*nq)]

	tmp <- qt(level, R - 1) * apply(Fitted[[2]], 1, sdtrim, trim = trim)
	Me.ci <- cbind(Me - tmp, Me + tmp)
	tmp <- qt(level, R - 1) * apply(Fitted[[3]] - Fitted[[1]], 1, sdtrim, trim = trim)
	IQR.ci <- cbind(IQR - tmp, IQR + tmp)
	
	IPR.ci <- Ap.ci <- Tp.ci <- list()
	for(j in 1:nq){
		ipr <- vecq[[j]] - vecp[[j]]
		tmp <- qt(level, R - 1) * apply(ipr, 1, sdtrim, trim = trim)
		IPR.ci[[j]] <- cbind(IPR[,j] - tmp, IPR[,j] + tmp)
		ap <- (vecq[[j]] - 2 * Fitted[[2]] + vecp[[j]])/ipr
		tmp <- qt(level, R - 1) * apply(ap, 1, sdtrim, trim = trim)
		Ap.ci[[j]] <- cbind(Ap[,j] - tmp, Ap[,j] + tmp)
		tp <- ipr/(Fitted[[3]] - Fitted[[1]])
		tmp <- qt(level, R - 1) * apply(tp, 1, sdtrim, trim = trim)
		Tp.ci[[j]] <- cbind(Tp[,j] - tmp, Tp[,j] + tmp)
	}
	
	CI <- list(Me = Me.ci, IQR = IQR.ci, IPR = IPR.ci, Ap = Ap.ci, Tp = Tp.ci)
	names(CI$IPR) <- names(CI$Ap) <- names(CI$Tp) <- probs
	val$CI <- CI
}

return(val)
}

plot.qlss <- function(x, z, whichp = NULL, interval = FALSE, type = "l", ...){

probs <- x$probs
if(!is.null(whichp)){
	if(length(whichp) > 1) stop(cat("Only one value for 'whichp' to choose from", probs, "\n"))
	sel <- which(whichp == probs)
	if(!(length(sel))) stop(cat("Choose 'whichp' from", probs, "\n"))
} else {
	sel <- 1
}

r <- order(z)
n <- length(z)
if(interval){
	if(is.null(x$CI)) stop("Use 'predict.qlss' with 'interval = TRUE' first.") else CI <- x$CI
}

if(interval){
	yl1 <- range(c(x$location$median,CI$Me))
	yl2 <- range(c(x$scale$IQR,CI$IQR))
	yl3 <- range(c(x$shape$skewness[,sel],CI$Ap[[sel]]))
	yl4 <- range(c(x$shape$shape[,sel],CI$Tp[[sel]]))
} else {
	yl1 <- range(c(x$location$median))
	yl2 <- range(c(x$scale$IQR))
	yl3 <- range(c(x$shape$skewness[,sel]))
	yl4 <- range(c(x$shape$shape[,sel]))
}

if(!all(is.finite(yl1))) yl1 <- NULL
if(!all(is.finite(yl2))) yl2 <- NULL
if(!all(is.finite(yl3))) yl3 <- NULL
if(!all(is.finite(yl4))) yl4 <- NULL

par(mfrow = c(2,2))
plot(z[r], x$location$median[r,], ylab = "Median", type = type, ylim = yl1, ...)
abline(h = 0, col = grey(.5))
if(interval){
lines(z[r], CI$Me[r,1], lty = 2, ...)
lines(z[r], CI$Me[r,2], lty = 2, ...)
}

plot(z[r], x$scale$IQR[r,], ylab = "IQR", type = type, ylim = yl2, ...)
abline(h = 0, col = grey(.5))
if(interval){
lines(z[r], CI$IQR[r,1], lty = 2, ...)
lines(z[r], CI$IQR[r,2], lty = 2, ...)
}

plot(z[r], x$shape$skewness[r,sel], ylab = "Skewness", type = type, ylim = yl3, ...)
abline(h = 0, col = grey(.5))
if(interval){
lines(z[r], CI$Ap[[sel]][r,1], lty = 2, ...)
lines(z[r], CI$Ap[[sel]][r,2], lty = 2, ...)
}

plot(z[r], x$shape$shape[r,sel], ylab = "Shape", type = type, ylim = yl4, ...)
abline(h = 0, col = grey(.5))
if(interval){
lines(z[r], CI$Tp[[sel]][r,1], lty = 2, ...)
lines(z[r], CI$Tp[[sel]][r,2], lty = 2, ...)
}

}

print.qlss <- function(x, ...){

if (!is.null(cl <- x$call)) {
        cat("call:\n")
        dput(cl)
        cat("\n")
}

n <- length(x$location$median)
txt <- if(n > 1) "Conditional Quantile-Based Location, Scale, and Shape" else "Unconditional Quantile-Based Location, Scale, and Shape"

if(n > 1){
cat(txt, "\n")
cat("-- Values are averaged over observations -- ", "\n")

cat("\n")

cat("** Location **", "\n")

cat("Median", "\n")
print(mean(x$location$median))

cat("** Scale **", "\n")

cat("Inter-quartile range (IQR)", "\n")
print(mean(x$scale$IQR))

cat("Inter-quantile range (IPR)", "\n")
print(colMeans(x$scale$IPR))

cat("**Shape**", "\n")

cat("Skewness index", "\n")
print(colMeans(x$shape$skewness))

cat("Shape index", "\n")
print(colMeans(x$shape$shape))

} else {

cat(txt, "\n")
cat("\n")

cat("** Location **", "\n")

cat("Median", "\n")
print(x$location$median)

cat("** Scale **", "\n")

cat("Inter-quartile range (IQR)", "\n")
print(x$scale$IQR)

cat("Inter-quantile range (IPR)", "\n")
print(x$scale$IPR)

cat("** Shape **", "\n")

cat("Skewness index", "\n")
print(x$shape$skewness)

cat("Shape index", "\n")
print(x$shape$shape)


}



}
