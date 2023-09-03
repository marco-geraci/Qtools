######################################################################
### Mid-distribution
######################################################################

# Sample mid-CDF and mid-QF

midecdf <- function(x, na.rm = FALSE){

if(!na.rm && any(is.na(x))) 
	return(NA)
if(na.rm && any(ii <- is.na(x))) 
	x <- x[!ii]
if(is.unsorted(x))
	x <- sort(x)

n <- length(x)
if (n < 1) 
	stop("'x' must have 1 or more non-missing values")
xo <- unique(x)

pmf <- as.numeric(table(x)/n)
val <- list()
val$call <- match.call()
val$x <- xo
val$y <- ecdf(x)(xo) - 0.5*pmf
val$fn <- approxfun(val$x, val$y, method = "linear", rule = 1)
val$data <- x
class(val) <- "midecdf"
return(val)

}

midquantile <- function(x, probs = 1:3/4, na.rm = FALSE){

if(!na.rm && any(is.na(x))) 
	return(NA)
if(na.rm && any(ii <- is.na(x))) 
	x <- x[!ii]
if(any(c(probs < 0, probs > 1)))
	stop("the probability must be between 0 and 1")
	
Fn <- midecdf(x)
Qn <- approxfun(Fn$y, Fn$x, method = "linear", rule = 2)
val <- list()
val$call <- match.call()
val$x <- probs
val$y <- Qn(probs)
val$fn <- Qn
val$data <- x
class(val) <- "midquantile"
return(val)
}

confint.midquantile <- function(object, parm = NULL, level = 0.95, ...){

x <- object$data
probs <- object$x

Fn <- midecdf(x, na.rm = TRUE)
k <- length(probs)
n <- length(x)
p <- table(x)/n
val <- dens <- rep(NA, k)
level <- level + (1-level)/2

for(i in 1:k){
	sel <- findInterval(probs[i], Fn$y)
	if(!sel %in% c(0, length(Fn$y))){
		lambda <- (Fn$y[sel+1] - probs[i])/(Fn$y[sel+1] - Fn$y[sel]);
		val[i] <- probs[i]*(1- probs[i]) - (1 - (lambda - 1)^2)*p[sel]/4 - (1 - lambda^2)*p[sel+1]/4
		dens[i] <- 0.5*(p[sel] + p[sel+1])/(Fn$x[sel+1] - Fn$x[sel])
	}
}
stderr <- sqrt(val/(n*dens^2))
LB <- object$y - qt(level, n - 1) * stderr
UB <- object$y + qt(level, n - 1) * stderr
val <- data.frame(midquantile = object$y, lower = LB, upper = UB)
rownames(val) <- paste0(probs*100, "%")
attr(val, "stderr") <- stderr
return(val)
}

# Conditional mid-CDF and mid-QF

cmidecdf <- function(formula, data, ecdf_est = "npc", npc_args = list(), theta = NULL, subset, weights, na.action, contrasts = NULL){

cl <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- quote(stats::model.frame)
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
intercept <- attr(terms.formula(formula, data = mf), "intercept") == 1
y <- model.response(mf, "numeric")
x <- model.matrix(mt, mf, contrasts)
w <- as.vector(model.weights(mf))
if (!is.null(w) && !is.numeric(w)) 
	stop("'weights' must be a numeric vector")

fit <- cmidecdf.fit(x = x, y = y, intercept = intercept, ecdf_est = ecdf_est, npc_args = npc_args, theta = theta)

return(fit)

}

cmidecdf.fit <- function(x, y, intercept, ecdf_est, npc_args = list(), theta = NULL){

	glm.ao <- function(x, y, theta){
	k <- length(theta)
	val <- rep(NA, k)
		for(i in 1:k){
			fit <- try(glm.fit(x, y, start = rep(0, ncol(x)), family = binomial(link = ao1(theta[i])))$deviance, silent = TRUE)
			if(!inherits(fit, "try-error")) val[i] <- fit
		}
	sel <- which.min(val)
	ans <- glm(y ~ x - 1, family = binomial(link = ao1(theta[sel])))
	ans$lambda <- theta[sel]
	return(ans)
	}

p <- ncol(x)
n <- length(y)
yo <- sort(unique(y))
K <- length(yo)
Z <- mapply(function(t, y) as.numeric(y <= t), yo, MoreArgs = list(y = y))
if(missing(intercept) & all(x[,1] == 1)) intercept <- TRUE

if(ecdf_est == "npc"){
	xnpc <- if(intercept) x[, 2:p, drop = FALSE] else x # remove intercept
	force_args <- list(xdat = xnpc, ydat = ordered(y), bandwidth.compute = is.null(npc_args$bws), oykertype = "wangvanryzin")
    npc_add <- setdiff(names(npc_args), names(force_args))
    npc_args <- c(force_args, npc_args[npc_add])
	bw <- try(fastDoCall(npcdistbw, npc_args), silent = TRUE)
	if(inherits(bw, "try-error")){
		ecdf_est <- "logit"
		warning("ecdf_est 'npc' failed, using 'logit' instead", "\n")
	} else {
		Fhat <- mapply(function(obj, x, y, n) npcdist(bws = obj, exdat = x, eydat = ordered(rep(y, n)))$condist, yo, MoreArgs = list(obj = bw, x = xnpc, n = n))
		Fse <- mapply(function(obj, x, y, n) npcdist(bws = obj, exdat = x, eydat = ordered(rep(y, n)))$conderr, yo, MoreArgs = list(obj = bw, x = xnpc, n = n))
	}
}

if(ecdf_est == "ao"){
	if(is.null(theta)) theta <- seq(0, 2, by = 0.05)
	fitbin <- apply(Z, 2, function(z, x, theta) suppressWarnings(glm.ao(x = x, y = z, theta = theta)), x = x, theta = theta)
	Fhat <- sapply(fitbin, predict, type = "response")
	Fse <- sapply(fitbin, function(x) predict(x, type = "response", se.fit = TRUE)$se.fit)
	bhat <- sapply(fitbin, coef) # p x K
	linkinv <- family(fitbin[[1]])$linkinv
}

if(ecdf_est %in% c("logit", "probit", "cloglog")){
	fitbin <- apply(Z, 2, function(z, x, link) suppressWarnings(glm(z ~ x - 1, family = binomial(link))), x = x, link = ecdf_est)
	Fhat <- sapply(fitbin, predict, type = "response")
	Fse <- sapply(fitbin, function(x) predict(x, type = "response", se.fit = TRUE)$se.fit)
	bhat <- sapply(fitbin, coef) # p x K
	linkinv <- family(fitbin[[1]])$linkinv
}

if(ecdf_est == "identity"){
	fitbin <- apply(Z, 2, function(z, x) suppressWarnings(lm(z ~ x - 1)), x = x)
	Fhat <- sapply(fitbin, predict, type = "response")
	Fse <- sapply(fitbin, function(x) predict(x, type = "response", se.fit = TRUE)$se.fit)
	bhat <- sapply(fitbin, coef) # p x K
	linkinv <- function(eta) eta
}

# rearrange if CDF is not monotone
for(j in 1:n){
	tmp <- Fhat[j,]
	if(any(diff(tmp) < 0)){
		sf <- rearrange(stepfun(yo, c(-Inf, tmp)))
		Fhat[j,] <- sf(yo)
	}
}


M <- apply(Fhat, 1, diff)
if(ncol(Fhat) > 2) M <- t(M)
G <- Fhat[,-1] - 0.5*M
G <- cbind(Fhat[,1]/2, G)
r <- c(max(G[,1]), min(G[,ncol(G)]))

attr(G, "range") <- r
ecdf_fit <- if(ecdf_est == "npc") bw else list(coef = bhat, linkinv = linkinv)

ans <- list(G = G, Fhat = Fhat, Fse = Fse, yo = yo, ecdf_fit = ecdf_fit, ecdf_est = ecdf_est)
class(ans) <- "cmidecdf"

return(ans)

}

midrqControl <- function(method = "Nelder-Mead", ecdf_est = "npc", npc_args = list()){

	list(method = method, ecdf_est = ecdf_est, npc_args = npc_args)

}

midrq <- function(formula, data, tau = 0.5, lambda = NULL, subset, weights, na.action, contrasts = NULL, offset, type = 3, midFit = NULL, control = list()){

cl <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action", 
	"offset"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- quote(stats::model.frame)
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
intercept <- attr(terms.formula(formula, data = mf), "intercept") == 1
y <- model.response(mf, "numeric")
x <- model.matrix(mt, mf, contrasts)
w <- as.vector(model.weights(mf))
if (!is.null(w) && !is.numeric(w)) 
	stop("'weights' must be a numeric vector")
offset <- as.vector(model.offset(mf))
if (!is.null(offset)){
	if (length(offset) != NROW(y)) 
		stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
			length(offset), NROW(y)), domain = NA)
}

if (is.null(names(control))) 
	control <- midrqControl()
else {
	control_default <- midrqControl()
	control_names <- intersect(names(control), names(control_default))
	control_default[control_names] <- control[control_names]
	control <- control_default
}

nq <- length(tau)
n <- nrow(x)
p <- ncol(x)
if (is.null(offset)) offset <- rep(0, n)
#x <- normalize(x) * sqrt(n)

# is the response binary?
binary <- setequal(sort(unique(y)), c(0,1))
if(binary){
	cat("Binary response detected", "\n")
}

# response on the scale of linear predictor
if(!is.null(lambda)){
	if(binary){
		hy <- ao(y, lambda)
	} else {
		if(any(y <= 0)) stop("The response must be strictly positive")
		hy <- bc(y, lambda)
	}
} else {hy <- y}

# fit model
if(p == 1 & intercept){
	midFit <- midquantile(y, probs = tau)
	fit <- list()
	for (j in 1:nq) {
		fit[[j]] <- list(par = midFit$y[j])
	}	
} else {

	# Estimate mid-CDF if not provided
	if(is.null(midFit)){
		midFit <- cmidecdf.fit(x = x, y = y, intercept = intercept, ecdf_est = control$ecdf_est, npc_args = control$npc_args, theta = seq(0, 2, by = 0.05))
	} else {control$ecdf_est <- midFit$ecdf_est}

	fit <- list()
	for (j in 1:nq) {
		fit[[j]] <- midrq.fit(x = x, y = y, offset = offset, lambda = lambda, binary = binary, midFit = midFit, type = type, tau = tau[j], method = control$method)
	}

}

names(fit) <- tau
betahat <- matrix(sapply(fit, function(x) x$par), nrow = p, ncol = nq)
rownames(betahat) <- colnames(x)
colnames(betahat) <- tau

yhat <- x %*% betahat + offset
if(!is.null(lambda)){
	if(binary){
		Fitted <- apply(yhat, 2, invao, lambda = lambda)
	} else {
		Fitted <- apply(yhat, 2, invbc, lambda = lambda)
	}
} else {
	Fitted <- yhat
}

fit$call <- cl
#fit$mf <- mf
fit$formula <- formula
mf <- match.call(expand.dots = FALSE)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("get_all_vars")
fit$data <- eval(mf, parent.frame())
fit$x <- x
fit$y <- y
fit$hy <- hy
fit$weights <- w
fit$offset <- offset
fit$tau <- tau
fit$lambda <- lambda
fit$binary <- binary
fit$intercept <- intercept
fit$type <- type
fit$control <- control
fit$coefficients <- betahat
fit$fitted.values <- Fitted
fit$residuals <- y - Fitted
fit$midFit <- midFit
fit$levels <- .getXlevels(mt, mf)
fit$terms <- mt
fit$term.labels <- colnames(x)

class(fit) <- "midrq"
return(fit)
}

midrq.fit <- function(x, y, offset, lambda, binary, midFit, type, tau, method){

n <- nrow(x)
p <- ncol(x)

# transform response to get starting values on the scale of linear predictor
if(!is.null(lambda)){
	if(binary){
		z <- ao(y, lambda) - offset
	} else {
		if(any(y <= 0)) stop("The response must be strictly positive")
		z <- bc(y, lambda) - offset
	}
} else {z <- y - offset}

if(type %in% c(1,2)){
	# starting values for beta
	b0 <- lm.fit(x, z)$coefficients

	if(p == 1){
		if(is.null(lambda)){
			fit <- optim(par = b0, fn = C_midrqLoss, G = midFit$G, x = x, yo = midFit$yo, offset = offset, type = type, tau = tau, n = n, p = p, k = length(midFit$yo), method = "Brent", lower = -Inf, upper = Inf)
			warning("'type = 3' is more reliable when no-intercept model is fitted")
		} else {
			if(binary){
				fit <- optim(par = b0, fn = C_midrqLoss_ao, G = midFit$G, x = x, yo = midFit$yo, offset = offset, type = type, tau = tau, lambda = lambda, n = n, p = p, k = length(midFit$yo), method = "Brent", lower = -Inf, upper = Inf)
			} else {
				fit <- optim(par = b0, fn = C_midrqLoss_bc, G = midFit$G, x = x, yo = midFit$yo, offset = offset, type = type, tau = tau, lambda = lambda, n = n, p = p, k = length(midFit$yo), method = "Brent", lower = -Inf, upper = Inf)
			}
		}
	} else {
		if(is.null(lambda)){
			fit <- optim(par = b0, fn = C_midrqLoss, G = midFit$G, x = x, yo = midFit$yo, offset = offset, type = type, tau = tau, n = n, p = p, k = length(midFit$yo), method = method)
		} else {
			if(binary){
				fit <- optim(par = b0, fn = C_midrqLoss_ao, G = midFit$G, x = x, yo = midFit$yo, offset = offset, type = type, tau = tau, lambda = lambda, n = n, p = p, k = length(midFit$yo), method = method)
			} else {
				fit <- optim(par = b0, fn = C_midrqLoss_bc, G = midFit$G, x = x, yo = midFit$yo, offset = offset, type = type, tau = tau, lambda = lambda, n = n, p = p, k = length(midFit$yo), method = method)
			}
		}
		fit$pseudoy <- NULL
	}
	
	fit$InitialPar <- b0
	
	} else if(type == 3){
	k <- length(midFit$yo)	
	up <- apply(tau - midFit$G, 1, function(x) which(x < 0)[1])
	FLAG <- sum(up == 1, na.rm = TRUE) + sum(is.na(up))
	low <- up - 1
	low[low == 0] <- 1
	low[is.na(low)] <- k
	up[is.na(up)] <- k
	Z <- cbind(midFit$yo[low], midFit$yo[up])
	PI <- t(apply(midFit$G, 1, function(x, p){
		up <- which(p - x < 0)[1]
		low <- up - 1
		low[low == 0] <- 1
		low[is.na(low)] <- length(x)
		up[is.na(up)] <- length(x)
		x[c(low, up)]
	}, p = tau))
	gamma <- (tau - PI[,1])/(PI[,2] - PI[,1])
	gamma[!is.finite(gamma)] <- 0
	B <- gamma*(Z[,2] - Z[,1]) + Z[,1]
	if(!is.null(lambda)){
		if(binary){
			B <- ao(B, lambda) - offset
		} else {
			B <- bc(B, lambda) - offset
		}
	} else {B <- B - offset}
	
	# is tau outside range?
	r <- attr(midFit$G, "range")
	if(any(tau < r[1]) | any(tau > r[2])){
		warning("tau = ", tau, " is outside mid-probabilities range ", "[", round(r[1], 3), ", " , round(r[2], 3), "] for ", FLAG, " out of ", n, " observations. See details for ?midrq")
	}

	fit <- list(par = qr.solve(x, B), pseudoy = B)

	}

return(fit)
}

fitted.midrq <- function(object, ...){

return(object$fitted.values)

}

residuals.midrq <- function(object, ...){

return(object$residuals)

}

predict.midrq <- function(object, newdata, offset, na.action = na.pass, type = "response", ...){

lambda <- object$lambda
if (!missing(newdata)) {
	mt <- terms(object)
	Terms <- delete.response(mt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	object$x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
	if (missing(offset)) 
		object$offset <- rep(0, nrow(object$x))
	object$data <- newdata
}
linpred <- drop(object$x %*% object$coefficients) + object$offset
if (type == "link") {
	return(linpred)
}
if (type == "response") {
	if(!is.null(lambda)){
		Fitted <- apply(linpred, 2, if(object$binary) invao else invbc, lambda = lambda)
	} else {
		Fitted <- linpred
	}
	return(Fitted)
}

}

coef.midrq <- coefficients.midrq <- function(object, ...){

return(object$coefficients)

}

vcov.midrq <- function(object, numerical = FALSE, robust = FALSE, ...){

	phi <- function(xnz, Fvec, nonZero, Z, w, n, k, tau, yo, offset, binary, lambda = NULL){
		Fvec[nonZero] <- xnz 
		Fhat <- matrix(Fvec, n, k)
		M <- apply(Fhat, 1, diff)
		if(ncol(Fhat) > 2) M <- t(M)
		G <- Fhat[,-1] - 0.5*M
		G <- cbind(Fhat[,1]/2, G)
			
		PI <- t(apply(G, 1, function(x, p){
			sel <- which(p - x < 0)[1]
			x[c(sel-1, sel)]
			}, p = tau))
		B <- (tau - PI[,1])/(PI[,2] - PI[,1])*(Z[,2] - Z[,1]) + Z[,1]

		if(!is.null(lambda)){
		if(binary){
			B <- ao(B, lambda) - offset
		} else {
			B <- bc(B, lambda) - offset
		}
		} else {B <- B - offset}

		ans <- qr.solve(w, B)
		return(ans)
	}
	
	huber <- function(x, a = 1.645){
		ifelse(abs(x) <= a, 0.5*x^2, a*(abs(x) - 0.5*a))
	}
	
	f <- function(b, args) {
		args$b <- b
		do.call(C_midrqLoss, args = args)
	}

tau <- object$tau
nq <- length(tau)
x <- object$x
n <- length(object$y)
p <- ncol(x)
G <- object$midFit$G
Gvec <- as.vector(G)
yo <- object$midFit$yo
K <- length(yo)

if(object$intercept & p == 1){
	V <- as.list(attr(confint(object$midFit), "stderr")^2)
	names(V) <- tau
	return(V)
}

if(numerical){
	FIT_ARGS <- list(G = G, x = x, yo = yo, offset = object$offset, type = object$type, n = n, p = p, k = K)

	V <- list()
	for(j in 1:nq){
		FIT_ARGS$tau <- tau[j]
		H <- hessian(func = f, x = object$coefficients[,j], method = "Richardson", args = FIT_ARGS)
		ans <- MASS::ginv(H)
		V[[j]] <- ans/n
	}
} else {
	xb <- as.matrix(predict(object, type = "link")) # xb includes the offset
	xx <- solve(crossprod(x))
	rate <- if(object$midFit$ecdf_est == "npc") prod(object$midFit$ecdf_fit$xbw)*n else 1
	Fvec <- as.vector(object$midFit$Fhat)
	r <- attr(G, "range")

	# variance of Fhat
	J1 <- object$midFit$Fse^2

	V <- list()
	for(j in 1:nq){
		res <- if(robust) huber(object$hy - xb[,j]) else (object$hy - xb[,j])^2
		V1 <- xx %*% t(x) %*% Diagonal(x = res) %*% x %*% xx

		up <- apply(tau[j] - G, 1, function(x) which(x < 0)[1])
		low <- up - 1
		if(any(low == 0) | any(is.na(low))) stop("Something went wrong. Perhaps tau is outside allowed range ", "[", round(r[1], 3), ", " , round(r[2], 3), "]")
		nonZero <- c((low - 1)*n + 1:n, (up - 1)*n + 1:n)
		
		J2 <- jacobian(func = phi, x = Fvec[nonZero], Fvec = Fvec, nonZero = nonZero, Z = cbind(yo[low], yo[up]), method = "simple", w = x, n = n, k = K, tau = tau[j], yo = yo, offset = object$offset, binary = object$binary, lambda = object$lambda, method.args = list(eps = 1e-6))
		V2 <- J2 %*% Diagonal(x = J1[nonZero]) %*% t(J2)

		V[[j]] <- as.matrix(V1 + rate*V2)
	}
}

names(V) <- tau
return(V)
}

summary.midrq <- function(object, alpha = 0.05, numerical = FALSE, robust = FALSE, ...){

tau <- object$tau
nq <- length(tau)
p <- ncol(object$x)
bhat <- object$coefficients

if(object$intercept & p == 1){
	tmp <- confint(object$midFit, level = 1 - alpha)
	SE <- matrix(attr(tmp, "stderr"), nrow = 1)
	lower <- matrix(tmp$lower, nrow = 1)
	upper <- matrix(tmp$upper, nrow = 1)
} else {
	SE <- sapply(vcov(object, numerical = numerical, robust = robust), function(x) sqrt(diag(x)))
	lower <- bhat - SE*qnorm(1 - alpha/2, 0, 1)
	upper <- bhat + SE*qnorm(1 - alpha/2, 0, 1)
}

if(nq == 1){
	tTable <- data.frame(bhat, SE, lower, upper)
	names(tTable) <- c("Estimate", "Std.Err", "Lower", "Upper")
} else {
	tTable <- list()
	for(i in 1:nq){
		tTable[[i]] <- data.frame(bhat[,i], SE[,i], lower[,i], upper[,i])
		dimnames(tTable[[i]]) <- list(object$term.labels, c("Estimate", "Std.Err", "Lower", "Upper"))
	}
	names(tTable) <- tau
}

object$tTable <- tTable
class(object) <- "summary.midrq"

return(object)

}

midq2q <- function(object, newdata, observed = FALSE, ...){

tau <- object$tau
nt <- length(tau)
ecdf_fit <- object$midFit$ecdf_fit
yo <- object$midFit$yo
x <- model.matrix(object$formula[-2], newdata)
n <- nrow(x)
p <- ncol(x)
xnpc <- if(object$intercept) x[, 2:p, drop = FALSE] else x
if(object$midFit$ecdf_est == "npc"){
	Fhat <- mapply(function(obj, x, y, n) npcdist(bws = obj, exdat = x, eydat = ordered(rep(y, n)))$condist, yo, MoreArgs = list(obj = ecdf_fit, x = xnpc, n = n))
} else {
	Fhat <- apply(x %*% ecdf_fit$coef, 2, ecdf_fit$linkinv)
}

# rearrange if CDF is not monotone
for(j in 1:n){
	tmp <- Fhat[j,]
	if(any(diff(tmp) < 0)){
		sf <- rearrange(stepfun(yo, c(-Inf, tmp)))
		Fhat[j,] <- sf(yo)
	}
}

M <- apply(Fhat, 1, diff)
if (ncol(Fhat) > 2) 
	M <- t(M)
Ghat <- Fhat[, -1] - 0.5 * M
Ghat <- cbind(Fhat[, 1]/2, Ghat)
Hhat <- predict(object, newdata = newdata)
csi <- tmp <- matrix(NA, n, nt)
for(j in 1:n){
	sel <- findInterval(Hhat[j,], yo, all.inside = TRUE)
	low <- yo[sel] 
	up <- yo[sel + 1]
	if(observed){
		csi[j,] <- ifelse(tau > Fhat[j,sel], up, low)
		tmp[j,] <- Fhat[j,sel]
	} else {
		gamma <- (Hhat[j,] - low)/(up - low)
		gamma <- pmax(gamma, 0)
		gamma <- pmin(gamma, 1)
		pstar <- as.numeric((1-gamma)*Ghat[j,sel] + gamma*Ghat[j,sel + 1])
		sel <- findInterval(pstar, Ghat[j,])
		csi[j,] <- ifelse(pstar > Fhat[j,sel], ceiling(Hhat[j,]), floor(Hhat[j,]))
		tmp[j,] <- Fhat[j,sel]
	}
}
colnames(csi) <- tau
rownames(csi) <- rownames(x)
attr(csi, "Fhat") <- tmp
class(csi) <- "midq2q"

return(csi)

}

fastDoCall <- function(what, args, quote = FALSE, envir = parent.frame()) {

# Source: Gmisc 1.11.0 (Max Gordon)

  if (quote) {
    args <- lapply(args, enquote)
  }

  if (is.null(names(args)) ||
    is.data.frame(args)) {
    argn <- args
    args <- list()
  } else {
    # Add all the named arguments
    argn <- lapply(names(args)[names(args) != ""], as.name)
    names(argn) <- names(args)[names(args) != ""]
    # Add the unnamed arguments
    argn <- c(argn, args[names(args) == ""])
    args <- args[names(args) != ""]
  }

  if ("character" %in% class(what)) {
    if (is.character(what)) {
      fn <- strsplit(what, "[:]{2,3}")[[1]]
      what <- if (length(fn) == 1) {
        get(fn[[1]], envir = envir, mode = "function")
      } else {
        get(fn[[2]], envir = asNamespace(fn[[1]]), mode = "function")
      }
    }
    call <- as.call(c(list(what), argn))
  } else if ("function" %in% class(what)) {
    f_name <- deparse(substitute(what))
    call <- as.call(c(list(as.name(f_name)), argn))
    args[[f_name]] <- what
  } else if ("name" %in% class(what)) {
    call <- as.call(c(list(what, argn)))
  }

  eval(call,
    envir = args,
    enclos = envir
  )
}

# Plot and print functions

plot.midq2q <- function(x, ..., xlab = "p", ylab = "Quantile", main = "Ordinary Quantile Function", sub = TRUE, verticals = TRUE, col.steps = "gray70", cex.points = 1, jumps = FALSE){

n <- nrow(x)
k <- ncol(x)
Fhat <- attr(x, "Fhat")

if(n > 1) par(mfrow = c(ceiling(n/3), min(c(2,n))))

for(j in 1:n){
	sf <- stepfun(Fhat[j,], c(x[j,], x[j,k]))
	subtext <- paste("id = ", j)
	plot.stepfun(sf, xlab = xlab, ylab = ylab, main = main, verticals = verticals, col = col.steps, cex.points = cex.points, col.points = col.steps, do.points = jumps, xlim = c(0,1), ...)
	if(sub) mtext(text = subtext, side = 3, line = 0.5, cex = 0.8)
}

}

plot.midecdf <- function(x, ..., ylab = "p", main = "Ordinary and Mid-ECDF", verticals = FALSE, col.01line = "gray70", col.steps = "gray70", col.midline ="black", cex.points = 1, lty.midline = 2, lwd = 1, jumps = FALSE){

z <- sort(x$data)
n <- length(z)
vals <- unique(z)
nv <- length(vals)
Fz <- c(0, cumsum(table(z)/n))

pch <- if(jumps) 19 else ""
plot.stepfun(ecdf(x$data), ylab = ylab, main = main, verticals = verticals, col = col.steps, pch = pch, cex.points = cex.points, col.points = col.steps, ...)

if(jumps){
	points(vals, Fz[-nv], cex = cex.points, col = col.steps)
}

lines(x$x, x$fn(x$x), lty = lty.midline, col = col.midline, lwd = lwd)
abline(h = c(0, 1), col = col.01line, lty = 2)

}

plot.midquantile <- function(x, ..., xlab = "p", ylab = "Quantile", main = "Ordinary and Mid-Quantiles", col.steps = "gray70", col.midline ="black", cex.points = 1, lty.midline = 2, lwd = 1, jumps = FALSE){

z <- sort(x$data)
n <- length(z)
vals <- unique(z)
nv <- length(vals)
Fz <- c(0, cumsum(table(z)/n))
plot(c(0,1), range(z), axes = TRUE, type = "n", xlab = xlab, ylab = ylab, main = main, ...)

segments(x0 = Fz[1:nv], y0 = vals, x1 = Fz[-1], y1 = vals, col = col.steps)
if(jumps){
	points(Fz[1:nv], vals, cex = cex.points, col = col.steps)
	points(1, vals[nv], cex = cex.points, col = col.steps)
	points(Fz[-c(1,nv+1)], vals[-nv], pch = 19, cex = cex.points, col = col.steps)
}

lines(x$x, x$fn(x$x), lty = lty.midline, col = col.midline, lwd = lwd)
}

print.midecdf <- function(x, ...){

cat("Empirical mid-ECDF", "\n")
if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
}

}

print.midquantile <- function(x, ...){

cat("Empirical mid-quantile function", "\n")
if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
}


}

print.cmidecdf <- function(x, ...){

ecdf_est <- switch(x$ecdf_est,
			npc = "Nonparametric (kernel)",
			ao = "Aranda-Ordaz",
			logit = "logit",
			probit = "probit",
			cloglog = "cloglog",
			identity = "linear")

cat("Empirical conditional mid-ECDF", "\n")
cat("Estimator:", ecdf_est)
cat("\n")
if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
}

}

print.midrq <- function(x, ...){

if (!is.null(cl <- x$call)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}

coef <- x$coefficients
cat("\nCoefficients linear predictor:\n")
print(coef, ...)
nobs <- length(x$y)
p <- ncol(x$x)
rdf <- nobs - p
cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")

}

print.summary.midrq <- function(x, ...){

if (!is.null(cl <- x$call)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}

cat("\nCoefficients linear predictor:\n")
print(x$tTable, ...)
nobs <- length(x$y)
p <- ncol(x$x)
rdf <- nobs - p
cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")

}

