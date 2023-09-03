##################################################
### Restricted quantiles
##################################################

rrq <- function(formula, tau, data, subset, weights, na.action, method = "fn", model = TRUE, contrasts = NULL, ...){

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval.parent(mf)
if (method == "model.frame")
	return(mf)
mt <- attr(mf, "terms")
weights <- as.vector(model.weights(mf))
y <- model.response(mf)
x <- model.matrix(mt, mf, contrasts)

eps <- .Machine$double.eps^(2/3)
nq <- length(tau)

if (nq > 1) {
        if (any(tau < 0) || any(tau > 1)) 
            stop("invalid tau:  taus should be >= 0 and <= 1")
        if (any(tau == 0)) 
            tau[tau == 0] <- eps
        if (any(tau == 1)) 
            tau[tau == 1] <- 1 - eps
}

fit.lad <- {
if (length(weights)) 
	rq.wfit(x, y, tau = 0.5, weights, method, ...)
	else rq.fit(x, y, tau = 0.5, method, ...)
}
r.lad <- fit.lad$residuals
r.abs <- abs(fit.lad$residuals)
beta <- fit.lad$coefficients
p <- length(beta)

#fit.lad <- rq(formula, tau = 0.5, data = data, method = method)
#data$r.lad <- fit.lad$residuals
#data$r.abs <- abs(fit.lad$residuals)
#beta <- fit.lad$coefficients

fit.lad <- {
if (length(weights)) 
	rq.wfit(x, r.abs, tau = 0.5, weights, method, ...)
	else rq.fit(x, r.abs, tau = 0.5, method, ...)
}
s.lad <- fit.lad$fitted.values
gamma <- fit.lad$coefficients

#fit.lad <- rq(stats::update.formula(formula, r.abs ~ .), tau = 0.5, data = data, method = method)
#data$s.lad <- fit.lad$fitted
#gamma <- fit.lad$coefficients

zeta <- rep(0, nq)
for (i in 1:nq) {
	zeta[i] <- {
	if (length(weights)) 
		rq.wfit(matrix(s.lad), matrix(r.lad), tau = tau[i], weights, method, ...)$coefficients
		else rq.fit(matrix(s.lad), matrix(r.lad), tau = tau[i], method, ...)$coefficients
	}
}

#zeta <- rq(r.lad ~ s.lad - 1, tau = tau, data = data, method = method)$coefficients

if (nq > 1){
	coef <- apply(outer(matrix(gamma, nrow = 1), zeta, "*"), 3, function(x, b) x + b, b = beta, simplify = TRUE)
	coef <- matrix(coef, nrow = p)
	taulabs <- paste0("tau = ", format(round(tau, 3)))
	dimnames(coef) <- list(dimnames(x)[[2]], taulabs)
} else {
	coef <- as.numeric(beta + zeta * gamma)
}

fit <- list(coefficients = coef, zeta = zeta, beta = beta, gamma = gamma, tau = tau)
fit$na.action <- attr(mf, "na.action")
fit$formula <- formula
fit$terms <- mt
fit$xlevels <- .getXlevels(mt, mf)
fit$call <- call
fit$weights <- weights
fit$method <- method
fit$x <- x
fit$y <- y
fit$fitted.values <- drop(x %*% coef)
fit$residuals <- drop(y - x %*% coef)
attr(fit, "na.message") <- attr(m, "na.message")
if (model)
	fit$model <- mf
class(fit) <- "rrq"
return(fit)
}

rrq.fit <- function(x, y, tau, method = "fn", ...){

if(length(tau) > 1) stop("only one quantile")

fit.lad <- rq.fit(x, y, tau = 0.5, method = method, ...)
r.lad <- fit.lad$residuals
r.abs <- abs(fit.lad$residuals)
beta <- fit.lad$coefficients

fit.lad <- rq.fit(x, r.abs, tau = 0.5, method = method, ...)
s.lad <- fit.lad$fitted.values
gamma <- fit.lad$coefficients

zeta <- rq.fit(s.lad, r.lad, tau = tau, method = method, ...)$coefficients

val <- beta + zeta * gamma

return(list(coefficients = val, zeta = zeta, beta = beta, gamma = gamma, tau = tau))
}

rrq.wfit <- function(x, y, tau, weights, method = "fn", ...){

if(length(tau) > 1) stop("only one quantile")

fit.lad <- rq.wfit(x, y, tau = 0.5, weights, method = method, ...)
r.lad <- fit.lad$residuals
r.abs <- abs(fit.lad$residuals)
beta <- fit.lad$coefficients

fit.lad <- rq.wfit(x, r.abs, tau = 0.5, weights, method = method, ...)
s.lad <- fit.lad$fitted.values
gamma <- fit.lad$coefficients

zeta <- rq.wfit(s.lad, r.lad, tau = tau, weights, method = method, ...)$coefficients

val <- beta + zeta * gamma

return(list(coefficients = val, zeta = zeta, beta = beta, gamma = gamma, tau = tau, weights = weights))
}

boot.rrq <- function(data, inds, object){

tau <- object$tau
nq <- length(tau)
all.obs <- rownames(object$x)
n <- length(all.obs)
nn <- dimnames(object$coefficients)[[1]]

fit <- stats::update(object, data = data[inds,])
val <- fit$coefficients

val <- as.vector(val)
names(val) <- rep(nn, nq)
return(val)

}

summary.rrq <- function(object, alpha = 0.05, se = "boot", R = 50, sim = "ordinary", stype = "i", ...){

call <- match.call(expand.dots = TRUE)

tau <- object$tau
nq <- length(tau)
ntot <- ncol(object$x)

if(se == "boot"){
	Args <- list()
	Args$data <- object$model
	Args$statistic <- boot.rrq
	Args$object <- object
	Args$R <- R
	Args$sim <- sim
	Args$stype <- stype
	nn <- c("strata","L","m","weights","ran.gen","mle","simple","parallel","ncpus","cl")
	nn <- nn[pmatch(names(call), nn, duplicates.ok = FALSE)]
	nn <- nn[!is.na(nn)]
	if(length(nn) > 0) {tmp <- as.list(call[[nn]]); names(tmp) <- nn; Args <- c(Args, tmp)}
	B <- do.call(boot, args = Args)
	ci <- mapply(boot.ci, index = 1:(ntot*nq), MoreArgs = list(boot.out = B, conf = 1 - alpha, type = "perc"))[4,]
	ci <- t(sapply(ci, function(x) x[4:5]))

	S <- cov(B$t, use = "complete.obs")
	val <- cbind(B$t0, apply(B$t, 2L, mean, na.rm=TRUE) - B$t0, sqrt(diag(S)), ci)
	nn <- c("Value", "Bias", "Std. Error", "Lower bound", "Upper bound")
	colnames(val) <- nn
	
	maxn <- seq(0, ntot*nq, by = ntot)[-1]
	minn <- seq(1, ntot*nq, by = ntot)
	ans <- list()
	for(j in 1:nq){
		ans[[j]] <- val[minn[j]:maxn[j], ]
	}
	names(ans) <- tau
	object$B <- B
} else {ans <- NULL}


object$coefficients <- ans
object$call <- call
class(object) <- c("summary.rrq", class(object))
return(object)

}

predict.rrq <- function(object, newdata, na.action = na.pass, ...){


tau <- object$tau
nq <- length(tau)
betahat <- object$coefficients

if(missing(newdata)) {x <- object$x} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
}

return(x %*% betahat)

}

print.rrq <- function(x, ...){

class(x) <- "rqs"
print.rqs(x, ...)

}

print.summary.rrq <- function(x, ...){

if (!is.null(cl <- x$call)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}

tau <- x$tau
nq <- length(tau)
mpar <- ncol(x$x)

cat("\nSummary for restricted regression quantiles\n")

for(i in 1:nq){
cat("\ntau = ", tau[i], "\n")

cat("\nCoefficients linear model:\n")
print(x$coefficients[[i]][1:mpar,], ...)
}

nobs <- length(x$y)
p <- ncol(x$x)
rdf <- nobs - p
cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
if (!is.null(attr(x, "na.message"))) 
	cat(attr(x, "na.message"), "\n")
invisible(x)

}
