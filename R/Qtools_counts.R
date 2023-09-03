##################################################
### QR for counts (Machado and Santos Silva)
##################################################

rq.counts <- function(formula, data = sys.frame(sys.parent()), tau = 0.5, subset, weights, na.action, contrasts = NULL, offset = NULL, method = "fn", M = 50, zeta = 1e-5, B = 0.999, cn = NULL, alpha = 0.05) 
{

# log-linear model
tsf <- "bc"
symm <- TRUE
dbounded <- FALSE
lambda <- 0 

nq <- length(tau)
if (nq > 1) 
	stop("One quantile at a time")

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
mf <- mf[c(1L, m)]
if (method == "model.frame")
	return(mf)
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")

x <- model.matrix(mt, mf, contrasts)
y <- model.response(mf, "numeric")

p <- ncol(x)
n <- nrow(x)
term.labels <- colnames(x)

offset <- model.offset(mf)
if(is.null(offset)) offset <- rep(0, n)
w <- as.vector(model.weights(mf))
if(is.null(w)) w <- rep(1, n)

Fn <- function(x, cn){
	xf <- floor(x)
	df <- x - xf
	if(df < cn & x >= 1){
		val <- xf - 0.5 + df/(2*cn)
	}
	if(any(cn <= df & df < (1 - cn), x < 1)){
		val <- xf
	}

	if(df >= (1 - cn)){
		val <- xf + 0.5 + (df - 1)/(2*cn)
	}

	return(val)
}

Fvec <- Vectorize(Fn)

# Add noise
Z <- replicate(M, addnoise(y, centered = FALSE, B = B))

# Transform Z
TZ <- apply(Z, 2, function(x, off, tsf, symm, lambda, tau, zeta){
	z <- ifelse((x - tau) > zeta, x - tau, zeta);
	switch(tsf,
		mcjI = mcjI(z, lambda, symm, dbounded = dbounded, omega = 0.001),
		bc = bc(z, lambda)) - off
	}, off = offset, tsf = tsf, symm = symm, lambda = lambda, tau = tau, zeta = zeta)

# Fit linear QR on TZ
fit <- apply(TZ, 2, function(y, x, weights, tau, method) 
	rq.wfit(x = x, y = y, tau = tau, weights = weights, method = method), x = x, tau = tau, weights = w, method = method)
	
# predicted values
yhat <- sapply(fit, function(obj, x) x %*% obj$coefficients, x = x)
yhat <- as.matrix(yhat)

# sweep offset back in
linpred <- sweep(yhat, 1, offset, "+")

# back-transform + offset tau
zhat <- matrix(NA, n, M)
for(i in 1:M){
zhat[,i] <- tau + switch(tsf,
	mcjI = invmcjI(linpred[,i], lambda, symm, dbounded = dbounded),
	bc = invbc(linpred[,i], lambda))
}
	
# covariance matrix
if(is.null(cn)) cn <- 0.5 * log(log(n))/sqrt(n)
F <- apply(zhat, 2, Fvec, cn = cn)
Fp <- apply(zhat + 1, 2, Fvec, cn = cn)

multiplier <- (tau - (TZ <= yhat))^2
a <- array(NA, dim = c(p, p, M))
for (i in 1:M) a[, , i] <- t(x * multiplier[, i]) %*% x/n

multiplier <- tau^2 + (1 - 2 * tau) * (y <= (zhat - 1)) + 
	((zhat - y) * (zhat - 1 < y & y <= zhat)) * (zhat - y - 
		2 * tau)
b <- array(NA, dim = c(p, p, M))
for (i in 1:M) b[, , i] <- t(x * multiplier[, i]) %*% x/n

multiplier <- (zhat - tau) * (F <= Z & Z < Fp)
d <- array(NA, dim = c(p, p, M))
sel <- rep(TRUE, M)
for (i in 1:M) {
	tmpInv <- try(solve(t(x * multiplier[, i]) %*% x/n), 
		silent = TRUE)
	if (!inherits(tmpInv, "try-error"))
		{d[, , i] <- tmpInv}
	else {sel[i] <- FALSE}
}
    
dad <- 0
dbd <- 0
for (i in (1:M)[sel]) {
	dad <- dad + d[, , i] %*% a[, , i] %*% d[, , i]
	dbd <- dbd + d[, , i] %*% b[, , i] %*% d[, , i]
}
    
m.n <- sum(sel)
if (m.n != 0) {
	V <- dad/(m.n^2) + (1 - 1/m.n) * dbd * 1/m.n
	V <- V/n ## CHECK V AND WEIGHTS
	stds <- sqrt(diag(V))
	} else {
	V <- NA
	stds <- NA
	warning("Standard error not available")
	}

betahat <- sapply(fit, function(x) x$coefficients)
betahat <- if (p == 1) mean(betahat) else rowMeans(betahat)

linpred <- if (p == 1) {
	mean(linpred[1, ])
} else {
	rowMeans(linpred)
}

Fitted <- tau + switch(tsf,
	mcjI = invmcjI(linpred, lambda, symm, dbounded = dbounded),
	bc = invbc(linpred, lambda))

lower <- betahat + qt(alpha/2, n - p) * stds
upper <- betahat + qt(1 - alpha/2, n - p) * stds
tP <- 2 * pt(-abs(betahat/stds), n - p)

ans <- cbind(betahat, stds, lower, upper, tP)
colnames(ans) <- c("Value", "Std. Error", "lower bound", "upper bound", "Pr(>|t|)")
rownames(ans) <- names(betahat) <- term.labels

fit <- list()
fit$call <- call
fit$method <- method
fit$mf <- mf
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("get_all_vars")
fit$data <- eval(mf, parent.frame())
fit$x <- x
fit$y <- y
fit$weights <- w
fit$offset <- offset
fit$tau <- tau
fit$lambda <- lambda
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
fit$coefficients <- betahat
fit$M <- M
fit$Mn <- m.n
fit$fitted.values <- Fitted
fit$tTable <- ans
fit$Cov <- V
fit$levels <- .getXlevels(mt, mf)
fit$terms <- mt
fit$term.labels <- term.labels
fit$rdf <- n - p

class(fit) <- "rq.counts"
	
return(fit)
}

coef.rq.counts <- coefficients.rq.counts <- function(object, ...){

tau <- object$tau
nq <- length(tau)
ans <- object$coefficients

if(nq == 1){
  names(ans) <- object$term.labels
}

return(ans)

}

fitted.rq.counts <- function(object, ...){

return(object$fitted.values)

}

predict.rq.counts <- function(object, newdata, offset, na.action = na.pass, type = "response", namevec = NULL, ...) 
{

tsf <- object$tsf
symm <- attributes(tsf)$symm
lambda <- object$lambda

if(!missing(newdata)){
	mt <- terms(object)
	Terms <- delete.response(mt)
	m <- object$mf <- model.frame(Terms, newdata, na.action = na.action, xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
	object$x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
	if(missing(offset)) object$offset <- rep(0, nrow(object$x))
	object$data <- newdata
}

linpred <- drop(object$x %*% object$coefficients) + object$offset

if(type == "link"){
	return(linpred)
}

if(type == "response"){
	Fitted <- object$tau + switch(tsf, mcjI = invmcjI(linpred, lambda, symm, dbounded =FALSE), bc = invbc(linpred, lambda))
	return(Fitted)
}

if(type == "maref"){
	if(is.null(namevec)) stop("When type = 'maref', the argument namevec must be provided")
	Fitted <- maref(object, namevec = namevec)
	return(Fitted)
}
 
}

residuals.rq.counts <- function(object, ...){

ans <- drop(object$y) - predict(object, type = "response", ...)
return(ans)

}

print.rq.counts <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    tau <- x$tau
    nq <- length(tau)
    cat("Call: ")
    dput(x$call)
    cat("\n")
    if (nq == 1) {
        cat(paste("Quantile", tau, "\n"))
        cat("\n")
        cat("Fixed effects:\n")
        printCoefmat(x$tTable, signif.stars = TRUE, P.values = TRUE)
    }
    else {
    NULL
	}
}

maref.rq.counts <- function(object, namevec){

tau <- object$tau
nq <- length(tau)

tsf <- object$tsf
symm <- attributes(tsf)$symm
dbounded <- attributes(tsf)$dbounded

betahat <- as.matrix(object$coefficients)
lambda <- object$lambda

# Model frames and matrices
mt <- terms(object)
all_vars <- get_all_vars(delete.response(mt), object$data)

# Work out expression of linear predictor for symbolic derivative
f <- terms2expr(mt)
var_labels <- all.vars(mt)[-1]
g <- parse(text = paste0("function(", paste(var_labels, collapse = ", "), ", ", paste(attr(f, "coefs"), collapse = ", "), "){}"))
d2 <- deriv(expr = f, namevec = as.character(namevec), function.arg = eval(g))
cat("The linear component of the marginal effect is calculated as derivative of", "\n", deparse(f), "\n with respect to", namevec, "\n")

# 

linpred <- object$x %*% betahat
n <- nrow(object$x)
dlinpred <- matrix(NA, n, nq)

for(j in 1:nq){
        argsLs <- as.list(betahat[attr(f, "labels"),j])
        names(argsLs) <- attr(f, "coefs")
        argsLs <- c(as.list(all_vars), argsLs)
        D <- do.call(d2, args = argsLs)
        dlinpred[,j] <- attr(D, "gradient")
}

val <- matrix(NA, n, nq)
for(j in 1:nq)(
val[,j] <- switch(tsf,
        mcjI = d1mcjI(linpred[,j], lambda, symm, dbounded),
        bc = d1bc(linpred[,j], lambda)
        )*dlinpred[,j]
)

return(val)

}

addnoise <- function(x, centered = TRUE, B = 0.999) 
{

	n <- length(x)
    if (centered) 
        z <- x + runif(n, -B/2, B/2)
    else z <- x + runif(n, 0, B)
	
    return(z)
}
