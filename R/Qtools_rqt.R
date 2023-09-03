######################################################################
### Transformation models
######################################################################

# base transformation functions

powerbase <- function(x, lambda){
(x^(lambda) - 1)/lambda
}

invpowerbase <- function(x, lambda, replace = TRUE){
sx <- if(replace) 0 else NA
val <- (lambda*x + 1)^(1/lambda)
val[(lambda*x + 1) <= 0] <- sx
return(val)
}

powrecbase <- function(x, lambda){
1/(2*lambda) * (x^lambda - x^(-lambda))
}

# Logit transformation
logit <- function(theta, omega = 0.001){
if(any(theta < 0) | any(theta > 1)) stop("theta outside interval")
theta[theta == 0] <- omega
theta[theta == 1] <- 1-omega
val <- log(theta/(1-theta))
return(val)
}

# Inverse logit transformation
invlogit <- function(x){
val <- exp(x)/(1 + exp(x))
val[val < 0] <- 0
val[val > 1] <- 1
return(val)
}

# c-log-log transformation
cloglog <- function(theta, omega = 0.001){
if(any(theta < 0) | any(theta > 1)) stop("theta outside interval")
theta[theta == 0] <- omega
theta[theta == 1] <- 1-omega
val <- log(-log(1 - theta))
return(val)
}

# inverse c-log-log transformation
invcloglog <- function(x){
val <- 1 - exp(-exp(x))
return(val)
}

# Proposal I (one parameter)

mcjI <- function(x, lambda, symm = TRUE, dbounded = FALSE, omega = 0.001){

if(dbounded){
	if(any(x < 0) | any(x > 1)) stop("x outside interval")
	x[x == 0] <- omega
	x[x == 1] <- 1 - omega
} else {
	if(any(x <= 0)) stop("x must be strictly positive")
}

if(dbounded){
	if(symm){
		x <- x/(1-x)
	} else {
		x <- -log(1-x)
	}
} else {
	if(!symm){
		x <- log(1+x)
	}
}

if(lambda != 0){
	val <-  powrecbase(x, lambda)
	} else {val <- log(x)}

return(val)

}

# Inverse proposal I (one parameter)

invmcjI <- function(x, lambda, symm = TRUE, dbounded = FALSE){

if(dbounded){
	if(symm){
		if(lambda != 0){
				x <- lambda*x
				y <- (x + sqrt(1 + x^2))^(1/lambda)
				val <- y/(1+y)
			} else {
			val <- invlogit(x)
		}
	} else {
		if(lambda != 0){
				x <- lambda*x
				val <- (x + sqrt(1 + x^2))^(1/lambda)
				val <- 1 - exp(-val)
			} else {
			val <- invcloglog(x)
		}
	}
}
else {
	if(lambda != 0){
		x <- lambda*x
		val <- (x + sqrt(1 + x^2))^(1/lambda)
	} else {val <- exp(x)}

	if(!symm){
		val <- exp(val) - 1
	}
}

return(val)
}

# Proposal II (two parameters)

mcjII <- function(x, lambda, delta, dbounded = FALSE, omega = 0.001){

if(dbounded){
	if(any(x < 0) | any(x > 1)) stop("x outside interval")
	x[x == 0] <- omega
	x[x == 1] <- 1 - omega
} else {
	if(any(x <= 0)) stop("x must be strictly positive")
}

if(lambda == 0){
	lambda <- 1e-10
}
if(delta == 0){
	delta <- 1e-10
}

if(dbounded){
	x <- ((1-x)^(-delta) - 1)/delta
	val <- powrecbase(x, lambda)
}
else {
	x <- x/(1+x)
	x <- ((1-x)^(-delta) - 1)/delta
	val <- powrecbase(x, lambda)
}


return(val)

}

# Inverse proposal II (two parameters)

invmcjII <- function(x, lambda, delta, dbounded = FALSE){

if(lambda == 0){
	lambda <- 1e-10
}
if(delta == 0){
	delta <- 1e-10
}

if(dbounded){
	x <- lambda*x
	val <- (x + sqrt(1 + x^2))^(1/lambda)
	val <- 1 - (val*delta + 1)^(-1/delta)
}
else {
	x <- lambda*x
	val <- (x + sqrt(1 + x^2))^(1/lambda)
	val <- 1 - (val*delta + 1)^(-1/delta)
	val <- val/(1-val)
}

return(val)
}

# Aranda-Ordaz transformation (symmetric and asymmetric)

ao <- function(theta, lambda, symm = TRUE, omega = 0.001){
if(any(theta < 0) | any(theta > 1)) stop("theta outside interval")
theta[theta == 0] <- omega
theta[theta == 1] <- 1 - omega

if(symm){
	if(lambda != 0){
	val <- (2/lambda)* ((theta^lambda) - (1-theta)^lambda)/((theta^lambda) + (1-theta)^lambda)
	} else {
	val <- logit(theta)
	}
} else {
	if(lambda != 0){
	val <- log(((1-theta)^(-lambda) - 1)/lambda)
	} else {
	val <- cloglog(theta)
	}
}


return(val)

}

# Inverse Aranda-Ordaz transformation (symmetric and asymmetric)

invao <- function(x, lambda, symm  = TRUE, replace = TRUE){
sx <- if(replace) 0 else NA
dx <- if(replace) 1 else NA

if(symm){
	if(lambda != 0){
		y <- (lambda*x/2)
		a <- (1 + y)^(1/lambda)
		b <- (1 - y)^(1/lambda)
		val <- rep(dx, length(x))
		val <- ifelse(abs(y) < 1, a/(a + b), val)
		val[y <= -1] <- sx
	} else {val <- invlogit(x)}
} else {
	if(lambda != 0){
		y <- lambda*exp(x)
		val <- ifelse(y > -1, 1 - (1 + y)^(-1/lambda), 1)
		
	} else {val <- invcloglog(x)}
}

return(as.numeric(val))
}

# Box-Cox transformation

bc <- function(x, lambda){
if(any(x <= 0)) stop("x must be strictly positive")
val <- if(lambda != 0) powerbase(x, lambda) else log(x)
return(val)
}

# Inverse Box-Cox transformation

invbc <- function(x, lambda, replace = TRUE){

val <- if(lambda != 0) invpowerbase(x, lambda, replace = replace) else exp(x)
return(val)
}

# Mapping from x.r[1],x.r[2] to 0,1

map <- function(x, x.r = NULL){
if(is.null(x.r)) x.r <- range(x, na.rm = TRUE)
theta <- (x - x.r[1])/(x.r[2] - x.r[1])
attr(theta, "range") <- x.r
return(theta)
}

# Mapping from 0,1 to x.r[1],x.r[2]

invmap <- function(x, x.r = NULL){

if(is.null(x.r)) x.r <- attr(x, "range")

(x.r[2] - x.r[1]) * x + x.r[1]

}

# L1-norm and residual cusum loss functions
l1Loss <- function(x, tau, weights){

ind <- ifelse(x < 0, 1, 0)
sum(weights * x * (tau - ind))/sum(weights)

}

rcLoss <- function(lambda, x, y, weights, tsf, symm = TRUE, dbounded = FALSE, tau = 0.5, method.rq = "fn"){

if(length(tau) > 1) stop("One quantile at a time")
n <- length(y)
out <- rep(NA, n)

z <- switch(tsf,
	mcjI = mcjI(y, lambda, symm, dbounded, omega = 0.001),
	bc = bc(y, lambda),
	ao = ao(y, lambda, symm, omega = 0.001)
	)

Rfun <- function(x, t, e) mean(apply(x, 1, function(xj,t) all(xj <= t), t = t) * e)

fit <- try(rq.wfit(x, z, tau = tau, weights = weights, method = method.rq), silent = T)

if(!inherits(fit, "try-error")){
	e <- as.numeric(fit$residuals <= 0)
	#out <- apply(x, 1, function(t, z, e) Rfun(z, t, e), z = x, e = tau - e)
	for(i in 1:n){
	#out[i] <- mean(apply(x, 1, function(x,t) all(x <= t), t = x[i,]) * (tau - e))
	out[i] <- mean(apply(t(x) <= x[i,], 2, function(x) all(x)) * (tau - e))	
	}
}

return(mean(out^2))

}

nlLoss <- function(theta, x, y, tau, tsf, symm = TRUE, dbounded = FALSE, smooth = FALSE, omicron = 0.001) {

if(any(is.na(theta))){
	ans <- Inf
	attr(ans, "grad") <- matrix(Inf, ncol(x))
	return(ans)
}

n <- length(y)
eta <- x %*% matrix(theta[-1])
res <- switch(tsf,
	bc = y - invbc(eta, lambda = theta[1]),
	ao = y - invao(eta, lambda = theta[1], symm = symm),
	mcjI = y - invmcjI(eta, lambda = theta[1], symm = symm, dbounded = dbounded)
)

d1 <- switch(tsf,
	bc = d1bc(eta, lambda = theta[1]),
	ao = d1ao(eta, lambda = theta[1], symm = symm),
	mcjI = d1mcjI(eta, lambda = theta[1], symm = symm, dbounded = dbounded)
)

d2 <- switch(tsf,
	bc = d2bc(eta, lambda = theta[1]),
	ao = d2ao(eta, lambda = theta[1], symm = symm),
	mcjI = d2mcjI(eta, lambda = theta[1], symm = symm, dbounded = dbounded)
)

if(smooth){
	s <- ifelse(res <= (tau - 1)*omicron, -1, ifelse(res >= tau*omicron, 1, 0))
	w <- as.numeric(1 - s^2)
	W <- diag(w, n, n)
	vs <- s*((2*tau - 1)*s + 1)/2
	cs <- sum(0.25*(1-2*tau)*omicron*s - 0.25*(1-2*tau+2*tau^2)*omicron*s^2)
	res <- matrix(res)
	ans <- as.numeric(0.5 * omicron * t(res) %*% W %*% res + t(vs) %*% res + cs)
	gradl <- -sum(1/omicron * W %*% (res * d2) + (vs * d2))
	gradb <- -t(x) %*% matrix(1/omicron * W %*% (res * d1) + (vs * d1))
	grad <- c(gradl, gradb)
	#hess <- 1/omicron * t(x) %*% (W * d1) %*% x
} else {
	ind <- tau - as.numeric(res < 0)
	ans <- as.numeric(sum(res*ind))
	grad <- c(-sum(ind*d2), -t(x) %*% (ind * d1))
	#hess <- matrix(0, ncol(x), ncol(x))
}

if(is.na(ans)){
	ans <- Inf
	grad <- matrix(Inf, ncol(x))
}

attr(ans, "grad") <- matrix(grad)

return(ans)
}

######################################################################
# One-parameter transformations (MCJI, Box-Cox, Aranda-Ordaz)
######################################################################

# Two-stage estimator

tsrq <- function(formula, data = sys.frame(sys.parent()), tsf = "mcjI", symm = TRUE, dbounded = FALSE, lambda = NULL, conditional = FALSE, tau = 0.5, subset, weights, na.action, contrasts = NULL, method = "fn"){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")

nq <- length(tau)

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
if (method == "model.frame") 
	return(mf)
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
x <- model.matrix(mt, mf, contrasts)
y <-  y.old <- model.response(mf, "numeric")
w <- if (missing(weights)) rep(1, length(y)) else as.vector(model.weights(mf))

if(tsf == "mcjII") stop("For two-parameter transformations, see 'tsrq2' and 'nlrq2'")
if(!tsf %in% c("bc","ao","mcjI")) stop("'tsf' not recognized")

isBounded <- (tsf == "mcjI" && dbounded)
isBounded <- tsf == "ao" || isBounded
if(isBounded) y <- map(y)

if(is.null(lambda) && !conditional){
	lambda <- if(tsf == "ao" & symm == FALSE) seq(-2, 2, by = 0.005)
		else seq(0, 2, by = 0.005)
}
nl <- length(lambda)

n <- length(y)
p <- ncol(x)

zhat <- res <- array(NA, dim = c(n, nq, nl))
matLoss <- rejected <- matrix(NA, nq, nl)
Ind <- array(NA, dim = c(n, nq, nl))

if(!conditional){
	# estimate linear QR for a sequence of lambdas
	for(i in 1:nl){

	# transform response
	newresponse <- switch(tsf,
		mcjI = mcjI(y, lambda[i], symm, dbounded, omega = 0.001),
		bc = bc(y, lambda[i]),
		ao = ao(y, lambda[i], symm, omega = 0.001)
		)

	# estimate linear QR for different taus
		for(j in 1:nq){
			fit <- try(do.call(rq.wfit, args = list(x = x, y = newresponse, tau = tau[j], weights = w, method = method)), silent = TRUE)
			if(!inherits(fit, "try-error")){
			zhat[,j,i] <- drop(x %*% fit$coefficients)
			
			Fitted <- switch(tsf,
				mcjI = invmcjI(zhat[,j,i], lambda[i], symm, dbounded),
				bc = invbc(zhat[,j,i], lambda[i]),
				ao = invao(zhat[,j,i], lambda[i], symm),
			)
			
			res <- y - Fitted
			
			if(tsf == "bc"){
				FLAG <- lambda[i]*zhat[,j,i] + 1 > 0
				Ind[,j,i] <- FLAG
				rejected[j,i] <- mean(!FLAG)
				}

			if(tsf == "ao" & symm == TRUE){
				FLAG <- abs(lambda[i]*zhat[,j,i]/2) - 1 < 0
				Ind[,j,i] <- FLAG
				rejected[j,i] <- mean(!FLAG)
				}
			
			matLoss[j,i] <- l1Loss(res, tau = tau[j], weights = w)
			}
		}
	}
	if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))
	# minimise for lambda
	lambdahat <- apply(matLoss, 1, function(x, lambda) lambda[which.min(x)], lambda = lambda)
} else {
	if(is.null(lambda)) stop("Must specify value for 'lambda' when 'conditional = TRUE'")
	if(length(lambda) != nq) stop("Length of 'lambda' must be the same as length of 'tau'")
	lambdahat <- lambda
}

betahat <- matrix(NA, p, nq)
colnames(betahat) <- tau
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau
fit <- list()
Rho <- function(u, tau) u * (tau - (u < 0))

for(j in 1:nq){
	# transform response with optimal lambda
	newresponse <- switch(tsf,
		mcjI = mcjI(y, lambdahat[j], symm, dbounded, omega = 0.001),
		bc = bc(y, lambdahat[j]),
		ao = ao(y, lambdahat[j], symm, omega = 0.001)
		)
	# fit final model
	z <- do.call(rq.wfit, args = list(x = x, y = newresponse, tau = tau[j], weights = w, method = method))

	betahat[,j] <- coefficients(z)
	tmp <- z$fitted.values
	Fitted[,j] <- switch(tsf,
		mcjI = invmcjI(tmp, lambdahat[j], symm, dbounded),
		bc = invbc(tmp, lambdahat[j]),
		ao = invao(tmp, lambdahat[j], symm)
	)
	class(z) <- "rq"
	z$na.action <- attr(mf, "na.action")
	z$formula <- stats::update(formula, newresponse ~ .)
	z$terms <- mt
	z$xlevels <- .getXlevels(mt, mf)
	z$call <- call
	z$tau <- tau[j]
	z$weights <- w
	z$residuals <- drop(z$residuals)
	z$rho <- sum(Rho(z$residuals, tau[j]))
	z$method <- method
	z$fitted.values <- drop(z$fitted.values)
	attr(z, "na.message") <- attr(m, "na.message")
	z$model <- mf
	fit[[j]] <- z
}

if(isBounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
names(lambdahat) <- paste("tau =", format(round(tau, 3)))

fit$call <- call
fit$method <- method
fit$mf <- mf
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("get_all_vars")
fit$data <- eval(mf, parent.frame())
fit$y <- y.old
if(isBounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$lambda <- lambdahat
fit$lambda.grid <- if(conditional) NULL else lambda
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
attr(fit$tsf, "dbounded") <- dbounded
attr(fit$tsf, "isBounded") <- isBounded
attr(fit$tsf, "npar") <- 1
attr(fit$tsf, "conditional") <- conditional
fit$objective <- matLoss
fit$optimum <- apply(matLoss, 1, function(x) x[which.min(x)])
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$rejected <- rejected
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 1
fit$fn <- "tsrq"
class(fit) <- "rqt"
return(fit)
}


# Cusum process estimator

rcrq <- function(formula, data = sys.frame(sys.parent()), tsf = "mcjI", symm = TRUE, dbounded = FALSE, lambda = NULL, tau = 0.5, subset, weights, na.action, contrasts = NULL, method = "fn"){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")

nq <- length(tau)

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
if (method == "model.frame") 
	return(mf)
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
x <- model.matrix(mt, mf, contrasts)
y <-  y.old <- model.response(mf, "numeric")
w <- if (missing(weights)) rep(1, length(y)) else as.vector(model.weights(mf))

if(tsf == "mcjII") stop("For two-parameter transformations, see 'tsrq2' and 'nlrq2'")
if(!tsf %in% c("bc","ao","mcjI")) stop("'tsf' not recognized")

isBounded <- (tsf == "mcjI" && dbounded)
isBounded <- tsf == "ao" || isBounded
if(isBounded) y <- map(y)

if(is.null(lambda)){
	lambda <- if(tsf == "ao" & symm == FALSE) seq(-2, 2, by = 0.05)
		else seq(0, 2, by = 0.05)
}

n <- length(y)
p <- ncol(x)
nl <- length(lambda)

matLoss <- rejected <- matrix(NA, nq, nl)

for(i in 1:nl){


# estimate linear QR for for sequence of lambdas

	for(j in 1:nq){
	matLoss[j,i] <- rcLoss(lambda[i], x, y, w, tsf, symm = symm, dbounded = dbounded, tau = tau[j], method.rq = method)
	}

}

if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))

# minimise for lambda
lambdahat <- apply(matLoss, 1, function(x, lambda) lambda[which.min(x)], lambda = lambda)


betahat <- matrix(NA, p, nq)
colnames(betahat) <- tau
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau
fit <- list()
Rho <- function(u, tau) u * (tau - (u < 0))

for(j in 1:nq){
	# transform response with optimal lambda
	newresponse <- switch(tsf,
		mcjI = mcjI(y, lambdahat[j], symm, dbounded, omega = 0.001),
		bc = bc(y, lambdahat[j]),
		ao = ao(y, lambdahat[j], symm, omega = 0.001)
		)
	# fit final model
	z <- do.call(rq.wfit, args = list(x = x, y = newresponse, tau = tau[j], weights = w, method = method))

	betahat[,j] <- coefficients(z)
	tmp <- z$fitted.values
	Fitted[,j] <- switch(tsf,
		mcjI = invmcjI(tmp, lambdahat[j], symm, dbounded),
		bc = invbc(tmp, lambdahat[j]),
		ao = invao(tmp, lambdahat[j], symm)
	)
	class(z) <- "rq"
	z$na.action <- attr(mf, "na.action")
	z$formula <- stats::update(formula, newresponse ~ .)
	z$terms <- mt
	z$xlevels <- .getXlevels(mt, mf)
	z$call <- call
	z$tau <- tau[j]
	z$weights <- w
	z$residuals <- drop(z$residuals)
	z$rho <- sum(Rho(z$residuals, tau[j]))
	z$method <- method
	z$fitted.values <- drop(z$fitted.values)
	attr(z, "na.message") <- attr(m, "na.message")
	z$model <- mf
	fit[[j]] <- z
}

if(isBounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
names(lambdahat) <- paste("tau =", format(round(tau, 3)))

fit$call <- call
fit$method <- method
fit$mf <- mf
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("get_all_vars")
fit$data <- eval(mf, parent.frame())
fit$y <- y.old
if(isBounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$lambda <- lambdahat
fit$lambda.grid <- lambda
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
attr(fit$tsf, "dbounded") <- dbounded
attr(fit$tsf, "isBounded") <- isBounded
attr(fit$tsf, "npar") <- 1
attr(fit$tsf, "conditional") <- FALSE
fit$objective <- matLoss
fit$optimum <- apply(matLoss, 1, function(x) x[which.min(x)])
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$rejected <- rejected
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 1
fit$fn <- "rcrq"
class(fit) <- "rqt"
return(fit)

}

# Nonlinear estimator

switch_check <- function(l0, l1, tol_ll, t0, t1, tol_theta, rule = "1"){
    deltal <- abs(l1/l0 - 1)
    deltat <- max(abs(t1/t0 - 1))
    switch(rule, `1` = deltal < tol_ll, `2` = deltal < tol_ll && 
        deltat < tol_theta)
}

nlControl <- function(tol_ll = 1e-05, tol_theta = 0.001, check_theta = FALSE, step = NULL, beta = 0.5, gamma = 1.25, reset_step = FALSE, maxit = 1000, smooth = FALSE, omicron = 0.001, verbose = FALSE) 
{
    if (beta > 1 || beta < 0) 
        stop("Beta must be a decreasing factor in (0,1)")
    if (gamma < 1) 
        stop("Beta must be a nondecreasing factor >= 1")
    if (maxit < 0) 
        stop("Number of iterations cannot be negative")
    list(tol_ll = tol_ll, tol_theta = tol_theta, 
        check_theta = check_theta, step = step, beta = beta, 
        gamma = gamma, reset_step = reset_step, maxit = as.integer(maxit),
        smooth = smooth, omicron = omicron, verbose = verbose)
}

nl.fit.rqt <- function(theta, x, y, tau, tsf, symm = TRUE, dbounded = FALSE, control){

step <- control$step
maxit <- control$maxit

theta_0 <- theta
ll_0 <- nlLoss(theta = theta_0, x = x, y = y, tau = tau, tsf = tsf, symm = symm, dbounded = dbounded, smooth = control$smooth, omicron = control$omicron)
eps <- .Machine$double.eps

for(i in 1:maxit) {
	if(control$verbose) cat(paste0("  (", i, ") logLik = ", round(ll_0,12), "\n"))
	# line search
	theta_1 <- theta_0 - attributes(ll_0)$grad*step
	ll_1 <- nlLoss(theta = theta_1, x = x, y = y, tau = tau, tsf = tsf, symm = symm, dbounded = dbounded, smooth = control$smooth, omicron = control$omicron)

	if(ll_1 > ll_0){
		if(control$verbose) cat("  Decreasing step...\n")
		step <- step*control$beta
	} else {
		rule <- if(control$check_theta) "2" else "1"
		check <- switch_check(ll_0, ll_1, control$tol_ll, theta_0, theta_1, control$tol_theta, rule = rule)
		if(check) break
		theta_0 <- theta_1
		ll_0 <- ll_1
		step <- if(control$reset_step) control$loop_step else step*control$gamma
	} 
}

list(par = as.numeric(theta_1), grad = attributes(ll_1)$grad, optimum = as.numeric(ll_1), CONVERGE = if(i==maxit) -1 else i)

}

nlrq1 <- function(formula, data = sys.frame(sys.parent()), tsf = "mcjI", symm = TRUE, dbounded = FALSE, start = NULL, tau = 0.5, subset, weights, na.action, contrasts = NULL, control = list()){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")

nq <- length(tau)

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
x <- model.matrix(mt, mf, contrasts)
y <-  y.old <- model.response(mf, "numeric")
w <- if (missing(weights)) rep(1, length(y)) else as.vector(model.weights(mf))

if(tsf == "mcjII") stop("For two-parameter transformations, see 'tsrq2' and 'nlrq2'")
if(!tsf %in% c("bc","ao","mcjI")) stop("'tsf' not recognized")

isBounded <- (tsf == "mcjI" && dbounded)
isBounded <- tsf == "ao" || isBounded
if(isBounded) y <- map(y)

n <- length(y)
p <- ncol(x)

# starting values
lambda_0 <- if(is.null(start)) 0 else start[1]

newresponse <- switch(tsf,
	mcjI = mcjI(y, lambda = lambda_0, symm = symm, dbounded = dbounded, omega = 0.001),
	bc = bc(y, lambda = lambda_0),
	ao = ao(y, lambda = lambda_0, symm = symm, omega = 0.001)
)

if(is.null(start)){
	start <- c(lambda_0, rq.fit(x, newresponse, tau = 0.5)$coef)
} else {
	if(length(start) != (p + 1)) stop("Length of vector of starting values does not match number of parameters")
}


if (is.null(names(control))) 
	control <- nlControl()
else {
	control_default <- nlControl()
	control_names <- intersect(names(control), names(control_default))
	control_default[control_names] <- control[control_names]
	control <- control_default
}
if (is.null(control$step)) 
	control$step <- sd(as.numeric(y))

# estimate nonlinear QR using gradient search

fit <- list()
betahat <- matrix(NA, p, nq)
lambdahat <- rep(NA, nq)
Fitted <- matrix(NA, n, nq)

for(j in 1:nq){
	fit[[j]] <- try(nl.fit.rqt(theta = start, x = x, y = y, tau = tau[j], tsf = tsf, symm = symm, dbounded = dbounded, control = control))
	if(!inherits(fit[[j]],"try-error")){
		tmp <- x%*%matrix(fit[[j]]$par[-1])
		betahat[,j] <- fit[[j]]$par[-1]
		lambdahat[j] <- fit[[j]]$par[1]
		if(!is.na(lambdahat[j])){
		Fitted[,j] <- switch(tsf,
			mcjI = invmcjI(tmp, lambdahat[j], symm, dbounded),
			bc = invbc(tmp, lambdahat[j]),
			ao = invao(tmp, lambdahat[j], symm)
		)}
	}
}

if(isBounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
names(lambdahat) <- paste("tau =", format(round(tau, 3)))

fit$call <- call
fit$method <- "gradient-search"
fit$mf <- mf
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("get_all_vars")
fit$data <- eval(mf, parent.frame())
fit$y <- y.old
if(isBounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$lambda <- lambdahat
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
attr(fit$tsf, "dbounded") <- dbounded
attr(fit$tsf, "isBounded") <- isBounded
attr(fit$tsf, "npar") <- 1
attr(fit$tsf, "conditional") <- FALSE
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 1
fit$fn <- "nlrq1"
fit$control <- control
class(fit) <- "rqt"
return(fit)
}

######################################################################
# Two-parameter transformations (MCJII)
######################################################################

# Two-stage estimator

tsrq2 <- function(formula, data = sys.frame(sys.parent()), dbounded = FALSE, lambda = NULL, delta = NULL, conditional = FALSE, tau = 0.5, subset, weights, na.action, contrasts = NULL, method = "fn"){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), 
	names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval.parent(mf)
if (method == "model.frame") 
	return(mf)
mt <- attr(mf, "terms")
x <- model.matrix(mt, mf, contrasts)
y <- y.old <- model.response(mf)
w <- if (missing(weights)) rep(1, length(y)) else as.vector(model.weights(mf))
if(dbounded) y <- map(y)

if(is.null(lambda) && !conditional){
	lambda <- seq(0, 2, by = 0.005)
}
if(is.null(delta) && !conditional){
	delta <- seq(0, 2, by = 0.005)
}
nl <- length(lambda)
nd <- length(delta)

n <- length(y)
p <- ncol(x)
nq <- length(tau)

matLoss <- array(NA, dim = c(nl, nd, nq), dimnames = list(lambda = 1:nl, delta = 1:nd, tau = tau))

if(!conditional){
	for(k in 1:nd){
		for(i in 1:nl){
		# transform response
		newresponse <- mcjII(y, lambda[i], delta[k], dbounded, omega = 0.001)
			for(j in 1:nq){
			fit <- try(do.call(rq.wfit, args = list(x = x, y = newresponse, tau = tau[j], weights = w, method = method)), silent = TRUE)
				
				if(!inherits(fit, "try-error")){
				Fitted <- invmcjII(drop(x %*% fit$coefficients), lambda[i], delta[k], dbounded)
				matLoss[i,k,j] <- l1Loss(y - Fitted, tau = tau[j], weights = w)
				}
			}
		}
	}
	if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))

	# minimise for lambda
	parhat <- apply(matLoss, 3, function(x, lambda, delta){
	m <- which(x == min(x, na.rm = TRUE), arr.ind = TRUE)[1,];
	return(c(lambda[m[1]], delta[m[2]]))}, lambda = lambda, delta = delta)
} else {
	if(is.null(lambda)) stop("Must specify value for 'lambda' when 'conditional = TRUE'")
	if(length(lambda) != nq) stop("Length of 'lambda' must be the same as length of 'tau'")
	if(is.null(delta)) stop("Must specify value for 'delta' when 'conditional = TRUE'")
	if(length(delta) != nq) stop("Length of 'delta' must be the same as length of 'tau'")
	parhat <- matrix(c(lambda, delta), ncol = nq, byrow = TRUE)
}

betahat <- matrix(NA, p, nq)
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau
fit <- list()
Rho <- function(u, tau) u * (tau - (u < 0))

for(j in 1:nq){
	# transform response with optimal lambda
	newresponse <- mcjII(y, parhat[1,j], parhat[2,j], dbounded, omega = 0.001)
	
	# fit final model
	z <- do.call(rq.wfit, args = list(x = x, y = newresponse, tau = tau[j], weights = w, method = method))

	betahat[,j] <- coefficients(z)
	tmp <- z$fitted.values
	Fitted[,j] <- invmcjII(z$fitted.values, parhat[1,j], parhat[2,j], dbounded)
	class(z) <- "rq"
	z$na.action <- attr(mf, "na.action")
	z$formula <- stats::update(formula, newresponse ~ .)
	z$terms <- mt
	z$xlevels <- .getXlevels(mt, mf)
	z$call <- call
	z$tau <- tau[j]
	z$weights <- w
	z$residuals <- drop(z$residuals)
	z$rho <- sum(Rho(z$residuals, tau[j]))
	z$method <- method
	z$fitted.values <- drop(z$fitted.values)
	attr(z, "na.message") <- attr(m, "na.message")
	z$model <- mf
	fit[[j]] <- z
}

if(dbounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
dimnames(parhat) <- list(c("lambda","delta"), paste("tau =", format(round(tau, 3))))

fit$call <- call
fit$method <- method
fit$mf <- mf
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("get_all_vars")
fit$data <- eval(mf, parent.frame())
fit$y <- y.old
if(dbounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$eta <- parhat
fit$lambda.grid <- if(conditional) NULL else lambda
fit$delta.grid <- if(conditional) NULL else delta
fit$tsf <- "mcjII"
attr(fit$tsf, "dbounded") <- dbounded
attr(fit$tsf, "isBounded") <- dbounded
attr(fit$tsf, "npar") <- 2
attr(fit$tsf, "conditional") <- conditional
fit$objective <- matLoss
fit$optimum <- if(conditional) NA else apply(matLoss, 3, function(x){m <- which(x == min(x, na.rm = TRUE), arr.ind = TRUE)[1,];return(x[m[1],m[2]])})
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 2
fit$fn <- "tsrq2"
class(fit) <- "rqt"
return(fit)
}

# Nelder-Mead optimization (joint estimation)

nlrq2 <- function(formula, data = sys.frame(sys.parent()), dbounded = FALSE, start = NULL, tau = 0.5, subset, weights, na.action, contrasts = NULL){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), 
	names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval.parent(mf)
mt <- attr(mf, "terms")
x <- model.matrix(mt, mf, contrasts)
y <- y.old <- model.response(mf)
w <- if (missing(weights)) rep(1, length(y)) else as.vector(model.weights(mf))
if(dbounded) y <- map(y)

n <- length(y)
p <- ncol(x)
nq <- length(tau)

if(is.null(start)) start <- rep(0, p + 2)

f <- function(theta, dataLs){

	if(theta[2] < 0) return(Inf)
	Fitted <- invmcjII(dataLs$x %*% matrix(theta[-c(1:2)]), lambda = theta[1], delta = theta[2], dbounded = dataLs$dbounded)
	return(l1Loss(dataLs$y - Fitted, tau = dataLs$tau, weights = dataLs$weights))
}

fit <- list()
betahat <- matrix(NA, p, nq)
parhat <- matrix(NA, 2, nq)
Fitted <- matrix(NA, n, nq)

for(j in 1:nq){
	fit[[j]] <- try(optim(par = start, fn = f, method = "Nelder-Mead", dataLs = list(x = x, y = y, dbounded = dbounded, tau = tau[j], weights = w)), silent = T)

	if(!inherits(fit[[j]], "try-error")){
		betahat[,j] <- fit[[j]]$par[-c(1:2)]
		parhat[,j] <- c(fit[[j]]$par[1], fit[[j]]$par[2])
		Fitted[,j] <- invmcjII(x %*% matrix(betahat[,j]), parhat[1,j], parhat[2,j], dbounded)
	}
}

if(dbounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
dimnames(parhat) <- list(c("lambda","delta"), paste("tau =", format(round(tau, 3))))


fit$call <- call
fit$method <- "Nelder-Mead"
fit$mf <- mf
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("get_all_vars")
fit$data <- eval(mf, parent.frame())
fit$y <- y.old
if(dbounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$eta <- parhat
fit$tsf <- "mcjII"
attr(fit$tsf, "dbounded") <- dbounded
attr(fit$tsf, "isBounded") <- dbounded
attr(fit$tsf, "npar") <- 2
attr(fit$tsf, "conditional") <- FALSE
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 2
fit$fn <- "nlrq2"
class(fit) <- "rqt"
return(fit)
}

######################################################################
# Print, summary, bootstrap, predict, fitted for class rqt
######################################################################

print.rqt <- function(x, ...){

if (!is.null(cl <- x$call)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}
type <- if(attr(x$tsf, "isBounded")) "(doubly bounded response)" else "(singly bounded response)"
tsf <- switch(x$tsf,
	mcjI = "Proposal I",
	bc = "Box-Cox",
	ao = "Aranda-Ordaz",
	mcjII = "Proposal II")
if(x$tsf %in% c("mcjI", "ao")){
	tsf <- paste(tsf, if(attr(x$tsf, "symm"))
	"symmetric" else "asymmetric")
}
tsf <- paste(tsf, "transformation", type)

cat(tsf, "\n")
cat("\nOptimal transformation parameter:\n")
if(x$tsf == "mcjII") print(x$eta) else print(x$lambda)

coef <- x$coefficients
cat("\nCoefficients linear model (transformed scale):\n")
print(coef, ...)

nobs <- length(x$y)
p <- ncol(x$x)
rdf <- nobs - p
cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
if (!is.null(attr(x, "na.message"))) 
	cat(attr(x, "na.message"), "\n")
invisible(x)
}

predict.rqt <- function(object, newdata, na.action = na.pass, type = "response", namevec = NULL, ...){

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
symm <- attributes(tsf)$symm
dbounded <- attributes(tsf)$dbounded
isBounded <- attributes(tsf)$isBounded

if(tsf == "mcjII"){
	etahat <- object$eta
} else {
	lambdahat <- object$lambda
}

if(!missing(newdata)) {
	mt <- terms(object)
	Terms <- delete.response(mt)
	object$mf <- model.frame(Terms, newdata, na.action = na.action, xlev = object$levels) # model frame
	if (!is.null(cl <- attr(Terms, "dataClasses")))
		.checkMFClasses(cl, object$mf)
	object$x <- model.matrix(Terms, object$mf, contrasts.arg = object$contrasts) # model matrix
	object$data <- newdata # input variables
}

linpred <- object$x %*% object$coefficients

Fitted <- matrix(NA, nrow(linpred), ncol(linpred))
if(type == "link"){
	return(linpred)
}

if(type == "response"){
	if(tsf == "mcjII"){
		for(j in 1:nq){
			Fitted[,j] <- invmcjII(x = linpred[,j], lambda = etahat['lambda',j], delta = etahat['delta',j], dbounded = dbounded)
		}
	} else {
		for(j in 1:nq){
			Fitted[,j] <- switch(tsf,
				mcjI = invmcjI(linpred[,j], lambdahat[j], symm, dbounded),
				bc = invbc(linpred[,j], lambdahat[j]),
				ao = invao(linpred[,j], lambdahat[j], symm))
		}
	}

	if(isBounded){
		Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(object$y))
	}
return(Fitted)
}

if(type == "maref"){
	if(is.null(namevec)) stop("When type = 'maref', the argument namevec must be provided")
	if(tsf == "mcjII") stop("Marginal effects not available for tsf = 'mcjII'")
	Fitted <- maref(object, namevec = namevec)
	return(Fitted)

}

}

terms2expr <- function(object){

if(!inherits(object, "terms")) stop("Only objects of class 'terms'")

Irm <- function(x){
	n <- nchar(x)
	flag <- substr(x, 1, 2) == "I(" & substr(x, n, n) == ")"
	if(flag) x <- substr(x, 2, n)
	return(x)
}

x <- object
variables <- as.list(attr(x, "variables"))[-1]
mt <- attr(x, "term.labels")
mt <- sapply(mt, Irm)
mt <- sub(":", "*", mt)
term.labels <- names(mt)

#coefs <- term.labels
#coefs <- gsub(pattern = "\\(", replacement = "", x = coefs)
#coefs <- gsub(pattern = "\\)", replacement = "", x = coefs)
#coefs <- gsub(pattern = "[[:punct:]]", replacement = "", x = coefs)
#coefs <- gsub(pattern = "[[:space:]]", replacement = "", x = coefs)
#coefs <- gsub(pattern = ":", replacement = "_", x = coefs)
#coefs <- paste0("beta.", coefs)
coefs <- paste0("beta", 1:length(term.labels))

val <- as.formula(paste(variables[1], "~", paste(paste(coefs, mt, sep = "*"), collapse = " + ")))
attr(val, "terms") <- as.vector(mt)
attr(val, "labels") <- term.labels
attr(val, "coefs") <- coefs
return(val)

}

maref.rqt <- function(object, namevec){

tau <- object$tau
nq <- length(tau)

tsf <- object$tsf
symm <- attributes(tsf)$symm
dbounded <- attributes(tsf)$dbounded

betahat <- object$coefficients
lambdahat <- object$lambda

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
	mcjI = d1mcjI(linpred[,j], lambdahat[j], symm, dbounded),
	bc = d1bc(linpred[,j], lambdahat[j]),
	ao = d1ao(linpred[,j], lambdahat[j], symm),
	)*dlinpred[,j]
)

return(val)

}

boot.rqt <- function(data, inds, object){

tau <- object$tau
nq <- length(tau)
all.obs <- rownames(object$x)
n <- length(all.obs)
flag <- !object$tsf %in% c("mcjII")
nn <- if(flag) c(object$term.labels, "lambda") else c(object$term.labels, "lambda", "delta")

if(nq == 1){
	fit <- stats::update(object, data = data[inds,])
	val <- fit$coefficients
	val <- if(flag) c(val, fit$lambda) else c(val, fit$eta)
} else {
	fit <- stats::update(object, data = data[inds,])
	val <- fit$coefficients
	val <- if(flag) rbind(val, fit$lambda) else rbind(val, fit$eta)
}

val <- as.vector(val)
names(val) <- rep(nn, nq)
return(val)

}

summary.rqt <- function(object, alpha = 0.05, se = "boot", R = 50, sim = "ordinary", stype = "i", conditional = FALSE, ...){

call <- match.call(expand.dots = TRUE)

tau <- object$tau
nq <- length(tau)
mpar <- ncol(object$x)
ntot <- mpar + attr(object$tsf, "npar")
if(mpar == 1) object$mf$intercept <- 1

flag <- (!conditional) && (se %in% c("iid","nid"))

if(object$tsf == "mcjII" && flag)  stop("Summary not available. Change to 'se = boot'.")
if(object$fn == "rcrq" && flag) stop("Summary not available. Change to 'se = boot'.")
if(conditional && object$fn == "nlrq2") stop("Conditional inference not available for objects from 'nlrq2'. Change to 'conditional = FALSE'.")


if(attr(object$tsf, "conditional")){
	if(!conditional) warning("Main call 'conditional = TRUE'")
	conditional <- TRUE
}

if(conditional){
	ans <- B <- list()
	for(j in 1:nq){
		Args <- list()
		Args$object <- object[[j]]
		Args$se <- se
		if(se == "boot") Args$R <- R
		nn <- c("covariance","hs","bsmethod","mofn","iid")
		nn <- nn[pmatch(names(call), nn, duplicates.ok = FALSE)]
		nn <- nn[!is.na(nn)]
		if(length(nn) > 0) {tmp <- as.list(call[[nn]]); names(tmp) <- nn; Args <- c(Args, tmp)}
		tmp <- do.call(summary.rq, args = Args)
		ans[[j]] <- tmp$coefficients
		if(!is.null(tmp$B)) B[[j]] <- tmp$B
		if(object$tsf == "mcjII") {
			tmp <- matrix(NA, nrow = 2, ncol = ncol(ans[[j]]))
			tmp[,1] <- object$eta[,j]
			rownames(tmp) <- c("lambda","delta")
		} else {
			tmp <- matrix(c(object$lambda[j], rep(NA, ncol(ans[[j]]) - 1)), nrow = 1)
			rownames(tmp) <- "lambda"
		}
		ans[[j]] <- rbind(ans[[j]], tmp)
	}
	object$B <- B
} else {
	if(se == "boot"){
		Args <- list()
		Args$data <- object$mf
		Args$statistic <- boot.rqt
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
	} else if(se %in% c("iid", "nid")) {
		ans <- list()
		S <- se_rqt(object, se = se)
		for(j in 1:nq){
			val <- c(object[[j]]$coefficients, lambda = object$lambda[j])
			val <- cbind(val, sqrt(diag(S[,,j])), val - sqrt(diag(S[,,j]))*qnorm(1-alpha/2), val + sqrt(diag(S[,,j]))*qnorm(1-alpha/2))
			nn <- c("Value", "Std. Error", "Lower bound", "Upper bound")
			colnames(val) <- nn
			ans[[j]] <- val
		}
		names(ans) <- tau
	} else ans <- NULL	
}

attr(ans, "conditional") <- conditional
object$coefficients <- ans
object$call <- call
class(object) <- c("summary.rqt", class(object))
return(object)

}

print.summary.rqt <- function(x, ...){

if (!is.null(cl <- x$call)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}

tau <- x$tau
nq <- length(tau)
mpar <- ncol(x$x)

type <- if(attr(x$tsf, "isBounded")) "(doubly bounded response)" else "(singly bounded response)"
tsf <- switch(x$tsf,
	mcjI = "Proposal I",
	bc = "Box-Cox",
	ao = "Aranda-Ordaz",
	mcjII = "Proposal II")
if (x$tsf %in% c("mcjI", "ao")){ 
	tsf <- paste(tsf, if (attr(x$tsf, "symm")) 
		"symmetric"
	else "asymmetric")
}
tsf <- paste(tsf, "transformation", type)

cat(tsf, "\n")

conditional <- if(attr(x$coefficients, "conditional")) "conditional" else "unconditional"
cat("\nSummary for", conditional, "inference\n")

for(i in 1:nq){
cat("\ntau = ", tau[i], "\n")

cat("\nOptimal transformation parameter:\n")
print(x$coefficients[[i]][-c(1:mpar),], ...)

cat("\nCoefficients linear model (transformed scale):\n")
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

fitted.rqt <- function(object, ...){

return(object$fitted.values)

}

residuals.rqt <- function(object, ...){

return(object$y - object$fitted.values)

}

coef.rqt <- coefficients.rqt <- function(object, all = FALSE, ...){

if(!inherits(object, "rqt")) stop("Class 'rqt' only")

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
nn <- object$term.labels

if(all){
	nn <- if(tsf %in% "mcjII") c(nn, "lambda", "delta") else c(nn, "lambda")
}

tpar <- if(tsf %in% "mcjII") object$eta else object$lambda

if(!all){
	ans <- object$coefficients
} else {
	ans <- if(nq == 1) c(object$coefficients, tpar)
		else rbind(object$coefficients, tpar)
}

if(nq == 1){
	names(ans) <- nn
} else {
	rownames(ans) <- nn
	colnames(ans) <- paste("tau =", tau)
}

return(ans)

}

##################################################
### Asymptotics
##################################################

d1bc <- function(x, lambda){
zero <- rep(0, length(x))
g1 <- deriv(~ (lambda*x + 1)^(1/lambda), "x", func = function(x,lambda){})
g2 <- deriv(~ exp(x), "x", func = function(x){})
    if (lambda != 0) {
        val <- as.numeric(attributes(g1(x, lambda))$gradient)
		val <- ifelse(lambda * x + 1 > 0, val, zero)
    }
    else {
        val <- as.numeric(attributes(g2(x))$gradient)
    }
    return(val)
}

d2bc <- function(x, lambda){
g1 <- deriv(~ (lambda*x + 1)^(1/lambda), "lambda", func = function(lambda,x){})
    if (lambda != 0) {
        val <- as.numeric(attributes(g1(lambda,x))$gradient)
    }
    else {
        val <- as.numeric(attributes(g1(0.00001,x))$gradient)
    }
    return(val)
}

d1mcjI <- function(x, lambda, symm, dbounded){

if(dbounded){
	if(symm){
	g1 <- deriv(~ (lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda)/(1 + (lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda)), "x", func = function(x, lambda){})
	g2 <- deriv(~ exp(x)/(1+exp(x)), "x", func = function(x){})
		if(lambda != 0){
			val <- as.numeric(attributes(g1(x, lambda))$gradient)
		} else {
			val <- as.numeric(attributes(g2(x))$gradient)
		}
	} else {
	g1 <- deriv(~ 1 - exp(-(lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda)), "x", func = function(x, lambda){})
	g2 <- deriv(~ 1 - exp(-exp(x)), "x", func = function(x){})
		if(lambda != 0){
			val <- as.numeric(attributes(g1(x, lambda))$gradient)
		} else {
			val <- as.numeric(attributes(g2(x))$gradient)
		}
	}
} else {
	if(symm){
	g1 <- deriv(~ (lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda), "x", func = function(x, lambda){})
	g2 <- deriv(~ exp(x), "x", func = function(x){})
		if (lambda != 0) {
			val <- as.numeric(attributes(g1(x, lambda))$gradient)
		} else {
			val <- as.numeric(attributes(g2(x))$gradient)
		}
	} else {
	g1 <- deriv(~ exp((lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda)) - 1, "x", func = function(x, lambda){})
	g2 <- deriv(~ exp(exp(x)) - 1, "x", func = function(x){})
		if (lambda != 0) {
			val <- as.numeric(attributes(g1(x, lambda))$gradient)
		} else {
			val <- as.numeric(attributes(g2(x))$gradient)
		}
	}
}

return(val)

}

d2mcjI <- function(x, lambda, symm, dbounded){

if(dbounded){
	if(symm){
	g1 <- deriv(~ (lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda)/(1 + (lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda)), "lambda", func = function(lambda, x){})
		if(lambda != 0){
			val <- as.numeric(attributes(g1(lambda, x))$gradient)
		} else {
			val <- as.numeric(attributes(g1(0.00001,x))$gradient)
		}
	} else {
	g1 <- deriv(~ 1 - exp(-(lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda)), "lambda", func = function(lambda, x){})
		if(lambda != 0){
			val <- as.numeric(attributes(g1(lambda, x))$gradient)
		} else {
			val <- as.numeric(attributes(g1(0.00001,x))$gradient)
		}
	}
} else {
	if(symm){
	g1 <- deriv(~ (lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda), "lambda", func = function(lambda, x){})
		if (lambda != 0) {
			val <- as.numeric(attributes(g1(lambda, x))$gradient)
		} else {
			val <- as.numeric(attributes(g1(0.00001,x))$gradient)
		}
	} else {
	g1 <- deriv(~ exp((lambda*x + sqrt(1 + (lambda*x)^2))^(1/lambda)) - 1, "lambda", func = function(lambda, x){})
		if (lambda != 0) {
			val <- as.numeric(attributes(g1(lambda, x))$gradient)
		} else {
			val <- as.numeric(attributes(g1(0.00001,x))$gradient)
		}
	}
}

return(val)

}

d1ao <- function(x, lambda, symm){
zero <- rep(0, length(x))
if(symm){
g1 <- deriv(~ (1 + lambda*x/2)^(1/lambda)/((1 + lambda*x/2)^(1/lambda) + (1 - lambda*x/2)^(1/lambda)), "x", func = function(x, lambda){})
g2 <- deriv(~ exp(x)/(1+exp(x)), "x", func = function(x){})
	if(lambda != 0){
		val <- as.numeric(attributes(g1(x, lambda))$gradient)
		val <- ifelse(abs(lambda * x/2) < 1, val, zero)
	} else {
		val <- as.numeric(attributes(g2(x))$gradient)
	}
} else {
g1 <- deriv(~ 1 - (1 + lambda*exp(x))^(-1/lambda), "x", func = function(x, lambda){})
g2 <- deriv(~ 1 - exp(-exp(x)), "x", func = function(x){})
	if(lambda != 0){
		val <- as.numeric(attributes(g1(x, lambda))$gradient)
		val <- ifelse(lambda * exp(x) > -1, val, zero)
		} else {
		val <- as.numeric(attributes(g2(x))$gradient)
	}
}

return(val)

}

d2ao <- function(x, lambda, symm){
zero <- rep(0, length(x))
if(symm){
g1 <- deriv(~ (1 + lambda*x/2)^(1/lambda)/((1 + lambda*x/2)^(1/lambda) + (1 - lambda*x/2)^(1/lambda)), "lambda", func = function(lambda, x){})
	if(lambda != 0){
		val <- as.numeric(attributes(g1(lambda, x))$gradient)
		val <- ifelse(abs(lambda * x/2) < 1, val, zero)
	} else {
		val <- as.numeric(attributes(g1(0.00001,x))$gradient)
	}
} else {
g1 <- deriv(~ 1 - (1 + lambda*exp(x))^(-1/lambda), "lambda", func = function(lambda, x){})
	if(lambda != 0){
		val <- as.numeric(attributes(g1(lambda, x))$gradient)
		val <- ifelse(lambda * exp(x) > -1, val, zero)
		} else {
		val <- as.numeric(attributes(g1(0.00001,x))$gradient)
	}
}

return(val)

}

se_rqt <- function(object, se = "nid"){


tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
symm <- attributes(tsf)$symm
dbounded <- attributes(tsf)$dbounded
isBounded <- attributes(tsf)$isBounded
lambdahat <- object$lambda

betahat <- as.matrix(object$coefficients)
linpred <- predict(object, type = "link")
fis <- sparsity.rqt(object, se = se)$density # density of 'u = y - hinv(xb)' 

x <- object$x
n <- nrow(x)
p <- ncol(x)

g1 <- g2 <- matrix(NA, n, nq)
dbl <- matrix(NA, p, nq)
V <- array(NA, dim = c(p + 1, p + 1, nq))

for(j in 1:nq){
	g1[,j] <- switch(tsf,
		mcjI = d1mcjI(linpred[,j], lambdahat[j], symm, dbounded),
		bc = d1bc(linpred[,j], lambdahat[j]),
		ao = d1ao(linpred[,j],lambdahat[j], symm)
		)

	g2[,j] <- switch(tsf,
		mcjI = d2mcjI(linpred[,j], lambdahat[j], symm, dbounded),
		bc = d2bc(linpred[,j], lambdahat[j]),
		ao = d2ao(linpred[,j], lambdahat[j], symm)
		)
	f0 <- fis[,j]
	
	dbl[,j] <- - solve(crossprod(sqrt(f0 * g1[,j]) * x)/n) %*% matrix(colMeans((f0 * g2[,j] * x)))

	A <- rbind(cbind(diag(p), matrix(0, p, p), rep(0, p)),
		c(rep(0,p),dbl[,j],1))

	d2 <- cbind(g1[,j] * x, g2[,j])
	d <- cbind(x,d2)
	H <- A %*% (t(f0*d) %*% d2)/n
	Hinv <- try(chol2inv(chol(H)), silent = TRUE)
	if(inherits(Hinv, "try-error")) Hinv <- try(solve(H), silent = TRUE)
	if(inherits(Hinv, "try-error")){
		Hinv <- matrix(NA, p + 1, p + 1)
		warning("Singular 'H' matrix")
	}
	L <- tau[j] * (1 - tau[j]) * A %*% (crossprod(d)/n) %*% t(A)

	V[,,j] <- Hinv %*% L %*% t(Hinv)/n
}

return(V)

}

sparsity.rqt <- function(object, se = "nid", hs = TRUE){
    mt <- terms(object)
    m <- model.frame(object)
    y <- model.response(m)
    x <- model.matrix(mt, m, contrasts = object$contrasts)
    wt <- model.weights(object$model)
    taus <- object$tau
    nt <- length(taus)
    eps <- .Machine$double.eps^(2/3)

    vnames <- dimnames(x)[[2]]
    residm <- sweep(- predict(object, type = "response"), 1, y, FUN = "+")
    n <- length(y)
    p <- length(coefficients(object, all = TRUE))
    rdf <- n - p
    if (!is.null(wt)) {
        residm <- residm * wt
        x <- x * wt
        y <- y * wt
    }
    if (is.null(se)) {
		se <- "nid"
    }

spar <- dens <- matrix(NA, n, nt)
for(i in 1:nt){

tau <- taus[i]

    if (se == "iid") {
		resid <- residm[,i]
        pz <- sum(abs(resid) < eps)
        h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = hs)))
        ir <- (pz + 1):(h + pz + 1)
        ord.resid <- sort(resid[order(abs(resid))][ir])
        xt <- ir/(n - p)
        spar[,i] <- rq(ord.resid ~ xt)$coef[2]
        dens[,i] <- 1/spar[,i]
    }
    else if (se == "nid") {
        h <- bandwidth.rq(tau, n, hs = hs)
        if (tau + h > 1) 
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
            stop("tau - h < 0:  error in summary.rq")
		call <- getCall(object)
		call$tau <- tau + h
		bhi <- eval(call, attr(terms(object), ".Environment"), parent.frame())
		call$tau <- tau - h
        blo <- eval(call, attr(terms(object), ".Environment"), parent.frame())
        dyhat <- predict(bhi, type = "response") - predict(blo, type = "response")
        if (any(dyhat <= 0)) 
            warning(paste(sum(dyhat <= 0), "non-positive fis"))
        f <- pmax(0, (2 * h)/(dyhat - eps))
        dens[,i] <- f
		spar[,i] <- 1/f
    }
    else if (se == "ker") {
        h <- bandwidth.rq(tau, n, hs = hs)
        if (tau + h > 1) 
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
            stop("tau - h < 0:  error in summary.rq")
        uhat <- c(residm[,i])
        h <- (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)), 
            (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
        f <- dnorm(uhat/h)/h
        dens[,i] <- f
		spar[,i] <- 1/f
    }
}# loop i

	colnames(dens) <- colnames(spar) <- taus
    return(list(density = dens, sparsity = spar, bandwidth = h))
}

sparsity.rq <- sparsity.rqs <-function(object, se = "nid", hs = TRUE){
    mt <- terms(object)
    m <- model.frame(object)
    y <- model.response(m)
    x <- model.matrix(mt, m, contrasts = object$contrasts)
    wt <- model.weights(object$model)
    taus <- object$tau
    nq <- length(taus)
    eps <- .Machine$double.eps^(2/3)

    vnames <- dimnames(x)[[2]]
    residm <- as.matrix(object$residuals)
    n <- length(y)
    p <- nrow(as.matrix(object$coef))
    rdf <- n - p
    if (!is.null(wt)) {
        residm <- residm * wt
        x <- x * wt
        y <- y * wt
    }
    if (is.null(se)) {
		se <- "nid"
    }

spar <- dens <- matrix(NA, n, nq)
for(i in 1:nq){

tau <- taus[i]

    if (se == "iid") {
		resid <- residm[,i]
		pz <- sum(abs(resid) < eps)
        h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = hs)))
        ir <- (pz + 1):(h + pz + 1)
        ord.resid <- sort(resid[order(abs(resid))][ir])
        xt <- ir/(n - p)
        spar[,i] <- rq(ord.resid ~ xt)$coef[2]
        dens[,i] <- 1/spar[,i]
    }
    else if (se == "nid") {
        h <- bandwidth.rq(tau, n, hs = hs)
        if (tau + h > 1) 
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
            stop("tau - h < 0:  error in summary.rq")
        bhi <- rq.fit.fnb(x, y, tau = tau + h)$coef
        blo <- rq.fit.fnb(x, y, tau = tau - h)$coef
        dyhat <- x %*% (bhi - blo)
        if (any(dyhat <= 0)) 
            warning(paste(sum(dyhat <= 0), "non-positive fis"))
        f <- pmax(0, (2 * h)/(dyhat - eps))
        dens[,i] <- f
	  spar[,i] <- 1/f
    }
    else if (se == "ker") {
        h <- bandwidth.rq(tau, n, hs = hs)
        if (tau + h > 1) 
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
            stop("tau - h < 0:  error in summary.rq")
        uhat <- c(residm[,i])
        h <- (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)), 
            (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
        f <- dnorm(uhat/h)/h
        dens[,i] <- f
	  spar[,i] <- 1/f
    }
}# loop i

	colnames(dens) <- colnames(spar) <- taus
    return(list(density = dens, sparsity = spar, bandwidth = h))
}

