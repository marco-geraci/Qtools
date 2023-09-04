##################################################
# Quantile ratio regression
##################################################

qrr <- function(formula, data, taus, start = "rq", tsf = "bc", symm = TRUE, dbounded = FALSE, linearize = TRUE, kernel = "Gaussian", maxIter = 10, epsilon = 1e-5, verbose = FALSE, method.rq = "fn", method.nlrq = "L-BFGS-B"){

cl <- match.call()

delta <- 1e-6
taus <- sort(taus)
if(length(taus) != 2) stop("Provide only two quantile levels")
if(diff(taus) == 0) stop("Provide two distinct quantile levels")
p2 <- taus[1]
p1 <- taus[2]

mf.old <- mf <- model.frame(formula, data)
intercept <- attr(terms(mf), "intercept") == 1
x <- model.matrix(formula, mf)
xn <- if(intercept) x[,c(-1),drop=FALSE] else x
y <- model.extract(mf, "response")
n <- length(y)

if(any(y <= 0)) stop("The response must be strictly positive")

H1 <- quantile(y, probs = p1, names = FALSE)

if(length(start) > 1){
	# start with user-provided values
	stopifnot(is.numeric(start))
	stopifnot(length(start) == n)
	H2 <- start
} else {
	if(start == "rq"){
		# start with rq
		fit2 <- try(rq.fit(x = x, y = y, tau = p2, method = "fn"), silent = TRUE)
		if(inherits(fit2, "try-error")) stop("rq failed with error: ", fit2, "Change start estimate")
		bhat2 <- fit2$coef
		H2 <- x%*%bhat2
		if(any(H2 <= 0)) warning("Some predictions are not strictly positive")
	}

	if(start == "tsrq"){
		# start with tsrq
		fit2 <- try(tsrq(formula, data = data, tau = p2, tsf = tsf, symm = symm, dbounded = dbounded, method = "fn"), silent = TRUE)
		if(inherits(fit2, "try-error")) stop("tsrq failed. Change start estimate")
		H2 <- predict(fit2, type = "response")
	}
	
	if(start == "conquer"){
		# start with conquer (convolution-type smoothing)
		fit2 <- try(conquer(X = xn, Y = y, tau = p2, kernel = kernel), silent = TRUE)
		if(inherits(fit2, "try-error")) stop("conquer failed. with error: ", fit2, "Change start estimate")
		bhat2 <- fit2$coeff
		H2 <- x%*%bhat2
		if(any(H2 <= 0)) warning("Some predictions are not strictly positive")
	}

	if(start == "llqr"){
		# start with llqr
		fit2 <- try(llqr(x = xn, y = y, tau = p2), silent = TRUE)
		if(inherits(fit2, "try-error")) stop("llqr failed. Change start estimate")
		H2 <- fit2$ll_est
	}
}

H2 <- as.numeric(H2)

tt <- colnames(x)
tt <- gsub(" ", "", tt, fixed = TRUE)
colnames(x) <- tt
tt <- tt[tt != "(Intercept)"]
tt <- if(intercept) c("b1", paste(tt, paste0("b", 2:ncol(x)), sep = "*")) else paste(tt, paste0("b", 1:ncol(x)), sep = "*")
lp <- paste(tt, collapse = " + ")

if(linearize){
	ff1 <- ff2 <- update.formula(formula, "z ~ .")
} else {
	ff1 <- as.formula(paste("z ~ 1 + exp(", lp, ")"))
	ff2 <- as.formula(paste("z ~ 1/(1 + exp(", lp, "))"))
	mf <- data.frame(mf, x)
}

iter <- 0
h1 <- h2 <- NULL

while(iter < maxIter){

if(verbose) cat("Iteration", iter + 1, "\n")

if(linearize){
	z <- y/H2 - 1
	z[z <= 0] <- delta
	z1 <- mf$z <- log(z)
	if(method.rq == "conquer"){
		fit1 <- conquer(X = xn, Y = mf$z, tau = p1, kernel = kernel)
		bhat1 <- fit1$coeff
		h1 <- fit1$bandwidth
	} else {
		fit1 <- rq(ff1, tau = p1, data = mf, method = method.rq)
		bhat1 <- coef(fit1)
	}
	zhat <- if(method.rq %in% c("pfn", "conquer")) x%*%bhat1 else predict(fit1)
	H1 <- (1 + exp(zhat))*H2
} else {
	mf$z <- y/H2
	fit1 <- nlrq(ff1, tau = p1, data = mf, start = as.list(bhat1), method = method.nlrq)
	bhat1 <- coef(fit1)
	H1 <- predict(fit1)*H2
}
if(verbose) cat("gamma1", bhat1, "\n")

if(linearize){
	z <- H1/y - 1
	z[z <= 0] <- delta
	z2 <- mf$z <- log(z)
	if(method.rq == "conquer"){
		fit2 <- conquer(X = xn, Y = mf$z, tau = 1 - p2, kernel = kernel)
		bhat2 <- fit2$coeff
		h2 <- fit2$bandwidth
	} else {
		fit2 <- rq(ff2, tau = 1 - p2, data = mf, method = method.rq)
		bhat2 <- coef(fit2)
	}
	zhat <- if(method.rq %in% c("pfn", "conquer")) x%*%bhat2 else predict(fit2)
	H2 <- H1/(1 + exp(zhat))
} else {
	mf$z <- y/H1
	fit2 <- nlrq(ff2, tau = p2, data = mf, start = as.list(bhat2), method = method.nlrq)
	bhat2 <- coef(fit2)
	H2 <- predict(fit2)*H1
}

if(verbose) cat("gamma2", bhat2, "\n")

if(max(abs(bhat1 - bhat2)) < epsilon) {
		cat("Algorithm converged", "\n")
		break
	}
iter <- iter + 1
}

if(iter == maxIter) warning("Algorithm reached maximum number of iterations")

if(linearize){
	if(method.rq == "conquer"){
		omega1 <- vcov_conquer(res = fit1$residual, x = x, tau = p1, h = h1, intercept = intercept)
		omega2 <- vcov_conquer(res = fit2$residual, x = x, tau = 1 - p2, h = h2, intercept = intercept)
	} else {
		omega1 <- summary(fit1, se = "nid", covariance = TRUE)$cov
		omega2 <- summary(fit2, se = "nid", covariance = TRUE)$cov
	}
} else {
	omega1 <- summary(fit1)$cov
	omega2 <- summary(fit2)$cov
}

names(bhat2) <- colnames(x)
ans <- list(call = cl, formula = formula, coef = bhat2, p1 = p1, p2 = p2, omega1 = omega1, omega2 = omega2, ff1 = ff1, ff2 = ff2, data = mf.old, x = x, y = y, H1 = as.numeric(H1), H2 = as.numeric(H2), method.rq = method.rq, intercept = intercept, taus = taus, start = start, tsf = tsf, symm = symm, dbounded = dbounded, linearize = linearize, kernel = kernel, maxIter = maxIter, epsilon = epsilon, verbose = verbose, method.rq = method.rq, method.nlrq = method.nlrq)
attr(ans, "linearize") <- linearize
class(ans) <- "qrr"
return(ans)
}

coef.qrr <- coefficients.qrr <- function(object, ...){

object$coef

}

predict.qrr <- function(object, newdata, na.action = na.pass, type = "response", ...){

tau <- object$tau
nq <- length(tau)
bhat <- coef(object)

if(missing(newdata)) {x <- object$x} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
}

linpred <- x %*% bhat

if(type == "link"){
	return(linpred)
}

if(type == "response"){
	return(1 + exp(linpred))
}


}

qprime <- function(bhat, object) {
	x <- object$x
	y <- object$y
	p1 <- object$p1
	p2 <- object$p2
	H1 <- object$H1
	H2 <- object$H2
	
	n <- length(y)
	p <- ncol(x)
	h <- max(((log(n)+p)/n)^0.4,0.05)
	delta <- 1e-6
	
	for(j in 1:100) {
		H2 <- H1/(1 + exp(x%*%bhat))
		H1 <- H2*(1 + exp(x%*%bhat))
	}

	# FLAG <- TRUE
	# iter <- 0
	# while(FLAG){
		# H2.n <- H1/(1 + exp(x%*%bhat))
		# H1.n <- H2*(1 + exp(x%*%bhat))
		# eps <- max(abs(c(H1.n, H2.n) - c(H1, H2)))
		# H1 <- H1.n
		# H2 <- H2.n
		# if(eps < 1e-6){
			# FLAG <- FALSE
		# } else {
			# iter <- iter + 1
			# if(iter > 99) FLAG <- FALSE
		# }
	# }

	z <- y/H2-1
	z[z<=0] <- delta
	res <- log(z)-x%*%bhat 
	der <- pnorm(-res/h)-p1
	ans <- t(x)%*%der/n

	z <- H1/y-1
	z[z<=0] <- delta
	res <- log(z)-x%*%bhat 
	der <- pnorm(-res/h)-p2
	ans <- ans + t(x)%*%der/n
	return(ans)
}

vcov.qrr <- function(object, method = "approximate", R = 200, update = TRUE, ...){

ans <- NULL

if(method == "approximate"){
	stopifnot(attr(object, "linearize"))

	p1 <- object$p1
	p2 <- object$p2
	y <- object$y
	n <- length(y)
	H1 <- object$H1
	H2 <- object$H2
	R <- H1/H2

	ja <- jacobian(qprime, x = object$coef, object = object)
	ans <- solve(ja)/n
}

if(method == "boot"){
	if(update){
		f <- function(data, inds, object) update(object, data = data[inds,])$coef
		B <- boot::boot(object$data, f, R = R, stype = "i", parallel = "no", object = object)$t
	} else {
		f <- function(data, inds, object) {
			object$data <- data[inds,]
			do.call(qrr, object[match(names(object$call)[-c(1)], names(object))])$coef
		}
		B <- boot::boot(object$data, f, R = R, stype = "i", parallel = "no", object = object)$t
	}
	ans <- var(B)
}

return(ans)

}

summary.qrr <- function(object, se = "approximate", R = 200, update = TRUE, ...){

V <- vcov(object, method = se, R = R, update = update)
se <- sqrt(diag(V))
n <- length(object$y)
bhat <- coef(object)
p <- length(bhat)

object$tTable <- data.frame(bhat, se, bhat/se, 2*(1 - pt(abs(bhat/se), df = n - p)))
names(object$tTable) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")

class(object) <- "summary.qrr"

return(object)
}


print.qrr <- function(x, ...){

if (!is.null(cl <- x$cl)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}

bhat <- coef(x)
cat("\nQuantile ratio regression", paste(x$p1,x$p2,sep=":"), "\n")
cat("\nCoefficients linear predictor:\n")
print(bhat, ...)
nobs <- length(x$y)
p <- ncol(x$x)
rdf <- nobs - p
cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
}

print.summary.qrr <- function(x, ...){

if (!is.null(cl <- x$cl)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}

cat("\nQuantile ratio regression", paste(x$p1,x$p2,sep=":"), "\n")
cat("\nCoefficients linear predictor:\n")
printCoefmat(x$tTable, ...)
nobs <- length(x$y)
p <- ncol(x$x)
rdf <- nobs - p
cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")

}

vcov_conquer <- function(res, x, tau, h, intercept){
	n <- length(res)
	xn <- if(intercept) x[,c(-1),drop=FALSE] else x
	Wh <- as.vector(pnorm(-res/h) - tau)^2
	Stau <- crossprod(sweep(x, 1, Wh, "*"), x)/n
	Dh <- crossprod(sweep(x, 1, dnorm(res/h), "*"), x)/(n*h)
	Dhinv <- solve(Dh)
	V <- Dhinv %*% Stau %*% Dhinv/n
	return(V)
}
