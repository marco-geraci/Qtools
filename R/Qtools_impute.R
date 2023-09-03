##################################################
### Multiple imputation
##################################################

mice.impute.rq <- function (y, ry, x, tsf = "none", symm = TRUE, dbounded = FALSE, lambda = NULL, x.r = NULL, par = NULL, conditional = TRUE, epsilon = 0.001, method.rq = "fn", ...) {

	if(!tsf %in% c("none", "mcjI", "bc", "ao")) stop(paste("Transformation"), tsf, "not available")
	if(is.null(x.r)) x.r <- range(y, na.rm = TRUE)
    isDbounded <- (tsf == "mcjI" && dbounded)
    isDbounded <- tsf == "ao" || isDbounded

    y.old <- y
	x <- cbind(1, as.matrix(x))

	m <- sum(ry)
	n <- sum(!ry)
    p <- ncol(x)
	sel <- sample(1:m, m, replace = TRUE)
	
    xobs <- x[ry, ]
	xobs <- xobs[sel,]
    xmis <- x[!ry, ]
	u <- round(runif(n, epsilon, 1 - epsilon) * 1000)
	u <- ifelse(u %in% c(1:4, 996:999), u/1000, (u - u%%5)/1000)
	taus <- unique(u)
	nt <- length(taus)
	if(tsf == "none") conditional <- TRUE
	
	if(conditional){
		if (is.null(lambda)) lambda <- 0
		if (tsf %in% c("mcjI", "bc", "ao")) {
			if (isDbounded) y <- map(y, x.r = x.r)
			z <- switch(tsf, mcjI = mcjI(y, lambda, symm, dbounded, 
            omega = 0.001), bc = bc(y, lambda), ao = ao(y, lambda, 
            symm, omega = 0.001))
		} else if(tsf == "none"){
        z <- y
		}
		yobs <- z[ry]
		yobs <- yobs[sel]
		fit <- matrix(NA, p, nt)
		for (j in 1:nt) {
			fit[, j] <- as.numeric(rq.fit(xobs, yobs, tau = taus[j], 
				method = method.rq)$coefficients)
		}
		ypred <- xmis %*% fit
		ypred <- diag(ypred[, match(u, taus), drop = FALSE])
		if(tsf %in% c("mcjI", "bc", "ao")) {
			val <- switch(tsf, mcjI = invmcjI(ypred, lambda, symm, 
				dbounded), bc = invbc(ypred, lambda), ao = invao(ypred, 
				lambda, symm))
			if (isDbounded) val <- invmap(val, x.r)
		} else {val <- ypred}
	} else {
		yobs <- y.old[ry]
		yobs <- yobs[sel]
		ypred <- matrix(NA, n, nt)
		fit <- nlrq1(yobs ~ xobs - 1, tsf = tsf, symm = symm, dbounded = dbounded, tau = taus)
		lambda <- fit$lambda
		lambda[is.na(lambda)] <- 0
		fit <- tsrq(yobs ~ xobs - 1, tsf = tsf, symm = symm, dbounded = dbounded, tau = taus, conditional = TRUE, lambda = lambda, method = method.rq)
		for(j in 1:nt) {
			ypred[,j] <- invmcjI(xmis%*%fit$coef[ , j, drop = FALSE], lambda = lambda[j], symm = symm, dbounded = dbounded)
		} # for
		val <- diag(ypred[, match(u, taus), drop = FALSE])
		if (isDbounded) val <- invmap(val, x.r)
	} # else

return(val)

}

# Impute using restricted quantiles

mice.impute.rrq <- function (y, ry, x, tsf = "none", symm = TRUE, dbounded = FALSE, lambda = NULL, epsilon = 0.001, method.rq = "fn", ...) 
{
    x <- cbind(1, as.matrix(x))
	y.old <- y
	isDbounded <- (tsf == "mcjI" && dbounded)
	isDbounded <- tsf == "ao" || isDbounded

	if(is.null(lambda))
		lambda <- 0

	if(tsf %in% c("mcjI","bc","ao")){
		if(isDbounded) y <- map(y)
		z <- switch(tsf,
			mcjI = mcjI(y, lambda, symm, dbounded, omega = 0.001),
			bc = bc(y, lambda),
			ao = ao(y, lambda, symm, omega = 0.001)
			)
	} else {z <- y}

	m <- sum(ry)
	n <- sum(!ry)
	p <- ncol(x)
	sel <- sample(1:m, m, replace = TRUE)

	xobs <- x[ry, ]
	xobs <- xobs[sel,]
	xmis <- x[!ry,]
	yobs <- z[ry]
	yobs <- yobs[sel]

	u <- round(runif(n, epsilon, 1 - epsilon)*1e3)
	u <- ifelse(u %in% c(1:4,996:999), u/1e3, (u - u %% 5)/1e3)
	taus <- unique(u)
	nt <- length(taus)

	fit <- matrix(NA, p, nt)
	for(j in 1:nt){
		fit[,j] <- as.numeric(rrq.fit(xobs, yobs, tau = taus[j], method = method.rq)$coef)
	}
	# n times nt matrix
	ypred <- xmis%*%fit
	# diagonal of n times n matrix
	ypred <- diag(ypred[, match(u, taus), drop = FALSE])
	
	if(tsf %in% c("mcjI","bc","ao")){
		val <- switch(tsf,
			mcjI = invmcjI(ypred, lambda, symm, dbounded),
			bc = invbc(ypred, lambda),
			ao = invao(ypred, lambda, symm));
		if(isDbounded) val <- invmap(val, range(y.old))
	} else {val <- ypred}

    return(val)
}
