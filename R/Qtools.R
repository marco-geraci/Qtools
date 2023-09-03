###            Qools: Utilities for quantiles
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed without any warranty,
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

".onAttach" <- function(lib, pkg) {
    if(interactive() || getOption("verbose"))
	packageStartupMessage(sprintf("Package %s (%s) loaded. Type citation(\"%s\") on how to cite this package\n", pkg,
		packageDescription(pkg)$Version, pkg))
}

##################################################
### Generics
##################################################

# Marginal effects

maref <- function(object, namevec) UseMethod("maref")

# Sparsity

sparsity <- function(object, se = "nid", hs = TRUE) UseMethod("sparsity")

######################################################################
### Confidence intervals for unconditional quantiles
######################################################################

qexact <- function(x, probs = 0.5, level = 0.95){

if(any(probs < 0) | any(probs > 1)) stop("Quantile index out of range: probs must be between 0 and 1")
if(any(is.na(x))) warning("Missing values will be omitted")
x <- as.numeric(na.omit(x))
n <- length(x)
z <- sort(x)
C <- gtools::combinations(n = n, r = 2, v = 1:n, repeats.allowed = FALSE)
nq <- length(probs)
out <- array(NA, dim = c(nq, 4), dimnames = list(paste0(probs*100, "%"), c("quantile", "lower", "upper", "conf.level")))
for(k in 1:nq){
	xi <- quantile(x, type = 1, probs = probs[k], names = FALSE)
	gamma <- pbinom(C[,2] - 1, size = n, prob = probs[k]) - pbinom(C[,1] - 1, size = n, prob = probs[k])
	gammac <- gamma - level
	if(all(gammac < 0)) stop(paste0("Quantile p = ", probs[k], ". Maximum confidence level for these data is ",	round(max(gamma), 3)))
	sel <- gammac == min(gammac[gammac >= 0])
	# if multiple solutions, take the narrowest CI
	if(sum(sel) > 1){
		sel <- which(sel)[which.min(z[C[sel,2]] - z[C[sel,1]])]
	}
	conf.level <- gamma[sel]
	CI <- z[C[sel,]]
	out[k,] <- c(xi, CI, conf.level)
	}
return(out)
}

##################################################
### Khmaladze and other tests
##################################################

KhmaladzeFormat <- function(object, epsilon){

if(!inherits(object, "KhmaladzeTest")) stop("class(object) must be 'KhmaladzeTest'")
tt <- get("KhmaladzeTable")
if(!(epsilon %in% unique(tt$epsilon))) stop("'epsilon' must be in c(0.05,0.10,0.15,0.20,0.25,0.30)")

p <- length(object$THn)
ans <- matrix(NA, p + 1, 2)
colnames(ans) <- c("Value", "Pr")

sel <- tt[tt$p == p & tt$epsilon == epsilon,3:5]
alpha <- c(0.01,0.05,0.1)

ans[1,1] <- object$Tn
ans[1,2] <- min(c(1,alpha[object$Tn > sel]))
ans[2:(p+1),1] <- object$THn

sig <- c("significant at 1% level", "significant at 5% level", "significant at 10% level", "not significant at 10% level")
null <- if(object$nullH == "location") "location-shift hypothesis" else "location-scale-shift hypothesis"

sel <- tt[tt$p == 1 & tt$epsilon == epsilon,3:5]
if(p==1){val <- min(c(1,alpha[object$THn > sel]))} else
{val <- rep(0,p); for(i in 1:p) val[i] <- min(c(1,alpha[object$THn[i] > sel]))}
ans[2:(p+1),2] <- val

mm <- match(ans[1,2], c(alpha,1))
nn <- match(ans[2:(p+1),2], c(alpha,1))

cat("Khmaladze test for the", null, "\n")
cat("Joint test is", sig[mm], "\n")
cat("Test(s) for individual slopes:", "\n")
for(i in 1:p){
cat(names(object$THn)[i], sig[nn][i], "\n")
}
invisible(ans)

}

normalize <- function(x){

n <- nrow(x)
p <- ncol(x)

if(p < 2 | n < p) stop("Provide n x p matrix with n > p > 1")

xx <- x

for (i in 2:p) {
	H <- c(crossprod(x[, i], xx[, 1:(i - 1)]))/diag(crossprod(xx[,1:(i - 1)]))
	xx[, i] <- x[, i] - matrix(xx[, 1:(i - 1)], nrow = n) %*% 
		matrix(H, nrow = i - 1)
}
H <- 1/sqrt(diag(crossprod(xx)))
xx <- t(t(xx) * H)

return(xx)

}

rcTest <- function(object, alpha = 0.05, B = 100, seed = NULL){

taus <- object$tau
nq <- length(taus)

x <- if(is.null(object[['x']])){
		do.call(model.matrix, args = list(object = as.formula(object$formula), data = object$model))
	} else {object[['x']]}
n <- nrow(x)
p <- ncol(x)
x <- normalize(x)*sqrt(n)
Rmat <- as.matrix(residuals(object))
if(!is.null(seed)) set.seed(seed)


Tn <- Tc <- pval <- vector()
for(j in 1:nq){
	tau <- taus[j]
	psi <- (Rmat[,j] > 0) * tau + (Rmat[,j] <= 0) * (tau - 1)
	omega <- replicate(B, sample(c(tau, -tau, 1-tau, tau-1), size = n, replace = TRUE, prob = c((1-tau)/2, (1-tau)/2, tau/2, tau/2)))

	ans <- C_rcTest(x, psi, omega, n, p, B)

	Tstar <- apply(ans$outstar, 3, function(x, n) eigen(x/n)$values[1], n = n)
	Tc[j] <- quantile(Tstar, 1 - alpha)
	Tn[j] <- eigen(ans$out/n)$values[1]
	pval[j] <- 1 - ecdf(Tstar)(Tn[j])
}

names(Tn) <- names(Tc) <- names(pval) <- taus

val <- list(Tn = Tn, Tcrit = Tc, p.value = pval, tau = taus)
class(val) <- "rcTest"
attr(val, "seed") <- if(!is.null(seed)) seed else NA
return(val)

}

GOFTest <- function(object, type = "cusum", alpha = 0.05, B = 100, seed = NULL){


val <- list()

val[[1]] <- switch(type,
	cusum = rcTest(object = object, alpha = alpha, B = B, seed = seed))

attr(val, "type") <- type
class(val) <- "GOFTest"
return(val)

}

print.GOFTest <- function (x, digits = max(3, getOption("digits") - 3), ...){

nt <- length(x)
type <- attributes(x)$type
txt <- vector()
txt[1] <- "Goodness-of-fit test for quantile regression based on the cusum process"
tau <- x[[1]]$tau
nq <- length(tau)

if(type == "cusum"){
	x <- x[[1]]
	cat(txt[1], "\n")
	cat("A large test statistic (small p-value) is evidence of lack of fit", "\n")
	for (j in 1:nq) {
		cat(paste0("Quantile ", tau[j], ": "))
		cat(paste0("Test statistic = ", round(x$Tn[j], digits), "; p-value = ", round(x$p.value[j],digits)), "\n")
	}
}

}

KhmaladzeTable <- structure(list(p = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 
11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 1L, 2L, 3L, 
4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 
18L, 19L, 20L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 
12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 1L, 2L, 3L, 4L, 
5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 
19L, 20L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 
13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 1L, 2L, 3L, 4L, 5L, 6L, 
7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 
20L), epsilon = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 
0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 
0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3), alpha01 = c(2.721, 4.119, 
5.35, 6.548, 7.644, 8.736, 9.876, 10.79, 11.81, 12.91, 14.03, 
15, 15.93, 16.92, 17.93, 18.85, 19.68, 20.63, 21.59, 22.54, 2.64, 
4.034, 5.267, 6.34, 7.421, 8.559, 9.573, 10.53, 11.55, 12.54, 
13.58, 14.65, 15.59, 16.52, 17.53, 18.46, 19.24, 20.21, 21.06, 
22.02, 2.573, 3.908, 5.074, 6.148, 7.247, 8.355, 9.335, 10.35, 
11.22, 12.19, 13.27, 14.26, 15.22, 16.12, 17.01, 17.88, 18.78, 
19.7, 20.53, 21.42, 2.483, 3.742, 4.893, 6.023, 6.985, 8.147, 
9.094, 10.03, 10.9, 11.89, 12.85, 13.95, 14.86, 15.69, 16.55, 
17.41, 18.19, 19.05, 19.96, 20.81, 2.42, 3.633, 4.737, 5.818, 
6.791, 7.922, 8.856, 9.685, 10.61, 11.48, 12.48, 13.54, 14.34, 
15.26, 16, 16.81, 17.59, 18.49, 19.4, 20.14, 2.32, 3.529, 4.599, 
5.599, 6.577, 7.579, 8.542, 9.413, 10.27, 11.15, 12.06, 12.96, 
13.82, 14.64, 15.46, 16.25, 17.04, 17.85, 18.78, 19.48), alpha05 = c(2.14, 
3.393, 4.523, 5.56, 6.642, 7.624, 8.578, 9.552, 10.53, 11.46, 
12.41, 13.34, 14.32, 15.14, 16.11, 16.98, 17.9, 18.83, 19.72, 
20.58, 2.102, 3.287, 4.384, 5.43, 6.465, 7.412, 8.368, 9.287, 
10.26, 11.17, 12.1, 13, 13.9, 14.73, 15.67, 16.56, 17.44, 18.32, 
19.24, 20.11, 2.048, 3.199, 4.269, 5.284, 6.264, 7.197, 8.125, 
9.044, 9.963, 10.85, 11.77, 12.61, 13.48, 14.34, 15.24, 16.06, 
16.93, 17.8, 18.68, 19.52, 1.986, 3.1, 4.133, 5.091, 6.07, 6.985, 
7.887, 8.775, 9.672, 10.52, 11.35, 12.22, 13.09, 13.92, 14.77, 
15.58, 16.43, 17.3, 18.09, 18.95, 1.923, 3, 4.018, 4.948, 5.853, 
6.76, 7.611, 8.51, 9.346, 10.17, 10.99, 11.82, 12.66, 13.46, 
14.33, 15.09, 15.95, 16.78, 17.5, 18.3, 1.849, 2.904, 3.883, 
4.807, 5.654, 6.539, 7.357, 8.211, 9.007, 9.832, 10.62, 11.43, 
12.24, 13.03, 13.85, 14.61, 15.39, 16.14, 16.94, 17.74), alpha1 = c(1.872, 
3.011, 4.091, 5.104, 6.089, 7.047, 7.95, 8.89, 9.82, 10.72, 11.59, 
12.52, 13.37, 14.28, 15.19, 16.06, 16.97, 17.84, 18.73, 19.62, 
1.833, 2.946, 3.984, 4.971, 5.931, 6.852, 7.77, 8.662, 9.571, 
10.43, 11.29, 12.2, 13.03, 13.89, 14.76, 15.65, 16.53, 17.38, 
18.24, 19.11, 1.772, 2.866, 3.871, 4.838, 5.758, 6.673, 7.536, 
8.412, 9.303, 10.14, 10.98, 11.86, 12.69, 13.48, 14.36, 15.22, 
16.02, 16.86, 17.7, 18.52, 1.73, 2.781, 3.749, 4.684, 5.594, 
6.464, 7.299, 8.169, 9.018, 9.843, 10.66, 11.48, 12.31, 13.11, 
13.91, 14.74, 15.58, 16.37, 17.17, 17.97, 1.664, 2.693, 3.632, 
4.525, 5.406, 6.241, 7.064, 7.894, 8.737, 9.517, 10.28, 11.11, 
11.93, 12.67, 13.47, 14.26, 15.06, 15.83, 16.64, 17.38, 1.602, 
2.602, 3.529, 4.365, 5.217, 6.024, 6.832, 7.633, 8.4, 9.192, 
9.929, 10.74, 11.51, 12.28, 13.05, 13.78, 14.54, 15.3, 16.05, 
16.79)), .Names = c("p", "epsilon", "alpha01", "alpha05", "alpha1"
), class = "data.frame", row.names = c(NA, -120L))

