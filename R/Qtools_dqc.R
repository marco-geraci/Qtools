##################################################
# Directional quantile classification
##################################################

dqcControl <- function(tau.range = c(0.001, 0.999), nt = 10, ndir = 50, seed = NULL){

list(tau.range = tau.range, nt = nt, ndir = ndir, seed = seed)

}

dqc <- function(formula, data, df.test, subset, weights, na.action, control = list(), fit = TRUE){

cl <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- quote(stats::model.frame)
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
intercept <- attr(terms.formula(formula, data = mf), "intercept") == 1

# train dataset
y <- model.response(mf)
x <- model.matrix(mt, mf)
w <- as.vector(model.weights(mf))
if (!is.null(w) && !is.numeric(w)) 
	stop("'weights' must be a numeric vector")
if(intercept){
	x <- x[,-c(1),drop = FALSE]
}

# test dataset
if(!is.null(df.test)){
	z <- model.matrix(formula(mt)[-2], df.test)
	if(intercept){
		z <- z[,-c(1),drop = FALSE]
	}
} else {
	z <- NULL
}

if(is.null(names(control)))
	control <- dqcControl()
else {
	control_default <- dqcControl()
	control_names <- intersect(names(control), names(control_default))
	control_default[control_names] <- control[control_names]
	control <- control_default
}

FIT_ARGS <- list(x = x, z = z, y = y, control = control)

if(!fit) return(FIT_ARGS)

res <- do.call(dqc.fit, FIT_ARGS)

res$call <- cl
res$nx <- nrow(x)
res$nz <- nrow(z)
res$p <- ncol(x)
res$control <- control
res$terms <- mt
res$term.labels <- colnames(x)
class(res) <- "dqc"
return(res)

}

dqc.fit <- function(x, z, y, control){

.checkfn <- function(x, p) x*(p - (x < 0))

if(!is.null(control$seed)) set.seed(control$seed)

y <- as.factor(y)
nx <- nrow(x)
nz <- nrow(z)
p <- ncol(x)

ndir <- control$ndir
nt <- control$nt
groups <- sort(levels(y))
ng <- length(groups)

if(is.null(row.names(x))) row.names(x) <- 1:nx
if(is.null(row.names(z))) row.names(z) <- 1:nz

# generate grid of taus
tau.range <- control$tau.range
if(any(is.na(tau.range))){
	stop("tau.range must not have NAs")
}
if(length(tau.range) == 1){
	taus <- tau.range
	nt <- 1
}
if(length(tau.range) == 2){
	tau.range <- sort(tau.range)
	taus <- seq(tau.range[1], tau.range[2], length = nt)
}
if(!length(tau.range) %in% c(1,2)){
	stop("I cannot understand 'tau.range'. It must be of length 1 or 2.")
}
if(any(taus <= 0) | any(taus >= 1)){
	stop("taus must be strictly in the unit interval (0,1)")
}

# order data
ord <- order(y)
x <- x[ord,]
y <- y[ord]
idx <- row.names(x)
idz <- row.names(z)

# marginal quantiles
csi <- lapply(split(data.frame(x), y), function(z, taus) apply(z, 2, quantile, probs = taus), taus = taus)
# element-wise signs for all combinations
c_groups <- gtools::combinations(n = ng, r = 2, v = groups, repeats.allowed = FALSE)
nc <- nrow(c_groups)
dir.sgn <- apply(c_groups, 1, function(i, x) sign(x[i][[1]] - x[i][[2]]), x = csi) # (nt x p) x nc

xu <- array(NA, dim = c(nx, ndir, nt), dimnames = list(obs = idx, dir = 1:ndir, tau = taus))
zu <- array(NA, dim = c(nz, ndir, nt), dimnames = list(obs = idz, dir = 1:ndir, tau = taus))

for (j in 1:nt) {
	# generate grid of directions uniformly over the p-dimensional unit sphere
	#theta <- matrix(runif(ndir * (p-1), 0, 2*pi), nrow = ndir, ncol = (p-1))
	#u <- mvmesh::Polar2Rectangular(r = rep(1, ndir), theta = theta) # ndir x p
	sgn.sel <- dir.sgn[seq(j, nt*p, by = nt),sample(1:nc, 1)]
	out <- C_projfun(x, z, sgn.sel, nx, nz, p, ndir)
	xu[,,j] <- out$xu
	zu[,,j] <- out$zu
}

B <- ndir*nt
xu <- matrix(as.numeric(xu), nrow = nx) # nx x B
zu <- matrix(as.numeric(zu), nrow = nz) # nx x B

ns <- as.integer(table(y))
minn <- c(0, cumsum(ns[-ng]))
maxn <- cumsum(ns)

Phi <- C_phifun(xu, zu, nx, nz, B, ndir, ng, taus, minn, maxn)

w <- colSums(Phi$out)
w <- -w/sqrt(sum(w^2))

dist.z <- matrix(0, nz, ng)
for (i in 1:ng){
	ss <- seq(i, ng*B, by = ng)
	dist.z[,i] <- Phi$Phi_z[, ss, drop = FALSE]%*%w
}
index <- apply(dist.z, 1, which.min)

ans <- data.frame(obs = idz, groups = factor(groups[index], levels = groups, labels = groups), value = apply(dist.z, 1, min))
list(ans = ans, groups = groups)

}

print.dqc <- function(x, ...){

z <- table(x$ans$groups)/x$nz*100

cat("Directional quantile classification", "\n")
if (!is.null(cl <- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	cat("Classification proportions (%) by class:\n")
	print(z)
	cat("\n")
}


}
