########################################################################
###Non_linear location parameter update to gevr_fit{eva pkg} function###
########################################################################


### gevr_fit_nloc function
## This function fits the gevR model to the data.
## Implemented new nonlinear location parameter capability



## Function to help deal with design matrix
adjScale <- function(x) {
  truemeans <- as.numeric(colMeans(x))
  truevars <- as.numeric(apply(x, 2, sd))
  adjmeans <- ifelse(truevars == 0, 0, truemeans)
  adjvars <- ifelse(truevars == 0, truemeans, truevars)
  if(ncol(x) == 1)
    adjmeans <- 0
  x <- t((t(x) - adjmeans) / adjvars)
  out <- list(x, truemeans, truevars, adjmeans, adjvars)
  names(out) <- c("mat", "truemeans", "truevars", "adjmeans", "adjvars")
  out
}

## Returns expected inverse fisher information matrix for GEVr or Gumbel distribution.
gevrFisher <- function(data, theta, gumbel = FALSE) {
  data <- as.matrix(data)
  R <- ncol(data)
  N <- nrow(data)
  loc <- theta[1]
  scale <- theta[2]
  if(!gumbel) {
    shape <- theta[3]
    gr1 <- gamma(R + shape + 1) / gamma(R)
    gr2 <- gamma(R + 2*shape + 1) / gamma(R)
    A <- (((1+shape)^2)*gr2)/((scale^2)*(1+2*shape))
    B <- (-1/((scale^2)*shape*(1+2*shape)))*(((1+shape)^2)*gr2 - (1+2*shape)*gr1)
    C <- (1/(scale*(shape^2)*(1+2*shape)))*(((1+2*shape)*gr1)*(shape*digamma(R+shape+1) + (shape^2+shape+1)/(1+shape)) - ((1+shape)^2)*gr2)
    D <- (1/((scale^2)*(shape^2)*(1+2*shape)))*(R*(1+2*shape) - 2*(1+2*shape)*gr1 + ((1+shape)^2)*gr2)
    E <- (1/(scale*((-shape)^3)*(1+2*shape)))*(R*(-shape)*(1+2*shape)*digamma(R+1) + (1+2*shape)*gr1*(shape*digamma(R+shape+1) + (1+(1+shape)^2)/(1+shape)) - ((1+shape)^2)*gr2 - R*(1+2*shape))
    F <- (1/(((-shape)^4)*(1+2*shape)))*((2*(1+2*shape)*gr1)*((-shape)*digamma(R+shape+1) - (shape^2+shape+1)/(1+shape)) + ((1+shape)^2)*gr2 + (R*(1+2*shape))*(1 + 2*shape*digamma(R+1) + (shape^2)*(1 + trigamma(R+1) + ((digamma(R+1))^2))))
    info <- solve(matrix(c(A, B, C, B, D, E, C, E, F), nrow=3, ncol=3)) / N
  } else {
    Br <- R * digamma(R + 1)
    Cr <- R * (digamma(R + 1)^2 + trigamma(R + 1) + 1)
    info <- (scale^2 / (N * (R * Cr - Br^2))) * matrix(c(Cr, Br, Br, R), nrow=2, ncol=2)
  }
  info
}





gevr_fit_nloc <- function(data, method = c("mle", "mps", "pwm"), information = c("expected", "observed"), locvars = NULL, scalevars = NULL,
                       shapevars = NULL, locform = ~ 1, scaleform = ~ 1, shapeform = ~ 1, loclink = identity, scalelink = identity,
                       shapelink = identity, gumbel = FALSE, start = NULL, opt = "Nelder-Mead", maxit = 10000, ...) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  method <- match.arg(method)
  information <- match.arg(information)
  if((gumbel) & ((!is.null(shapevars)) | (shapeform != ~ 1)))
    stop("Cannot specify covariates in shape parameter for Gumbel distribution!")
  if(!is.null(locvars))
    if(nrow(locvars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(!is.null(scalevars))
    if(nrow(scalevars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(!is.null(shapevars))
    if(nrow(shapevars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(((locform != ~ 1) & is.null(locvars)) | ((scaleform != ~ 1) & is.null(scalevars)) | ((shapeform != ~ 1) & is.null(shapevars)))
    stop("Need to specify covariates!")
  if(R > 1 & (method == "mps" | method == "pwm"))
    stop("If R > 1, MLE must be used")
  if((!is.null(locvars) | !is.null(scalevars) | !is.null(shapevars)) & method == "pwm")
    stop("Probability weighted moments can only be fitted for stationary data")
  
  if(locform == ~ 1)
    locvars <- as.data.frame(rep(1, n))
  
  if(scaleform == ~ 1)
    scalevars <- as.data.frame(rep(1, n))
  
  if(shapeform == ~ 1)
    shapevars <- as.data.frame(rep(1, n))
  
  locvars.model <- model.matrix(locform, data = locvars)
  locnames <- colnames(locvars.model)
  loccheck <- adjScale(locvars.model)
  if((rankMatrix(locvars.model)[1] < ncol(locvars.model)) | (rankMatrix(locvars.model)[1] > nrow(locvars.model)))
    stop("Location design matrix is singular")
  locvars.model <- loccheck$mat
  loctrans1 <- loccheck$adjmeans
  loctrans2 <- loccheck$adjvars
  
  scalevars.model <- model.matrix(scaleform, data = scalevars)
  scalenames <- colnames(scalevars.model)
  scalecheck <- adjScale(scalevars.model)
  if((rankMatrix(scalevars.model)[1] < ncol(scalevars.model)) | (rankMatrix(scalevars.model)[1] > nrow(scalevars.model)))
    stop("Scale design matrix is singular")
  scalevars.model <- scalecheck$mat
  scaletrans1 <- scalecheck$adjmeans
  scaletrans2 <- scalecheck$adjvars
  
  shapevars.model <- model.matrix(shapeform, data = shapevars)
  shapenames <- colnames(shapevars.model)
  shapecheck <- adjScale(shapevars.model)
  if((rankMatrix(shapevars.model)[1] < ncol(shapevars.model)) | (rankMatrix(shapevars.model)[1] > nrow(shapevars.model)))
    stop("Shape design matrix is singular")
  shapevars.model <- shapecheck$mat
  shapetrans1 <- shapecheck$adjmeans
  shapetrans2 <- shapecheck$adjvars
  
  trans1 <- c(loctrans1, scaletrans1, shapetrans1)
  trans2 <- c(loctrans2, scaletrans2, shapetrans2)
  
  locvars.model.orig <- t((t(locvars.model) * loctrans2) + loctrans1)
  scalevars.model.orig <- t((t(scalevars.model) * scaletrans2) + scaletrans1)
  shapevars.model.orig <- t((t(shapevars.model) * shapetrans2) + shapetrans1)
  
  ## Probability Weighted Moments
  ## Also use this as the intial estimates for other methods
  x <- sort(as.vector(data[, 1]))
  solve_shape <- function(sh, mom0, mom1, mom2) {
    (3^sh - 1) / (2^sh - 1) - (3 * mom2 - mom0) / (2 * mom1 - mom0)
  }
  moments <- rep(0, 3)
  for(j in 1:n) {
    moments[1] <- moments[1] + (x[j] / n)
    moments[2] <- moments[2] + ((j - 1) / (n - 1)) * (x[j] / n)
    moments[3] <- moments[3] + (((j - 1) * (j - 2)) / ((n - 1) * (n - 2))) * (x[j] / n)
  }
  mom0 <- moments[1]
  mom1 <- moments[2]
  mom2 <- moments[3]
  
  if(!gumbel) {
    shape0 <- uniroot(solve_shape, interval = c(-5, +5), mom0 = mom0, mom1 = mom1, mom2 = mom2)$root
    scale0 <- ((2 * mom1 - mom0) * shape0) / (gamma(1 - shape0) * (2^shape0 - 1))
    loc0 <- mom0 + scale0 * (1 - gamma(1 - shape0)) / shape0
    theta0 <- c(loc0, scale0, shape0)
  } else {
    scale0 <- (mom0 - 2 * mom1) / - log(2)
    loc0 <- mom0 + scale0 * digamma(1)
    theta0 <- c(loc0, scale0)
  }
  
  if(is.null(start)) {
    locinit <- c(loc0, rep(0, ncol(locvars.model) - 1))
    scaleinit <- c(scale0, rep(0, ncol(scalevars.model) - 1))
    if(!gumbel) {
      shapeinit <- c(shape0, rep(0, ncol(shapevars.model) - 1))
      init <- c(locinit, scaleinit, shapeinit)
    } else {
      init <- c(locinit, scaleinit)
    }
  } else {
    init <- start
    locinit <- init[1:ncol(locvars.model)]
    scaleinit <- init[(ncol(locvars.model) + 1):(ncol(locvars.model) + ncol(scalevars.model))]
    shapeinit <- init[(ncol(locvars.model) + ncol(scalevars.model) + 1):length(init)]
  }
  
  parnum <- c(ncol(locvars.model), ncol(scalevars.model), ncol(shapevars.model))
  if(gumbel) parnum[3] <- 0
  
  if(method == "pwm") {
    names(theta0) <- c("Location", "Scale", "Shape")[1:length(theta0)]
    out <- list(n = n, data = data, par.ests = theta0, par.ses = NA, varcov = NA,
                converged = NA, nllh.final = NA, R = R,
                stationary = TRUE, parnum = parnum,
                par.sum = theta0, gumbel = gumbel,
                covars = list(locvars.model.orig, scalevars.model.orig, shapevars.model.orig),
                covars.orig = list(locvars, scalevars, shapevars),
                links = list(loclink, scalelink, shapelink),
                forms = list(locform, scaleform, shapeform),
                method = method, information = information)
  }
  
  negloglik <- function(vars, locvars1, scalevars1, shapevars1, x) {
    loc <- vars[1:length(locinit)]
    if(any(loc[2:length(loc)] < 0)) return(NA) #added in
    scale <- vars[(length(locinit) + 1):(length(locinit) + length(scaleinit))]
    if(!gumbel) shape <- vars[(length(locinit) + length(scaleinit) + 1):length(vars)] else shape <- 0
    locmat <- t(loc * t(locvars1))
    scalemat <- t(scale * t(scalevars1))
    shapemat <- t(shape * t(shapevars1))
    locvec <- loclink(rowSums(locmat))
    scalevec <- scalelink(rowSums(scalemat))
    shapevec <- shapelink(rowSums(shapemat))
    w <- matrix(((x - locvec) / scalevec), ncol = R)
    z <- pmax(matrix(w * shapevec, ncol = R), -1)
    cond1 <- any(scalevec <= 0)
    cond2 <- min(1 + w * shapevec) <= 0
    log.density <- ifelse(shapevec == 0, rowSums(-log(pmax(scalevec, 0)) - w) - exp(-w[,R]),
                          rowSums(-log(pmax(scalevec, 0)) - ((1/shapevec) + 1) * log1p(z)) -
                            exp((-1/shapevec) * log1p(z[,R])))
    log.density[(is.nan(log.density) | is.infinite(log.density))] <- 0
    if(cond1 | cond2) {
      abs(sum(log.density)) + 1e6
    } else {
      - sum(log.density)
    }
  }
  
  mpsobj <- function(vars, locvars1, scalevars1, shapevars1, x) {
    loc <- vars[1:length(locinit)]
    scale <- vars[(length(locinit) + 1):(length(locinit) + length(scaleinit))]
    if(!gumbel) shape <- vars[(length(locinit) + length(scaleinit) + 1):length(vars)] else shape <- 0
    locmat <- t(loc * t(locvars1))
    scalemat <- t(scale * t(scalevars1))
    shapemat <- t(shape * t(shapevars1))
    locvec <- loclink(rowSums(locmat))
    scalevec <- scalelink(rowSums(scalemat))
    shapevec <- shapelink(rowSums(shapemat))
    w <- as.vector((x - locvec) / scalevec)
    z <- pmax(w * shapevec, -1)
    cond1 <- any(scalevec <= 0)
    cond2 <- min(1 + w * shapevec) <= 0
    cdf <- ifelse(shapevec == 0, exp(-exp(-w)), exp(-exp((-1/shapevec)*log1p(z))))
    cdf[(is.nan(cdf) | is.infinite(cdf))] <- 0
    cdf <- c(0, cdf, 1)
    D <- diff(cdf)
    cond3 <- any(D < 0)
    ## Check if any differences are zero due to rounding and adjust
    D <- ifelse(D <= 0, .Machine$double.eps, D)
    if(cond1 | cond2 | cond3) {
      abs(sum(log(D))) + 1e6
    } else {
      - sum(log(D))
    }
  }
  
  
  if(method != "pwm") {
    if(method == "mle") {
      fit <- optim(init, negloglik, hessian = FALSE, method = opt, control = list(maxit = maxit, ...),
                   locvars1 = locvars.model, scalevars1 = scalevars.model, shapevars1 = shapevars.model,
                   x = data)
    } else {
      data.order <- order(data)
      data.sort <- as.matrix(sort(data))
      locvars.model.sort <- apply(locvars.model, 2, function(x) x[data.order])
      scalevars.model.sort <- apply(scalevars.model, 2, function(x) x[data.order])
      shapevars.model.sort <- apply(shapevars.model, 2, function(x) x[data.order])
      locvars.model.orig.sort <- apply(locvars.model.orig, 2, function(x) x[data.order])
      scalevars.model.orig.sort <- apply(scalevars.model.orig, 2, function(x) x[data.order])
      shapevars.model.orig.sort <- apply(shapevars.model.orig, 2, function(x) x[data.order])
      fit <- optim(init, mpsobj, hessian = FALSE, method = opt, control = list(maxit = maxit, ...),
                   locvars1 = locvars.model.sort, scalevars1 = scalevars.model.sort,
                   shapevars1 = shapevars.model.sort, x = data.sort)
    }
    
    if(fit$convergence)
      warning("optimization may not have succeeded")
    
    loc.ests <- fit$par[1:length(locinit)] / loctrans2
    scale.ests <- fit$par[(length(locinit) + 1):(length(locinit) + length(scaleinit))] / scaletrans2
    if(!gumbel)
      shape.ests <- fit$par[(length(locinit) + length(scaleinit) + 1):length(fit$par)] / shapetrans2
    
    loc.ests <- ifelse(loccheck$truevars == 0, loc.ests - sum(loc.ests * loctrans1), loc.ests)
    scale.ests <- ifelse(scalecheck$truevars == 0, scale.ests - sum(scale.ests * scaletrans1), scale.ests)
    if(!gumbel)
      shape.ests <- ifelse(shapecheck$truevars == 0, shape.ests - sum(shape.ests * shapetrans1), shape.ests)
    
    if(!gumbel) par.ests <- c(loc.ests, scale.ests, shape.ests) else par.ests <- c(loc.ests, scale.ests)
    
    if((information == "observed") | (locform != ~ 1) | (scaleform != ~ 1) | (shapeform != ~ 1)) {
      if(method == "mle") {
        varcov <- NA # solve(optimHess(par.ests, negloglik, locvars1 = locvars.model.orig,
        #                          scalevars1 = scalevars.model.orig, shapevars1 = shapevars.model.orig,
        #                          x = data))
      } else {
        varcov <- solve(optimHess(par.ests, mpsobj, locvars1 = locvars.model.orig.sort,
                                  scalevars1 = scalevars.model.orig.sort, shapevars1 = shapevars.model.orig.sort,
                                  x = data.sort))
      }
    } else {
      varcov <- gevrFisher(data, par.ests, gumbel)
    }
    par.ses <- diag(length(par.ests))
    
    if(!gumbel) {
      names(par.ests) <- c(paste('Location', colnames(locvars.model.orig), sep = ' '),
                           paste('Scale', colnames(scalevars.model.orig), sep = ' '),
                           paste('Shape', colnames(shapevars.model.orig), sep = ' '))
    } else {
      names(par.ests) <- c(paste('Location', colnames(locvars.model.orig), sep = ' '),
                           paste('Scale', colnames(scalevars.model.orig), sep = ' '))
    }
    
    names(par.ses) <- names(par.ests)
    
    par.sum <- data.frame(par.ests, par.ses, par.ests / par.ses, 2 * pnorm(abs(par.ests / par.ses), lower.tail = FALSE))
    colnames(par.sum) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    par.sum$codes <- ifelse(par.sum[, 4] < 0.001, '***',
                            ifelse(par.sum[, 4] < 0.01, '**',
                                   ifelse(par.sum[, 4] < 0.05, '*',
                                          ifelse(par.sum[, 4] < 0.1, '.', ' '))))
    colnames(par.sum) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "")
    
    if(method == "mle") {
      out <- list(n = n, data = data, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                  converged = fit$convergence, nllh.final = fit$value, R = R,
                  stationary = ((locform == ~ 1) & (scaleform == ~ 1) & (shapeform == ~ 1)),
                  parnum = parnum, par.sum = par.sum, gumbel = gumbel,
                  covars = list(locvars.model.orig, scalevars.model.orig, shapevars.model.orig),
                  covars.orig = list(locvars, scalevars, shapevars),
                  links = list(loclink, scalelink, shapelink),
                  forms = list(locform, scaleform, shapeform),
                  method = method, information = information)
    } else {
      out <- list(n = n, data = as.matrix(data), par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                  converged = fit$convergence, moran = fit$value, R = R,
                  stationary = ((locform == ~ 1) & (scaleform == ~ 1) & (shapeform == ~ 1)),
                  parnum = parnum, par.sum = par.sum, gumbel = gumbel,
                  covars = list(locvars.model.orig, scalevars.model.orig, shapevars.model.orig),
                  covars.orig = list(locvars, scalevars, shapevars),
                  links = list(loclink, scalelink, shapelink),
                  forms = list(locform, scaleform, shapeform),
                  method = method, information = information)
    }
    
  }
  
  class(out) <- "gevr_fit_nloc"
  out
}


