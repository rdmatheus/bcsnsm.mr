mr_mle <- function(y, X = NULL, Z = NULL, association = nonassociative(),
                   mu.link = "log", sigma.link = "log",
                   copula = "gaussian", eta = NULL, margin = "bcno",
                   control = control_fit(...), ...)
{

  ### Specific optim parameters
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian
  inits <- control$inits

  lambda_id <- !grepl("Log-", as.bcs(margin)$name)
  nu_id <- as.bcs(margin)$extrap

  control$method <- control$hessian <- control$inits <- NULL

  ### Data setting -------------------------------------------------------------
  n <- length(y)

  if (is.null(X)) X <- matrix(1, nrow = n)
  if (is.null(Z)) Z <- matrix(1, nrow = n)

  k1 <- ncol(X)
  k2 <- ncol(Z)

  ### Copula -------------------------------------------------------------------
  copula <- match.arg(copula, c("gaussian", "t", "slash", "hyp"))
  mcopula <- make_copula(copula, eta / (1 - eta))
  
  qPSI <- mcopula$qPSI
  dgf <- mcopula$dgf
  
  
  # Parameter indexation -------------------------------------------------------
  par_id <- list(beta = 1 : k1,
                   kappa = 1 : k2 +  k1,
                   lambda = if(lambda_id) 1 + k1 + k2 else NULL,
                   nu = if(any(nu_id)) 1:sum(as.numeric(nu_id)) + k1 + k2 + 1 * as.numeric(lambda_id) else NULL,
                   gamma = if(association$npar > 0){
                     1:association$npar + sum(as.numeric(nu_id)) + k1 + k2 + 1 * as.numeric(lambda_id)
                   } else {
                     NULL
                   })

  ## Initial values ------------------------------------------------------------
  gamma_inits <- association$start(y)

  links <- paste0(mu.link, "_", sigma.link)

  if (margin %in% c("bcpe", "bchp", "bcloi", "lpe", "lhp", "lloi")) {
    faux <- gamlss.dist::BCPE
    faux2 <- function(z, nu) gamlss.dist::pPE(z, 0, 1, nu)
  } else if (margin %in% c("bcno", "lno")) {
    faux <- gamlss.dist::BCCG
    faux2 <- function(z, nu) stats::pnorm(z)
  } else {
    faux <- gamlss.dist::BCT
    faux2 <- function(z, nu) stats::pt(z, nu)
  }


  gamlss_fit <- switch (links,
                        log_log = {
                          suppressWarnings(gamlss::gamlss(y ~ X + 0, ~ Z + 0,
                                                          family = faux(mu.link = "log", sigma.link = "log",
                                                                        nu.link = "identity"),
                                                          trace = FALSE))#, tau.fix = TRUE, tau.start = 4))
                        },

                        log_identity = {
                          suppressWarnings(gamlss::gamlss(y ~ X + 0, ~ Z + 0,
                                                          family = faux(mu.link = "log",sigma.link = "identity",
                                                                        nu.link = "identity"),
                                                          trace = FALSE))#, tau.fix = TRUE, tau.start = 4))
                        },

                        identity_log = {
                          suppressWarnings(gamlss::gamlss(y ~ X + 0, ~ Z + 0,
                                                          family = faux(mu.link = "identity", sigma.link = "log",
                                                                        nu.link = "identity"),
                                                          trace = FALSE))#, tau.fix = TRUE, tau.start = 4))
                        },

                        identity_identity = {
                          suppressWarnings(gamlss::gamlss(y ~ X + 0, ~ Z + 0,
                                                          family = faux(mu.link = "identity", sigma.link = "identity",
                                                                        nu.link = "identity"),
                                                          trace = FALSE))#, tau.fix = TRUE, tau.start = 4))
                        })

  inits <- c(stats::coef(gamlss_fit, "mu"),
             stats::coef(gamlss_fit, "sigma"),
             if(any(lambda_id)) stats::coef(gamlss_fit, "nu") else NULL,
             if(any(nu_id)) min(gamlss_fit$tau.fv[1], 50) else NULL,
             gamma_inits)

  ## Log-likelihood ------------------------------------------------------------
  EPS <- .Machine$double.eps^(1/1.6)

  ll <- function(theta){

    ### Parameter setting ------------------------------------------------------
    beta <- theta[par_id$beta]
    kappa <- theta[par_id$kappa]
    gamma <- if (association$npar > 0) theta[par_id$gamma] else NULL

    # Marginal parameters
    mu <- c(stats::make.link(mu.link)$linkinv(X%*%beta))
    sigma <- c(stats::make.link(sigma.link)$linkinv(Z%*%kappa))
    lambda <- if (lambda_id) theta[par_id$lambda] else NA
    nu <- if (nu_id) theta[par_id$nu] else NA

    ### Association matrix
    Gamma <- association$Gamma(gamma, n)
    if (!is.null(Gamma)) Gamma <- as.matrix(Gamma)
    
    dec <- tryCatch(Rfast::cholesky(Gamma), error = function(e) e)
    if (inherits(dec, "error")) dec <- NULL

    ### Out
    if (any(!is.finite(mu)) | any(!is.finite(sigma)) | any(!is.finite(lambda[lambda_id])) |
        any(mu < 0) | any(sigma < 0) | any(nu[nu_id] < 0) | is.null(dec)) {

      -Inf

    }else {

      epsilon <- matrix(qPSI(pmin(pmax(get(paste0("p", margin))(q = y,
                                                         mu = mu,
                                                         sigma = sigma,
                                                         lambda = lambda,
                                                         nu = nu), EPS), 1 - EPS)),
                        nrow = 1)

      log_f <- get(paste0("d", margin))(x = y,
                                        mu = mu,
                                        sigma = sigma,
                                        lambda = lambda,
                                        nu = nu, log = TRUE)


      tmp <- backsolve(dec, t(epsilon), transpose = TRUE)
      rss <- colSums(tmp^2)
      - sum(log(diag(dec))) + dgf(rss, n, log = TRUE) -
        sum(dgf(epsilon^2, 1L, log = TRUE)) + sum(log_f)

    }

  }


  ## Estimates -----------------------------------------------------------------
  opt <- stats::optim(par = inits,
                      fn = ll,
                      method = method,
                      control = control,
                      hessian = hessian)

  opt$par_id <- par_id
  opt$start <- inits

  # Convergence status
  if (opt$convergence > 0)
    warning(cat("optimization failed to converge\n"))

  # Coefficients and marginal parameters
  beta <- opt$par[par_id$beta]
  names(beta) <- colnames(X)
  kappa <- opt$par[par_id$kappa]
  names(kappa) <- colnames(Z)

  mu <- c(stats::make.link(mu.link)$linkinv(X%*%beta))
  sigma <- c(stats::make.link(sigma.link)$linkinv(Z%*%kappa))
  lambda <- if (lambda_id) opt$par[par_id$lambda] else NULL
  nu <- if (nu_id) opt$par[par_id$nu] else NULL

  # Association matrix parameters
  gamma <- if (association$npar) opt$par[par_id$gamma] else NULL

  # Assymptotic covariance estimates
  if (hessian){

    vcov <- try(chol2inv(Rfast::cholesky(-opt$hessian)), silent = TRUE)
    vcov <- if (unique(grepl("Error", vcov))) {
      matrix(NA, nrow = length(opt$par), ncol = length(opt$par)) 
    } else {
      vcov
    }

    colnames(vcov) <- rownames(vcov) <- c(colnames(X), colnames(Z),
                                          if (lambda_id) "lambda" else NULL,
                                          if (nu_id) "nu" else NULL,
                                          if (length(gamma)) names(gamma) else NULL)
  }

  out <- list(coefficients = list(beta = beta,
                                  kappa = kappa),
             fitted.values = list(mu = mu,
                                  sigma = sigma,
                                  lambda = lambda,
                                  nu = nu),
             links = list(mu.link = mu.link, sigma.link = sigma.link),
             margin = margin,
             logLik = opt$value,
             vcov = vcov,
             copula = copula,
             eta = eta,
             association = association,
             gamma = gamma,
             nobs = n,
             optim_params = opt,
             y = y, X = X, Z = Z)

  out

}




