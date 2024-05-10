#' @name association
#'
#' @title Association Structures for the BCS-NSM Marginal Regression Models
#'
#' @description Association structures for the class of the BCS-NSM marginal regression models.
#'
#' @param p,q order of the autoregressive and moving average components, respectively, for the
#'     \emph{ARMA(p, q)} association structure for time series.
#' @param id 	subject id for longitudinal/clustered data. This is a vector of
#'     the same length of the number of observations. Please note that data
#'     must be sorted in way that observations from the same cluster are
#'     contiguous.
#' @param type  a character string specifying the correlation structure among
#'     groups for longitudinal/clustered data. At the
#'     moment, the following are implemented:
#'     \tabular{ll}{
#'     \code{uniform}  \tab uniform or exchangeable. \cr
#'     \code{unstructured}  \tab Unstructured \cr
#'     \code{ar1}  \tab autoregressive of order 1. \cr
#'     \code{ma1}  \tab moving average of order 1. \cr
#'     }
#' @param D matrix with values of the distances between pairs of data
#'     locations for spatial data.
#' @param smt value of the shape parameter of the Matérn association structure.
#'     The default \code{smt = 0.5} corresponds to an exponential correlation
#'     model.
#'
#' @details
#' 
#' The association matrix plays an important role in describing dependence in the class of the 
#'     BCS-NSM marginal regression models, and the \code{"association"} objects contains all the
#'     essential information about it for model fitting in the \code{bcsnsm.mr} package.    
#'
#' The functions are a direct adaptation of the functions of the \code{\link[gcmr]{gcmr}} package.
#'     The documentation of the original functions of the package can be seen
#'     at \code{\link[gcmr]{cormat.gcmr}} documentation. 
#'     
#' The current available association matrices are:
#'  \tabular{ll}{
#'  \bold{Function}  \tab \bold{Dependence structure}\cr
#'  \code{nonassociative()}  \tab non-associative structure, which provides independence only under
#'                                the Gaussian copula \cr
#'  \code{arma(p, q)}  \tab autoregressive moving average process. \cr
#'  \code{cluster(id, type)}  \tab longitudinal or clustered data. \cr
#'  \code{matern(D, smt)}  \tab spatial data. \cr
#'  }
#' 
#' The \code{nonassociative()} function can be specified with the Gaussian copula to model 
#'     independent data. For time series analysis, the association structure can be specified with 
#'     the function \code{arma(p, q)}, where \code{p} and \code{q} are the orders of the
#'     \emph{ARMA(p, q)} process. Spatial data can be described with the function
#'     \code{matern(D, shape)}, where \code{D} is a matrix with values of the distances between
#'     observations and \code{shape} is the value of the shape parameter of the Matérn correlation 
#'     function, which is, in general, fixed in the applications. The default \code{shape = 0.5}
#'     corresponds to an exponential association model between pairs of data locations. Finally,
#'     clustered or longitudinal data can be modeled with the function \code{cluster(id, type)},
#'     where \code{id} is a vector of the same length as the number of observations and \code{type}
#'     is a character which informs the type of correlation within each group with the following
#'     options: \code{"ar1"}, \code{"ma1"}, \code{"uniform"} and \code{"unstructured"}. Here, the 
#'     data must be ordered such that observations from the same individual (or group) are contiguous.
#'
#'
#' @return
#'  
#'  An \code{"association"} object, which consists of a list with the following components:
#'  \itemize{
#'  \item{\code{npar}:}{ the number of parameters related to the association structure.}
#'  \item{\code{start}:}{ a function (\code{function(y)}) that returns initial values for the 
#'      parameters of the association matrix, based on the data \code{y}.}
#'  \item{\code{Gamma}:}{ a function (\code{function(gamma, n)}) that returns the association matrix
#'      as a function of the parameter vector \code{gamma} and sample size \code{n}.}
#'  \item{\code{name}:}{ a character with the name of the association structure.}
#'  }
#'
#' @references Guido Masarotto, & Cristiano Varin (2017). Gaussian Copula
#'     Regression in R. \emph{Journal of Statistical Software}, \bold{77}, 1--26.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export

## Non-associative ---------------------------------------------------
#' @rdname association
#' @export
nonassociative <- function() {
  ans <- list()
  
  ans$npar <- 0L
  
  ans$start <- function(x) NULL
  
  ans$Gamma <- function(gamma = NULL, n) diag(n)
  
  ans$name <- "Non-associative"
  
  ans
}


## ARMA(p,q) working correlation for time-series ---------------------------------------------------
#' @rdname association
#' @export
arma <- function(p = 0L, q = 0L) {

  ar <- if(p) 1:p else NULL
  ma <- if(q) 1:q + p else NULL

  ans <- list()

  ans$npar <- p + q

  ans$start <- function(y = NULL) {

    n <- length(y)

    if (p | q){
      gamma <- rep(0, p + q)
      names(gamma) <- c(if (p) paste0("ar", 1:p) else NULL, if (q) paste0("ma", 1:q) else NULL)
    } else {
      gamma <- NULL
    }

    gamma
  }

  ans$Gamma <- function(gamma, n){

    if (p | q){

      if((p && any(Mod(polyroot(c(1,-gamma[ar]))) < 1.01) ) || (q && any(Mod(polyroot(c(1, gamma[ma]))) < 1.01)))
        return(NULL)

      gamma <- stats::ARMAacf(gamma[ar], gamma[ma], n - 1)
      r <- seq(1, n)
      Matrix::as.matrix(outer(r, r, function(i, j) gamma[1 + abs(i - j)]))

    } else {

      diag(n)

    }


  }

  ans$name <- paste0("ARMA(", p, ", ", q,")")

  ans
}

## Matern working correlation for spatial data -----------------------------------------------------

## D is a distance matrix, smt is the smoothing parameter
#' @rdname association
#' @export
matern <- function(D, smt = 0.5) {
  ans <- list()
  ans$npar <- 1
  ans$start <- function(y = NULL) {
    gamma <- stats::median(D)
    names(gamma) <- c("gamma")
    attr(gamma,"lower") <- sqrt(.Machine$double.eps)
    gamma
  }

  ans$Gamma <- function(gamma, n){
    S <- try(Matrix::forceSymmetric(.matern(D, gamma, smt)), silent = TRUE)
    if(inherits(S, "try-error") | gamma < 0) NULL else S
  }

  ans$name <- "Matern spatial"
  ans
}

# Adapted from 'geoR' package
.matern <- function (u, phi, kappa)
{
  if (is.vector(u))
    u <- matrix(u)

  dimnames(u) <- list(NULL, NULL)

  uphi <- u/phi

  out <- Matrix::Matrix(NaN, dim(u)[1], dim(u)[2])

  id1 <- which(u > 0 & phi > 0 & kappa > 0, arr.ind = TRUE)
  id2 <- which(u <= 0 & phi > 0 & kappa > 0, arr.ind = TRUE)

  out[id1] <- exp(-(kappa - 1) * log(2) - lgamma(kappa) + kappa * log(uphi[id1]) + log(besselK(uphi[id1], kappa)))
  out[id2] <- 1
  out[which(u > 600 * phi, arr.ind = TRUE)] <- 0

  out
}

### Longitudinal/Clustered data working correlation ----------------------------
# It assumes that it is not possible that all the observations inside a cluster
# can be missing
#' @rdname association
#' @export
cluster <- function(id, type = c("uniform", "unstructured", "ar1", "ma1")) {

  type <- match.arg(type, c("uniform", "unstructured", "ar1", "ma1"))

  if(!length(rle(id)$values) == length(unique(id)))
    stop("data must be sorted in way that observations from the same cluster are contiguous")

  ng <- 1:length(unique(id))
  if (!(length(ng) > 1)) stop("only one strata")

  
  ans <- list(type = type, id = id)
  ans$npar <- if (type != "unstructured") 1 else choose(max(table(id)), 2)
  data <- data.frame(id = id)
  fn <- switch(type,
               "ar1" = function(g) nlme::corAR1(g, form = ~ 1 | id),
               "ma1" = function(g) nlme::corARMA(g, form = ~ 1 |id, p = 0, q = 1),
               "uniform" = function(g) nlme::corCompSymm(g, form = ~ 1 | id),
               "unstructured" = function(g) nlme::corSymm(g, form = ~ 1 | id))
  ans$start <- function(y = NULL) {
    np <-  if(type != "unstructured") 1 else choose(max(table(id)), 2)
    alpha <- rep(0, np)
    names(alpha) <- switch(type, "ar1" = "ar1", "ma1" = "ma1",
                           "uniform" = "alpha",
                           "unstructured" = paste("alpha", 1:ans$npar, sep = ""))
    eps <- sqrt(.Machine$double.eps)
    attr(alpha, "lower") <- rep(-1 + eps, np)
    attr(alpha, "upper") <- rep(1 - eps, np)
    alpha
  }
  ans$Gamma <- function(alpha, n) {
    q <- try(nlme::corMatrix(nlme::Initialize(fn(alpha), data = data)),
             silent = TRUE)
    if (inherits(q, "try-error")) return(NULL)
    g <- split(rep(TRUE, n), id)
    q <- try(lapply(ng, function(i) q[[i]][g[[i]],g[[i]]]), silent = TRUE)
    if (inherits(q, "try-error") ) NULL else as.matrix(Matrix::bdiag(q))
  }
  ans$name <- "Clustered"
  ans
}

