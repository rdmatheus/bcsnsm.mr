% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_association_matrices.R
\name{association}
\alias{association}
\alias{nonassociative}
\alias{arma}
\alias{matern}
\alias{cluster}
\title{Association Structures for the BCS-NSM Marginal Regression Models}
\usage{
nonassociative()

arma(p = 0L, q = 0L)

matern(D, smt = 0.5)

cluster(id, type = c("uniform", "unstructured", "ar1", "ma1"))
}
\arguments{
\item{p, q}{order of the autoregressive and moving average components, respectively, for the
\emph{ARMA(p, q)} association structure for time series.}

\item{D}{matrix with values of the distances between pairs of data
locations for spatial data.}

\item{smt}{value of the shape parameter of the Matérn association structure.
The default \code{smt = 0.5} corresponds to an exponential correlation
model.}

\item{id}{subject id for longitudinal/clustered data. This is a vector of
the same length of the number of observations. Please note that data
must be sorted in way that observations from the same cluster are
contiguous.}

\item{type}{a character string specifying the correlation structure among
groups for longitudinal/clustered data. At the
moment, the following are implemented:
\tabular{ll}{
\code{uniform}  \tab uniform or exchangeable. \cr
\code{unstructured}  \tab Unstructured \cr
\code{ar1}  \tab autoregressive of order 1. \cr
\code{ma1}  \tab moving average of order 1. \cr
}}
}
\value{
An \code{"association"} object, which consists of a list with the following components:
 \itemize{
 \item{\code{npar}:}{ the number of parameters related to the association structure.}
 \item{\code{start}:}{ a function (\code{function(y)}) that returns initial values for the 
     parameters of the association matrix, based on the data \code{y}.}
 \item{\code{Gamma}:}{ a function (\code{function(gamma, n)}) that returns the association matrix
     as a function of the parameter vector \code{gamma} and sample size \code{n}.}
 \item{\code{name}:}{ a character with the name of the association structure.}
 }
}
\description{
Association structures for the class of the BCS-NSM marginal regression models.
}
\details{
The association matrix plays an important role in describing dependence in the class of the 
    BCS-NSM marginal regression models, and the \code{"association"} objects contains all the
    essential information about it for model fitting in the \code{bcsnsm.mr} package.    

The functions are a direct adaptation of the functions of the \code{\link[gcmr]{gcmr}} package.
    The documentation of the original functions of the package can be seen
    at \code{\link[gcmr]{cormat.gcmr}} documentation. 
    
The current available association matrices are:
 \tabular{ll}{
 \bold{Function}  \tab \bold{Dependence structure}\cr
 \code{nonassociative()}  \tab non-associative structure, which provides independence only under
                               the Gaussian copula \cr
 \code{arma(p, q)}  \tab autoregressive moving average process. \cr
 \code{cluster(id, type)}  \tab longitudinal or clustered data. \cr
 \code{matern(D, smt)}  \tab spatial data. \cr
 }

The \code{nonassociative()} function can be specified with the Gaussian copula to model 
    independent data. For time series analysis, the association structure can be specified with 
    the function \code{arma(p, q)}, where \code{p} and \code{q} are the orders of the
    \emph{ARMA(p, q)} process. Spatial data can be described with the function
    \code{matern(D, shape)}, where \code{D} is a matrix with values of the distances between
    observations and \code{shape} is the value of the shape parameter of the Matérn correlation 
    function, which is, in general, fixed in the applications. The default \code{shape = 0.5}
    corresponds to an exponential association model between pairs of data locations. Finally,
    clustered or longitudinal data can be modeled with the function \code{cluster(id, type)},
    where \code{id} is a vector of the same length as the number of observations and \code{type}
    is a character which informs the type of correlation within each group with the following
    options: \code{"ar1"}, \code{"ma1"}, \code{"uniform"} and \code{"unstructured"}. Here, the 
    data must be ordered such that observations from the same individual (or group) are contiguous.
}
\references{
Guido Masarotto, & Cristiano Varin (2017). Gaussian Copula
    Regression in R. \emph{Journal of Statistical Software}, \bold{77}, 1--26.
}
\author{
Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
}
