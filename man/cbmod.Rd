% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelling.R
\name{cbmod}
\alias{cbmod}
\title{Fitting models for Secondary Attack Rate with Chain Binomial response}
\usage{
cbmod(
  y,
  s0,
  x = NULL,
  i0 = 1,
  generations = Inf,
  link = "identity",
  optim_method = "BFGS"
)
}
\arguments{
\item{y}{numeric, the number of infected cases.}

\item{s0}{numeric, the number of initial susceptibles.}

\item{x}{matrix of predictors (design matrix).}

\item{i0}{numeric, number of initial infected. Default is 1.}

\item{generations}{numeric.}

\item{link}{Link function. Default is 'identity'.}

\item{optim_method}{Optimization method used by optim.}
}
\value{
A list of class \code{cbmod} with the following components:
\itemize{
\item \code{parameters} The point estimate of the regression coefficients.
\item \code{se} Standard error of the regression coefficient estimates.
\item \code{vcov} Variance-Covariance matrix of the regression coefficient estimates.
\item \code{p_values} P-values of the null hypothesis that the regression regression coefficient estimate is 0.
\item \code{loglikelihood} the log likelihood value at the point estimate.
\item \code{npar} Number of parameters.
\item \code{sar_hat} Vector of fitted secondary attack rates.
\item \code{fitted_values} Vector of expected outbreak size (final attack rate).
\item \code{link} Link function used by the regression model.
\item \code{null_model} = Null model, fitted with \code{\link[=estimate_sar]{estimate_sar()}}. This is equivalent to an intercept only model.
\item \code{warnings} Warning_messages,
\item \code{est_time}: Time used to fit the model.
\item \code{omitted_values} Vector indicating data points that were ignored during estimation because of missing values.
}
}
\description{
Fitting models for Secondary Attack Rate with Chain Binomial response
}
\details{
The following link functions are available: \code{identity}, \code{log}, \code{logit}, and \code{cloglog}.
}
\examples{
set.seed(234)
mydata <- data.frame(infected = rchainbinom(n = 15, s0 = 5, sar = 0.2,
   i0 = 1, generations = Inf),
   s0 = 5, i0 = 1, generations = Inf)
xmat <- model.matrix(~ 1, data = mydata)
res <- cbmod(y = mydata$infected, s0 = mydata$s0, x = xmat, i0 = mydata$i0,
   generations = mydata$generations)
summary(res)


}
\seealso{
Methods for \code{cbmod} objects:
\itemize{
\item \code{\link[=summary.cbmod]{summary.cbmod()}}
\item \code{\link[=predict.cbmod]{predict.cbmod()}}
\item \code{\link[=coef.cbmod]{coef.cbmod()}}
\item \code{\link[=confint.cbmod]{confint.cbmod()}}
\item \code{\link[=vcov.cbmod]{vcov.cbmod()}}
\item \code{\link[=tidy.cbmod]{tidy.cbmod()}}
\item \code{\link[=glance.cbmod]{glance.cbmod()}}
}
}
