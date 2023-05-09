

# Inverse logit function. Used to map real numbe line to (0,1), to help optim.
inv_logit <- function(x){
  1 / (1+exp(-x))
}


# The negative chain binomial log-likelihood, for use with optim.
negloglok_cb <- function(sar, infected, s0, i0 , generations, transform_inv_logit = TRUE){

  # The parameter should be transformed when using this function for optimization (estimation),
  # but not when computing the hessian (for standard errors) or when the parameter
  # has already been transformed (in regression modelling).
  if (transform_inv_logit){
    sar <- inv_logit(sar)
  } else {
    if (any(sar < 0) | any(sar > 1)){
      return(Inf)
    }
  }

  -sum(log(dchainbinom(x = infected, s0 = s0, sar = sar, i0 = i0, generations = generations)))

}


#' Estimate Secondary Attack Rate of the Chain Binomial Model
#'
#' Given data on the number of infected after a number of generation, initial
#' number of susceptible, and initial number of infected, estimate the
#' secondary attack rate (SAR) using maximum likelihood.
#'
#' @param infected numeric.
#' @param s0 numeric.
#' @param i0 numeric.
#' @param generations numeric.
#'
#'@export
estimate_sar <- function(infected, s0, i0 = 1, generations=Inf){

  stopifnot(all(infected <= s0),
            all(infected >= 0))


  sar_init <- 0.5
  optim_res <- optim(par = sar_init,
                     fn = negloglok_cb, method = 'L-BFGS-B',
                     hessian = FALSE,
                     infected = infected, s0 = s0, i0 = i0, generations = generations)


  sar_hat <- inv_logit(optim_res$par)

  he <- numDeriv::hessian(negloglok_cb, x = sar_hat,
                          infected = infected, s0 = s0, i0 = i0,
                          generations = generations,
                          transform_inv_logit = FALSE)
  sar_se <- sqrt(as.numeric(solve(he)))


  res <- list(sar_hat = sar_hat,
              se = sar_se)

  class(res) <- 'sar'

  return(res)
}




confint.sar <- function(object, level = 0.95){

  upr <- min(object$sar_hat + (object$se * qnorm((1-level)/2, lower.tail = FALSE)) , 1)
  lwr <- max(object$sar_hat + (object$se * qnorm((1-level)/2, lower.tail = TRUE)), 0)

  cn <- sprintf('%.1f %%', c(100 * (1-level)/2, 100 * (1-((1-level)/2)) ))

  res <- c(lwr, upr)
  names(res) <- cn
  return(res)
}


# my_estimate <- list(sar_hat = 0.15, se = 0.03)
# class(my_estimate) <- 'sar'
# confint(my_estimate)





