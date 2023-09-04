

# Inverse logit function. Used to map real number line to (0,1), to help optim.
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

  nll <- -sum(log(dchainbinom(x = infected, s0 = s0, sar = sar, i0 = i0, generations = generations)))

  if (is.nan(nll)){
    return(Inf)
  } else {
    return(nll)
  }
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
estimate_sar <- function(infected, s0, i0 = 1, generations=Inf, se = TRUE){

  stopifnot(all(infected <= s0),
            all(infected >= 0))


  sar_init <- 0.5
  optim_res <- optim(par = sar_init,
                     fn = negloglok_cb, method = 'L-BFGS-B',
                     hessian = FALSE,
                     infected = infected, s0 = s0, i0 = i0, generations = generations)


  sar_hat <- inv_logit(optim_res$par)

  if (se){
    he <- numDeriv::hessian(negloglok_cb, x = sar_hat,
                            infected = infected, s0 = s0, i0 = i0,
                            generations = generations,
                            transform_inv_logit = FALSE)
    sar_se <- sqrt(as.numeric(solve(he)))
  } else {
    sar_se <- NA
  }


  inp <- list(infected = infected,
                    s0 = s0,
                    i0 = i0,
                    generations = generations)

  res <- list(sar_hat = sar_hat,
              se = sar_se,
              loglikelihood = -optim_res$value,
              data = inp)

  class(res) <- 'sar'

  return(res)
}



# Function used by uniroot to find the value of the -2*log-likelihood-ratio
# that corresponds to the desired critical value of the chi-square distribution.
#
# x = parameter value
# infected, s0, i0, generations = data for the chain binomial likelihood.
# max_loglik = log-likelihood at the ML estimate.
# critical value = The value the -2*log-likelihood-ratio should be compared with,
#                    intended to be the critical value from the chi-square distribution.
obj_ci_wilks <- function(x, infected, s0, i0, generations, max_loglik, critical_value){

  # Log likelihood at the parameter value x.
  x_loglik <-  -negloglok_cb(sar = x, infected = infected, s0 = s0, i0 = i0,
                             generations = generations, transform_inv_logit = FALSE)

  # -2 log-likelihood ratio.
  nll2 <-  -2 * (x_loglik - max_loglik)

  critical_value - nll2

}


# Find the search interval for use with uniroot.
find_intervall_upr <- function(sh){

  if (sh >= 0.99){
    search_int <- c(0, 1)
  } else {
    search_int <- c(sh, 1)
  }

  return(search_int)

}

find_intervall_lwr <- function(sh){

  if (sh <= 0.01){
    search_int <- c(0, 1)
  } else {
    search_int <- c(0, sh)
  }

  return(search_int)

}


#'@export
confint.sar <- function(object, method = 'chisq', level = 0.95){

  stopifnot(method %in% c('chisq', 'normal'))

  if (method == 'normal'){

    # Compute Standard error if needed.
    if (is.na(object$se)){

      he <- numDeriv::hessian(negloglok_cb, x = object$sar_hat,
                              infected = object$data$infected,
                              s0 = object$data$s0,
                              i0 = object$data$i0,
                              generations = object$data$generations,
                              transform_inv_logit = FALSE)
      sar_se <- sqrt(as.numeric(solve(he)))

    } else{
      sar_se <- object$se
    }

    # compute CI using normal approximation.
    ci_lwr <- max(object$sar_hat + (sar_se * qnorm((1-level)/2, lower.tail = TRUE)), 0)
    ci_upr <- min(object$sar_hat + (sar_se * qnorm((1-level)/2, lower.tail = FALSE)) , 1)

  } else if (method == 'chisq'){

    plwr <- (1-level)/2
    pupr <- 1 - plwr


    critical_value_lower <- qchisq(plwr, df = 1, lower.tail = FALSE)

    uniroot_res_lwr <- uniroot(f = obj_ci_wilks, interval = find_intervall_lwr(object$sar_hat),
                       infected = object$data$infected, s0 = object$data$s0,
                       i0 = object$data$i0, generations = object$data$generations,
                       max_loglik = object$loglikelihood,
                       critical_value = critical_value_lower,
                       tol = 0.0000001)

    critical_value_upper <- qchisq(pupr, df = 1, lower.tail = TRUE)

    # browser()

    # -negloglok_cb(sar = 1, infected = object$data$infected, s0 = object$data$s0,
    #               i0 = object$data$i0, generations = object$data$generations, transform_inv_logit = FALSE)
    #
    # obj_ci_wilks(x = 0.9999, infected = object$data$infected, s0 = object$data$s0,
    #              i0 = object$data$i0, generations = object$data$generations,
    #              max_loglik = object$loglikelihood, critical_value_upper)

    uniroot_res_upr <- uniroot(f = obj_ci_wilks, interval = find_intervall_upr(object$sar_hat),
                       infected = object$data$infected, s0 = object$data$s0,
                       i0 = object$data$i0, generations = object$data$generations,
                       max_loglik = object$loglikelihood,
                       critical_value = critical_value_upper,
                       tol = 0.0000001)

    ci_lwr <- min(object$sar_hat, uniroot_res_lwr$root)
    ci_upr <- max(object$sar_hat, uniroot_res_upr$root)

  }


  cn <- sprintf('%.1f %%', c(100 * (1-level)/2, 100 * (1-((1-level)/2)) ))

  res <- c(ci_lwr, ci_upr)
  names(res) <- cn
  return(res)
}












