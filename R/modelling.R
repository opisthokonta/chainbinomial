

# Link functions.
cb_link <- function(x, link){

  if (link == 'identity'){
    x
  }else if (link == 'logit'){
    log(x / (1-x))
  } else if (link == 'log'){
    log(x)
  } else if (link == 'cloglog'){
    log(-log(1-x))
  }

}

# Inverse link functions.
cb_invlink <- function(x, link){

  if (link == 'identity'){
    x
  }else if (link == 'logit'){
    1 / (1+exp(-x))
  } else if (link == 'log'){
    exp(x)
  } else if (link == 'cloglog'){
    1 - exp(-exp(x))
  }

}

# Objective function used with optim.
cb_reg_obj <- function(par, xmat, y, s0, i0, generations, link){
  sar_hat <- cb_invlink(as.numeric(xmat %*% par), link = link)
  negloglok_cb(sar = sar_hat, infected = y, s0 = s0, i0 = i0, generations = generations,
               transform_inv_logit = FALSE)
}


# Initial parameters (regression coefficients).
# Assumes that the first coefficient is the intercept.
initial_params <- function(y, s0, x, link){

  yy <- cb_link(mean(y/s0)*0.8, link = link)
  npar <- ncol(x)

  params <- c(yy, rep(0, npar-1))

  return(params)
}


#' Fitting models for Secondary Attack Rate with Chain Binomial response
#'
#' @param y numeric, the number of infected cases.
#' @param s0 numeric, the number of initial suceptibles.
#' @param x matrix of predictors (design matrix).
#' @param i0 numeric, number of initial infected. Default is 1.
#' @param generations numeric.
#' @param link Link function. Default is 'identitity'.
#'
#'
#'
#' @export
cbmod <- function(y, s0, x = NULL, i0 = 1, generations = Inf, link = 'identity'){

  stopifnot(length(y) == length(s0),
            is.null(x) | is.matrix(x),
            all(y >= 0),
            all(s0 >= 1))

  if (is.null(x)){
    x <- cbind(rep(1, length(y)))
  } else if (is.matrix(x)){
    stopifnot(nrow(x) == length(y),
              ncol(x) <= length(y))
  }

  warning_messages <- character(0)

  start_time <- Sys.time()

  par_init <- initial_params(y = y, s0 = s0, x = x, link = link)

  optim_res <- optim(par = par_init, fn = cb_reg_obj, method = 'BFGS', hessian=TRUE,
                     xmat = x, y = y, s0 = s0, i0 = i0, generations = generations, link = link)


  if (optim_res$convergence != 0){
    wmsg_convergence <- 'Did not converge (optim). Parameter estimates are unreliable.'
    warning_messages <- append(warning_messages, wmsg_convergence)
    warning(wmsg_convergence)
  }


  beta_hat <- optim_res$par
  names(beta_hat) <- colnames(x)

  vcov <- solve(optim_res$hessian)
  beta_se <- sqrt(diag(vcov))

  sar_hat <- cb_invlink(as.numeric(x %*% optim_res$par), link = link)



  end_time <- Sys.time()
  est_time <- difftime(end_time, start_time, units='secs')

  res <- list(parameters = beta_hat,
              se = beta_se,
              vcov = vcov,
              loglikelihood = -optim_res$value,
              npar = length(optim_res$par),
              sar_hat = sar_hat,
              link = link,
              warnings = warning_messages,
              est_time = est_time)

  class(res) <- 'cbmod'

  return(res)

}


#' @export
summary.cbmod <- function(object, ...){

  if (length(object$warnings) == 0){
    cat(sprintf('Model sucsessfully fitted in %.2f seconds\n\n', object$est_time))
  } else {
    cat(sprintf('Model fitted with warnings (%d). Results may be unreliable.\n\n', length(object$warnings)))
  }

  term_names_length <- max(nchar(names(object$parameters)))

  cat('Coefficients:\n')
  cat(sprintf('%-*s %8s %9s\n', term_names_length, '', 'Estimate', 'Std.error'))

  for (ii in 1:length(object$parameters)){
    cat(sprintf('%-*s % 8.3f % 9.3f\n', term_names_length, names(object$parameters)[ii], object$parameters[ii], object$se[ii]))
  }

}



