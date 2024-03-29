

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

  if (any(sar_hat > 1 | sar_hat < 0)){
    return(Inf)
  }

  nll <- negloglok_cb(sar = sar_hat, infected = y, s0 = s0, i0 = i0, generations = generations,
                     transform_inv_logit = FALSE)

  return(nll)
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
#' @param s0 numeric, the number of initial susceptibles.
#' @param x matrix of predictors (design matrix).
#' @param i0 numeric, number of initial infected. Default is 1.
#' @param generations numeric.
#' @param link Link function. Default is 'identity'.
#' @param optim_method Optimization method used by optim.
#'
#' @details
#' The following link functions are available: `identity`, `log`, `logit`, and `cloglog`.
#'
#' @returns A list of class `cbmod`.
#'
#' @export
cbmod <- function(y, s0, x = NULL, i0 = 1, generations = Inf, link = 'identity',
                  optim_method = 'BFGS'){

  stopifnot(length(y) == length(s0),
            is.null(x) | is.matrix(x),
            all(y >= 0, na.rm = TRUE),
            all(s0 >= 1, na.rm = TRUE),
            all(i0 >= 1, na.rm = TRUE),
            all(generations >= 1, na.rm = TRUE),
            length(link) == 1,
            is.character(link),
            length(optim_method) == 1,
            optim_method %in% c('BFGS', 'L-BFGS-B', 'Nelder-Mead'))

  if (is.null(x)){
    x <- cbind(rep(1, length(y)))
  } else if (is.matrix(x)){
    stopifnot(nrow(x) == length(y),
              ncol(x) <= length(y))
  }

  if (length(i0) > 1){
    stopifnot(length(i0) == length(y))
  } else {
    i0 <- rep(i0, length(y))
  }

  if (length(generations) > 1){
    stopifnot(length(generations) == length(y))
  }else {
    generations <- rep(generations, length(y))
  }

  # Check for missing values.
  is_na_y <- is.na(y) | is.na(s0) | is.na(i0) | is.na(generations)
  is_na_x <- apply(x, FUN = function(x){any(is.na(x))}, MARGIN = 1)
  is_na <- is_na_y | is_na_x

  if (any(is_na)){
    y <- y[!is_na]
    s0 <- s0[!is_na]
    i0 <- i0[!is_na]
    generations <- generations[!is_na]
    x <- x[!is_na,]
  }

  warning_messages <- character(0)

  start_time <- Sys.time()

  par_init <- initial_params(y = y, s0 = s0, x = x, link = link)

  optim_res <- optim(par = par_init, fn = cb_reg_obj, method = optim_method, hessian = TRUE,
                     xmat = x, y = y, s0 = s0, i0 = i0, generations = generations, link = link)


  if (optim_res$convergence != 0){
    wmsg_convergence <- 'Did not converge (optim). Parameter estimates are unreliable.'
    warning_messages <- append(warning_messages, wmsg_convergence)
    warning(wmsg_convergence)
  }


  beta_hat <- optim_res$par
  names(beta_hat) <- colnames(x)

  vcov <- solve(optim_res$hessian)
  colnames(vcov) <- colnames(x)
  rownames(vcov) <- colnames(x)

  beta_se <- sqrt(diag(vcov))
  names(beta_se) <- colnames(x)

  sar_hat <- cb_invlink(as.numeric(x %*% optim_res$par), link = link)

  # Null model, useful for testing.
  sar_hat_0 <- estimate_sar(infected = y, s0 = s0, i0 = i0, generations = generations)

  # Fitted values (aka estimate of expected final attack rate).
  yhat <- echainbinom(s0 = s0, i0 = i0, sar = sar_hat, generations = generations)

  # p-values
  p_values <- pnorm(abs(beta_hat), mean = 0, sd = beta_se, lower.tail = FALSE) * 2


  end_time <- Sys.time()
  est_time <- difftime(end_time, start_time, units='secs')

  res <- list(parameters = beta_hat,
              se = beta_se,
              vcov = vcov,
              p_values = p_values,
              loglikelihood = -optim_res$value,
              npar = length(optim_res$par),
              sar_hat = sar_hat,
              fitted_values = yhat,
              link = link,
              null_model = sar_hat_0,
              warnings = warning_messages,
              est_time = est_time,
              omitted_values = which(is_na))

  class(res) <- 'cbmod'

  return(res)

}


# Internal function for comparing the fitted model to the intercept-only model.
# using likelihood ratio test.
cbmod_lr_test <- function(object){

  if (object$npar > 1){
    lr_test_stat <- 2*(object$loglikelihood - object$null_model$loglikelihood)
    lr_df <- object$npar-1
    lr_pv <- pchisq(lr_test_stat, df = lr_df, lower.tail = FALSE)

  } else {
    lr_test_stat <- NA
    lr_df <- 0
    lr_pv <- NA
  }

  res <- list(lr = lr_test_stat, df = lr_df, p.value = lr_pv)
  return(res)
}



#' @export
summary.cbmod <- function(object, ...){

  cat(sprintf('Chain Binomial model with %s link.\n', object$link))

  if (length(object$warnings) == 0){
    cat(sprintf('Model successfully fitted in %.2f seconds\n\n', object$est_time))
  } else {
    cat(sprintf('Model fitted with warnings (%d). Results may be unreliable.\n\n', length(object$warnings)))
  }


  cat(sprintf('Model log-likelihood: %17.1f\n', object$loglikelihood))
  cat(sprintf('Null log-likelihood: %18.1f\n', object$null_model$loglikelihood))

  lr_res <- cbmod_lr_test(object)

  cat(sprintf('Chisq (df = %d): %23.3f\n', lr_res$df, lr_res$lr))
  cat(sprintf('p-value: %30.3f\n', lr_res$p.value))

  cat('\n')

  term_names_length <- max(nchar(names(object$parameters)))

  cat('Coefficients:\n')
  cat(sprintf('%-*s %8s %9s %8s\n', term_names_length, '', 'Estimate', 'Std. Error', 'P-value'))

  for (ii in 1:length(object$parameters)){
    cat(sprintf('%-*s % 8.3f % 9.3f % 9.3f\n', term_names_length, names(object$parameters)[ii], object$parameters[ii], object$se[ii], object$p_values[ii]))
  }

}

#' @export
confint.cbmod <- function(object, level = 0.95){

  upr <- object$parameters + (object$se * qnorm((1-level)/2, lower.tail = FALSE))
  lwr <- object$parameters + (object$se * qnorm((1-level)/2, lower.tail = TRUE))

  cn <- sprintf('%.1f %%', c(100 * (1-level)/2, 100 * (1-((1-level)/2)) ))

  res <- cbind(lwr, upr)
  colnames(res) <- cn

  return(res)

}

#' @export
vcov.cbmod <- function(object){
  object$vcov
}

#' @export
coef.cbmod <- function(object){
  object$parameters
}



#' Predict Method for cbmod Fits
#'
#' @param object a fitted object of class inheriting from "cbmod".
#' @param x matrix of predictors (design matrix). Must have the same column names and order as the x matrix used to fit the model.
#' @param type the type of prediction, either 'link' (default) or 'sar'. The default is on the scale
#' of the linear predictors. 'sar' gives the predicted secondary attack rate, by transforming the linear
#'  predictors by the inverse link function used in the model fit.
#'
#' @return a vector of predictions.
#' @export
predict.cbmod <- function(object, x = NULL, type = 'link'){

  if(!identical(colnames(x), names(object$parameters))){
    stop('The column names and order in newdata must match those of names(object$parameters).')
  }

  stopifnot(length(type) == 1)
  stopifnot(type %in% c('link', 'response', 'sar'))


  eta_hat <- as.numeric(x %*% object$parameters)

  if (type %in% c('response', 'sar')){
    sar_hat <- cb_invlink(eta_hat, link = object$link)

    if (any(sar_hat < 0) | any(sar_hat > 1)){
      warning('predicted SAR smaller outside 0-1 range.')
    }

  }


  if (type == 'link'){
    return(eta_hat)
  } else if (type %in% c('response', 'sar')){
    return(sar_hat)
  }

}




