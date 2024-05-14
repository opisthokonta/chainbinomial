


#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy a cbmod object
#'
#' @param x A `cbmod` object.
#' @param conf.int Logical indicating whether or not to include a confidence interval in the tidied output. Defaults to FALSE.
#' @param conf.level The confidence level to use for the confidence interval if conf.int = TRUE. Defaults to 0.95, which corresponds to a
#'   95 percent confidence interval.
#' @param ... Unused.
#'
#' @export
tidy.cbmod <- function(x, conf.int=FALSE, conf.level = 0.95, ...){

  stopifnot(conf.level > 0,
            conf.level < 1)

  res <- tibble::tibble(term = names(x$parameters),
                        estimate = x$parameters,
                        std.error = x$se,
                        p.value = x$p_values)

  if (conf.int) {
    ci <- confint.cbmod(x, level = conf.level)
    res$conf.low <- ci[,1]
    res$conf.high <- ci[,2]
  }

  return(res)

}



#' @importFrom generics glance
#' @export
generics::glance

#' Glance at a cbmod object
#'
#' @param x A `cbmod` object.
#' @param ... Unused.
#'
#' @export
glance.cbmod <- function(x, ...){


  res <- tibble::tibble(logLik = x$loglikelihood,
                        npar = x$npar)

  return(res)

}

