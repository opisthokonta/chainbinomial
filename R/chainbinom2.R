

# Longini Koopman Final Size Distribution -----

# For the final size distribution.
# Compute the p(I,I) factors from 0 to maxi, for a given sar and cpi
compute_factors_lk <- function(maxi, prob, cpi){

  stopifnot(length(maxi) == 1,
            maxi >= 0,
            length(prob) == 1,
            length(cpi) == 1)


  # The p(i,i) factors
  pii <- numeric(maxi+1)

  pii[1] <- 1 # p(0,0)

  if (maxi == 0){
    return(pii)
  }

  pii[2] <- cpi # p(1,1)

  if (maxi == 1){
    return(pii)
  }

  # Loop over all possible number of susceptible from 2 to maxi
  for (ss in 2:maxi){

    pvec <- numeric(ss+1)

    pvec[1] <- (1 - cpi)^ss # prob of 0 infected.

    # Loop over all number of infected, from 1 to ss-1.
    for (ii in 1:(ss-1)){
      idx <- ii+1
      pvec[idx] <- choose(ss, ii) * pii[ii+1] * ((1 - cpi)^(ss - ii)) * ((1-prob)^(ii * (ss - ii)))
    }

    pii[ss+1] <- max(0, 1 - sum(pvec))

  }

  return(pii)

}


#' @rdname dchainbinom
#' @export
dchainbinom2 <- function(x, s0, prob, cpi){

  # Check input.
  stopifnot(is.numeric(x) | is.logical(x),
            is.numeric(s0) | is.logical(s0),
            is.numeric(prob) | is.logical(prob),
            is.numeric(cpi) | is.logical(cpi),
            all(prob >= 0, na.rm = TRUE),
            all(prob <= 1, na.rm = TRUE),
            all(cpi >= 0, na.rm = TRUE),
            all(cpi <= 1, na.rm = TRUE),
            all(x >= 0, na.rm = TRUE),
            all(s0 >= 0, na.rm = TRUE))

  # Coerce to numericals, in case of logicals being provided.
  x <- as.numeric(x)
  s0 <- as.numeric(s0)
  prob <- as.numeric(prob)
  cpi <- as.numeric(cpi)

  # Combine input to matrix, to expand/recycle input data.
  inp <- as.matrix(cbind(x, s0, prob, cpi))

  n <- nrow(inp)
  res <- numeric(n)

  for (ii in 1:n){

    if (any(is.na(inp[ii,]))){

      res[ii] <- NA

    } else {

      if (inp[ii, 'x'] > inp[ii, 's0']){
        res[ii] <- 0
        next
      }

      pii <- compute_factors_lk(maxi = inp[ii, 'x'], prob = inp[ii, 'prob'], cpi = inp[ii, 'cpi'])

      BB <- (1 - inp[ii, 'cpi']) ^ (inp[ii, 's0'] - inp[ii, 'x'])
      QQ <- (1 - inp[ii, 'prob']) ^ (inp[ii, 'x'] * (inp[ii, 's0'] - inp[ii, 'x']))

      binom_coef <- choose(inp[ii, 's0'], inp[ii, 'x'])
      res[ii] <- binom_coef * pii[length(pii)] * BB * QQ

    }

  }

  return(res)

}



