

# Chain probabilities ----

# These should work, and return a single number (a probability).
cp_1 <- chain_prob(x = c(1,2,3), s0 = 6, sar = 0.1)
cp_2 <- chain_prob(x = c(2,1,0), s0 = 4, sar = 0.1)

cp_1b <- chain_prob(x = c(1,2,3), s0 = 6, sar = 0.8)
cp_2b <- chain_prob(x = c(2,1,0), s0 = 4, sar = 0.8)

# Impossible chain, Should return probability 0.
cp_0a <- chain_prob(x = c(0,1,0), s0 = 4, sar = 0.8)
cp_0b <- chain_prob(x = c(1,0,1), s0 = 4, sar = 0.8)


# Chain of length 1 (index case only) is not defined. Shoudl give NA.
cp_na1 <- chain_prob(x = c(0), s0 = 4, sar = 0.8)
cp_na2 <- chain_prob(x = c(1), s0 = 4, sar = 0.8)


test_that("Chain probabilities", {
  expect_true(is.numeric(cp_1))
  expect_true(length(cp_1) == 1)

  expect_true(is.numeric(cp_2))
  expect_true(length(cp_2) == 1)

  expect_true(is.numeric(cp_1b))
  expect_true(length(cp_2b) == 1)

  expect_true(is.numeric(cp_0a))
  expect_true(length(cp_0a) == 1)
  expect_true(cp_0a == 0)

  expect_true(is.numeric(cp_0b))
  expect_true(length(cp_0b) == 1)
  expect_true(cp_0b == 0)

})



# Generating all scenarios ----

sc1 <- all_scenarios(target = 5, steps = 2)
sc2 <- all_scenarios(target = 5, steps = 1)
sc3 <- all_scenarios(target = 10, steps = 7)

test_that("All scenarions", {

  expect_true(nrow(sc1) == 2) # nrow == steps.
  expect_true(all(colSums(sc1) == 5)) # columns sum to target.

  expect_true(nrow(sc2) == 1)
  expect_true(all(colSums(sc2) == 5))

  expect_true(nrow(sc3) == 7)
  expect_true(all(colSums(sc3) == 10))

  # Same number of columns when the number of steps >= target.
  expect_true(ncol(all_scenarios(target = 5, steps = 5)) == ncol(all_scenarios(target = 5, steps = 6)))
  expect_true(ncol(all_scenarios(target = 2, steps = 2)) == ncol(all_scenarios(target = 2, steps = 4)))

  # Some checks that the scenarios sum to the correct number.
  expect_true(all(colSums(all_scenarios(target = 5, steps = 1)) == 5))
  expect_true(all(colSums(all_scenarios(target = 5, steps = 4)) == 5))
  expect_true(all(colSums(all_scenarios(target = 5, steps = 5)) == 5))
  expect_true(all(colSums(all_scenarios(target = 5, steps = 8)) == 5))

  expect_true(all(colSums(all_scenarios(target = 1, steps = 1)) == 1))
  expect_true(all(colSums(all_scenarios(target = 1, steps = 5)) == 1))

  expect_true(all(colSums(all_scenarios(target = 0, steps = 1)) == 0))
  expect_true(all(colSums(all_scenarios(target = 0, steps = 5)) == 0))


})



# Chain binomial PMF. ----

dcb_1_g1 <- dchainbinom(x = 0:5, s0 = 5, sar = 0.11, generations = 1)
dcb_1_g2 <- dchainbinom(x = 0:5, s0 = 5, sar = 0.11, generations = 2)
dcb_1_g3 <- dchainbinom(x = 0:5, s0 = 5, sar = 0.11, generations = 3)
dcb_1_g4 <- dchainbinom(x = 0:5, s0 = 5, sar = 0.11, generations = 4)
dcb_1_g5 <- dchainbinom(x = 0:5, s0 = 5, sar = 0.11, generations = 5)
dcb_1_g6 <- dchainbinom(x = 0:5, s0 = 5, sar = 0.11, generations = 6)
dcb_1_ginf <- dchainbinom(x = 0:5, s0 = 5, sar = 0.11, generations = Inf)

tol_sum_to_1 <- 2e-15

test_that("PMF is ok", {

  expect_true(abs(sum(dcb_1_g1) - 1) < tol_sum_to_1)
  expect_true(abs(sum(dcb_1_g2) - 1) < tol_sum_to_1)
  expect_true(abs(sum(dcb_1_g3) - 1) < tol_sum_to_1)
  expect_true(abs(sum(dcb_1_g4) - 1) < tol_sum_to_1)
  expect_true(abs(sum(dcb_1_g5) - 1) < tol_sum_to_1)
  expect_true(abs(sum(dcb_1_g6) - 1) < tol_sum_to_1)
  expect_true(abs(sum(dcb_1_ginf) - 1) < tol_sum_to_1)


  expect_true(all(dcb_1_g1 >= 0))
  expect_true(all(dcb_1_g2 >= 0))
  expect_true(all(dcb_1_g3 >= 0))
  expect_true(all(dcb_1_g4 >= 0))
  expect_true(all(dcb_1_g5 >= 0))
  expect_true(all(dcb_1_g6 >= 0))
  expect_true(all(dcb_1_ginf >= 0))


  expect_true(all(dcb_1_g5 == dcb_1_g6))
  expect_true(all(dcb_1_g5 == dcb_1_ginf))

  # Probability of more infected than s0 should be 0
  expect_true(dchainbinom(x = 0:6, s0 = 5, sar = 0.11, generations = 1)[7] == 0)
  expect_true(dchainbinom(x = 0:6, s0 = 5, sar = 0.11, generations = 3)[7] == 0)
  expect_true(dchainbinom(x = 0:6, s0 = 5, sar = 0.11, generations = Inf)[7] == 0)
  expect_true(dchainbinom(x = 0:6, s0 = 5, sar = 0.11, generations = Inf)[2] != 0)

})


# Additional check of sum to 1


check_sum_to_1 <- function(s0, sar, g, i0 = 1){

  ss <- sum(dchainbinom(x = 0:s0, s0 = s0, i0=i0, sar= sar, g = g))

  if (ss == 1){
    return(TRUE)
  } else {
    abs(ss-1) < 1e-15
  }
}


test_that("PMF sum to 1", {

  expect_true(check_sum_to_1(s0 = 3, sar= 0.1, g = 1))
  expect_true(check_sum_to_1(s0 = 8, sar= 0.3, g = 1))
  expect_true(check_sum_to_1(s0 = 2, sar= 0.1, g = 1))
  expect_true(check_sum_to_1(s0 = 1, sar= 0.2, g = 1))
  expect_true(check_sum_to_1(s0 = 0, sar= 0.1, g = 1))

  expect_true(check_sum_to_1(s0 = 3, sar = 0.1, g = 3))
  expect_true(check_sum_to_1(s0 = 8, sar= 0.3, g = 3))
  expect_true(check_sum_to_1(s0 = 2, sar= 0.1, g = 3))
  expect_true(check_sum_to_1(s0 = 1, sar= 0.2, g = 3))
  expect_true(check_sum_to_1(s0 = 0, sar= 0.1, g = 3))

  expect_true(check_sum_to_1(s0 = 3, sar = 0.1, g = 8))
  expect_true(check_sum_to_1(s0 = 8, sar= 0.3, g = 8))
  expect_true(check_sum_to_1(s0 = 2, sar= 0.1, g = 8))
  expect_true(check_sum_to_1(s0 = 1, sar= 0.2, g = 8))
  expect_true(check_sum_to_1(s0 = 0, sar= 0.1, g = 8))


  # with i0 = 2.
  expect_true(check_sum_to_1(s0 = 3, sar= 0.1, g = 1, i0=2))
  expect_true(check_sum_to_1(s0 = 8, sar= 0.3, g = 1, i0=2))
  expect_true(check_sum_to_1(s0 = 2, sar= 0.1, g = 1, i0=2))
  expect_true(check_sum_to_1(s0 = 1, sar= 0.2, g = 1, i0=2))
  expect_true(check_sum_to_1(s0 = 0, sar= 0.1, g = 1, i0=2))

  expect_true(check_sum_to_1(s0 = 3, sar = 0.1, g = 3, i0=2))
  expect_true(check_sum_to_1(s0 = 8, sar= 0.3, g = 3, i0=2))
  expect_true(check_sum_to_1(s0 = 2, sar= 0.1, g = 3, i0=2))
  expect_true(check_sum_to_1(s0 = 1, sar= 0.2, g = 3, i0=2))
  expect_true(check_sum_to_1(s0 = 0, sar= 0.1, g = 3, i0=2))

  expect_true(check_sum_to_1(s0 = 3, sar = 0.1, g = 8, i0=2))
  expect_true(check_sum_to_1(s0 = 8, sar= 0.3, g = 8, i0=2))
  expect_true(check_sum_to_1(s0 = 2, sar= 0.1, g = 8, i0=2))
  expect_true(check_sum_to_1(s0 = 1, sar= 0.2, g = 8, i0=2))
  expect_true(check_sum_to_1(s0 = 0, sar= 0.1, g = 8, i0=2))

})


# PMF when g = 1. ----

# Function to test that when the numbers of genereations = 1,
# it should be the same as the ordinary binomial model.
compare_g1_binom <- function(s0, sar, tol = 0.0000001){

  probvec1 <- dchainbinom(x = 0:s0, s0 = s0, sar = sar, generations = 1)
  probvec1_binom <- dbinom(x = 0:s0, size = s0, prob = sar)

  all((probvec1 - probvec1_binom < tol))

}


test_that("PMF when g = 1", {

  expect_true(compare_g1_binom(s0 = 3, sar = 0.1))
  expect_true(compare_g1_binom(s0 = 8, sar = 0.3))
  expect_true(compare_g1_binom(s0 = 2, sar = 0.1))
  expect_true(compare_g1_binom(s0 = 1, sar = 0.2))
  expect_true(compare_g1_binom(s0 = 0, sar = 0.1))

  expect_true(dchainbinom(x = 1, s0 = 0, sar = 0.1, generations = 1) == 0)

})



# Prob for x = 0 ----
# The probability for x = 0 should be the same as the ordinary binomial, when i0=1,
# but not neccecarily if i0 > 1.

# Explicit formula for x=0, also implemented in the dchainbinom function.
prob0 <- function(sar, s0, i0){
  (1 - sar)^(i0*s0)
}


test_that("P(x=0)", {

  expect_true(prob0(sar = 0.1, s0=1, i0 = 1) == dbinom(x = 0, size = 1, prob = 0.1))
  expect_true(prob0(sar = 0.5, s0=1, i0 = 1) == dbinom(x = 0, size = 1, prob = 0.5))

  expect_true(prob0(sar = 0.1, s0=2, i0 = 1) == dbinom(x = 0, size = 2, prob = 0.1))
  expect_true(prob0(sar = 0.5, s0=2, i0 = 1) == dbinom(x = 0, size = 2, prob = 0.5))

  expect_true(prob0(sar = 0.1, s0=4, i0 = 1) == dbinom(x = 0, size = 4, prob = 0.1))
  expect_true(prob0(sar = 0.5, s0=4, i0 = 1) == dbinom(x = 0, size = 4, prob = 0.5))

  expect_true(prob0(sar = 0.1, s0=2, i0 = 1) == dchainbinom(x = 0, s0 = 2, sar = 0.1, i0 = 1, generations = Inf))
  expect_true(prob0(sar = 0.1, s0=2, i0 = 1) == dchainbinom(x = 0, s0 = 2, sar = 0.1, i0 = 1, generations = 1))

  expect_true(prob0(sar = 0.1, s0=4, i0 = 1) == dchainbinom(x = 0, s0 = 4, sar = 0.1, i0 = 1, generations=Inf))
  expect_true(prob0(sar = 0.1, s0=4, i0 = 1) == dchainbinom(x = 0, s0 = 4, sar = 0.1, i0 = 1, generations=2))

  expect_true(dchainbinom(x = 0, s0 = 4, sar = 0.1, i0 = 1, generations = 2) == dchainbinom(x = 0, s0 = 4, sar = 0.1, i0=  1, generations = Inf))


  expect_false(prob0(sar = 0.1, s0=1, i0 = 2) == dbinom(x = 0, size = 1, prob = 0.1))
  expect_false(prob0(sar = 0.5, s0=1, i0 = 2) == dbinom(x = 0, size = 1, prob = 0.5))

  expect_false(prob0(sar = 0.1, s0=3, i0 = 2) == dbinom(x = 0, size = 3, prob = 0.1))
  expect_false(prob0(sar = 0.5, s0=3, i0 = 2) == dbinom(x = 0, size = 3, prob = 0.5))

  expect_false(dchainbinom(x = 0, s0 = 4, sar = 0.1, i0 = 2, generations = 2) == dchainbinom(x = 0, s0 = 4, sar = 0.1, i0=  1, generations = Inf))


})




# Compare PMF with simulation ----

# Function to approximate the chain binomial distribution after g generations
# using simulations.
#
# example:
# simulate_distribution(s0 = 5, sar = 0.1, g=3)
simulate_distribution <- function(s0, sar, g, i0 = 1){

  probs <- numeric(s0+1)
  simulated_cb <- rchainbinom(500000, s0 = s0, sar = sar, i0 = i0, g = g)
  sim_prop <- table(simulated_cb) / length(simulated_cb)
  probs[as.numeric(names(sim_prop)) + 1 ] <- sim_prop
  return(probs)
}


# Compares the PMF with simulated data.
compare_pmf_vs_simulation <- function(s0, sar, g, i0 = 1, tol = 0.001, maxtries = 3){

  for (ii in 1:maxtries){

    probvec <- dchainbinom(x = 0:s0, s0 = s0, sar = sar, i0 = i0, generations = g)
    probvec_sim <- simulate_distribution(s0 = s0, sar = sar, i0 = i0, g = g)

    if (all((probvec - probvec_sim < tol))){
      return(TRUE)
    } else {
      next
    }
  }

  return(FALSE)

}

if (FALSE){
  test_that("PMF vs simulation", {

    expect_true(compare_pmf_vs_simulation(s0 = 3, sar = 0.1, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 8, sar = 0.3, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 2, sar = 0.1, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 1, sar = 0.2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 0, sar = 0.1, g = 1))

    expect_true(compare_pmf_vs_simulation(s0 = 3, sar = 0.1, g = 3))
    expect_true(compare_pmf_vs_simulation(s0 = 8, sar = 0.3, g = 3))
    expect_true(compare_pmf_vs_simulation(s0 = 2, sar = 0.1, g = 3))
    expect_true(compare_pmf_vs_simulation(s0 = 1, sar = 0.2, g = 3))
    expect_true(compare_pmf_vs_simulation(s0 = 0, sar = 0.1, g = 3))

    expect_true(compare_pmf_vs_simulation(s0 = 3, sar = 0.1, g = 8))
    expect_true(compare_pmf_vs_simulation(s0 = 8, sar = 0.3, g = 8))
    expect_true(compare_pmf_vs_simulation(s0 = 2, sar = 0.1, g = 8))
    expect_true(compare_pmf_vs_simulation(s0 = 1, sar = 0.2, g = 8))
    expect_true(compare_pmf_vs_simulation(s0 = 0, sar = 0.1, g = 8))


    # with i0 = 2
    expect_true(compare_pmf_vs_simulation(s0 = 3, sar = 0.1, i0 = 2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 8, sar = 0.3, i0 = 2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 2, sar = 0.1, i0 = 2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 1, sar = 0.2, i0 = 2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 0, sar = 0.1, i0 = 2, g = 1))

    expect_true(compare_pmf_vs_simulation(s0 = 3, sar = 0.1, i0 = 2, g = 3))
    expect_true(compare_pmf_vs_simulation(s0 = 8, sar = 0.3, i0 = 2, g = 3))
    expect_true(compare_pmf_vs_simulation(s0 = 2, sar = 0.1, i0 = 2, g = 3))
    expect_true(compare_pmf_vs_simulation(s0 = 1, sar = 0.2, i0 = 2, g = 3))
    expect_true(compare_pmf_vs_simulation(s0 = 0, sar = 0.1, i0 = 2, g = 3))

    expect_true(compare_pmf_vs_simulation(s0 = 3, sar = 0.1, i0 = 2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 8, sar = 0.3, i0 = 2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 2, sar = 0.1, i0 = 2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 1, sar = 0.2, i0 = 2, g = 1))
    expect_true(compare_pmf_vs_simulation(s0 = 0, sar = 0.1, i0 = 2, g = 1))

  } )


}


# Expected value ----

# Expeced value should be same as for ordinary binomial when g=1.
ecb_vs_eb_g1 <- function(s0, sar){
  ecb <- echainbinom(s0 = s0, i0 = 1, sar = sar, generations = 1)
  eb <- s0*sar #sum(0:s0 * dbinom(x=0:s0, size=s0, prob = sar))

  abs(ecb - eb) < 0.00000001
}

test_that("Expected value g=1", {

  expect_true(ecb_vs_eb_g1(s0 = 3, sar=0.1))
  expect_true(ecb_vs_eb_g1(s0 = 3, sar=0.6))

  expect_true(ecb_vs_eb_g1(s0 = 2, sar=0.1))
  expect_true(ecb_vs_eb_g1(s0 = 2, sar=0.6))

  expect_true(ecb_vs_eb_g1(s0 = 9, sar=0.1))
  expect_true(ecb_vs_eb_g1(s0 = 9, sar=0.6))

})

# chain binomial expected value should be greater than ordinary binomial
# when g > 1.

ecb_vs_eb <- function(s0, sar, g){
  ecb <- echainbinom(s0 = s0, i0 = 1, sar = sar, generations = g)
  eb <- s0*sar

  ecb > eb
}


test_that("Expected value g>1", {

  expect_true(ecb_vs_eb(s0 = 3, sar=0.1, g=2))
  expect_true(ecb_vs_eb(s0 = 3, sar=0.6, g=2))
  expect_true(ecb_vs_eb(s0 = 2, sar=0.1, g=2))
  expect_true(ecb_vs_eb(s0 = 2, sar=0.6, g=2))
  expect_true(ecb_vs_eb(s0 = 9, sar=0.1, g=2))
  expect_true(ecb_vs_eb(s0 = 9, sar=0.6, g=2))

  expect_true(ecb_vs_eb(s0 = 3, sar=0.1, g=5))
  expect_true(ecb_vs_eb(s0 = 3, sar=0.6, g=5))
  expect_true(ecb_vs_eb(s0 = 2, sar=0.1, g=5))
  expect_true(ecb_vs_eb(s0 = 2, sar=0.6, g=5))
  expect_true(ecb_vs_eb(s0 = 9, sar=0.1, g=5))
  expect_true(ecb_vs_eb(s0 = 9, sar=0.6, g=5))

  expect_true(ecb_vs_eb(s0 = 3, sar=0.1, g=Inf))
  expect_true(ecb_vs_eb(s0 = 3, sar=0.6, g=Inf))
  expect_true(ecb_vs_eb(s0 = 2, sar=0.1, g=Inf))
  expect_true(ecb_vs_eb(s0 = 2, sar=0.6, g=Inf))
  expect_true(ecb_vs_eb(s0 = 9, sar=0.1, g=Inf))
  expect_true(ecb_vs_eb(s0 = 9, sar=0.6, g=Inf))

})

# The expected value should never be greater than s0.
ecb_le_s0 <- function(s0, sar, g){
  ecb <- echainbinom(s0 = s0, i0 = 1, sar = sar, generations = g)
  ecb < s0
}

test_that("Expected value < s0", {
  expect_true(ecb_le_s0(s0 = 3, sar=0.1, g=2))
  expect_true(ecb_le_s0(s0 = 3, sar=0.6, g=2))
  expect_true(ecb_le_s0(s0 = 2, sar=0.1, g=2))
  expect_true(ecb_le_s0(s0 = 2, sar=0.6, g=2))
  expect_true(ecb_le_s0(s0 = 9, sar=0.1, g=2))
  expect_true(ecb_le_s0(s0 = 9, sar=0.6, g=2))

  expect_true(ecb_le_s0(s0 = 3, sar=0.1, g=5))
  expect_true(ecb_le_s0(s0 = 3, sar=0.6, g=5))
  expect_true(ecb_le_s0(s0 = 2, sar=0.1, g=5))
  expect_true(ecb_le_s0(s0 = 2, sar=0.6, g=5))
  expect_true(ecb_le_s0(s0 = 9, sar=0.1, g=5))
  expect_true(ecb_le_s0(s0 = 9, sar=0.6, g=5))

  expect_true(ecb_le_s0(s0 = 3, sar=0.1, g=Inf))
  expect_true(ecb_le_s0(s0 = 3, sar=0.6, g=Inf))
  expect_true(ecb_le_s0(s0 = 2, sar=0.1, g=Inf))
  expect_true(ecb_le_s0(s0 = 2, sar=0.6, g=Inf))
  expect_true(ecb_le_s0(s0 = 9, sar=0.1, g=Inf))
  expect_true(ecb_le_s0(s0 = 9, sar=0.6, g=Inf))
})

# The expected value should equal s0 when sar = 1.
ecb_eq_s0_sar_eq_1 <- function(s0, i0, g){
  ecb <- echainbinom(s0 = s0, i0 = i0, sar = 1, generations = g)
  ecb == s0
}


test_that("Expected value == s0", {

  expect_true(ecb_eq_s0_sar_eq_1(s0  = 3, i0 = 1, g = 1))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 3, i0 = 1, g = 2))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 3, i0 = 1, g = 3))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 3, i0 = 1, g = Inf))

  expect_true(ecb_eq_s0_sar_eq_1(s0  = 3, i0 = 2, g = 1))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 3, i0 = 2, g = 2))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 3, i0 = 2, g = 3))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 3, i0 = 2, g = Inf))

  expect_true(ecb_eq_s0_sar_eq_1(s0  = 9, i0 = 2, g = 1))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 9, i0 = 2, g = 2))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 9, i0 = 2, g = 3))
  expect_true(ecb_eq_s0_sar_eq_1(s0  = 9, i0 = 2, g = Inf))
})







