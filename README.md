chainbinomial
================

The Chain Binomial model for infectious disease spread is especially
suitable for modelling of small outbreaks, such as outbreaks in
households. This package contains tools for analyzing data using the
Chain Binomial model and estimating the transmission probability. The
household transmission probability is defined as the probability that an
infected household member infects a susceptible household member.

Keywords: Infectious disease modelling, Chain Binomial model, Secondary
attack rate, Final attack rate, Final size distribution.

## Installation

On this package you will find the repository for the development version
of the package. It is recommended that you install and use the version
that is available from CRAN.
([link](https://cran.r-project.org/package=chainbinomial)). You can
install with this command:

``` r
install.packages("chainbinomial")
```

To install the latest development version:

``` r
install.packages("devtools")
devtools::install_github("opisthokonta/chainbinomial")
```

## Chain Binomial probabilities

Consider a household with 4 persons. A single household member becomes
infected by a contagious disease outside of the household, and the other
3 household members are susceptible to the disease. Assuming a
transmission probability of 0.23, we can compute the probability that 2
of the 3 susceptible household members becomes infected using the
`dchainbinom` function. The `dchainbinom` functions works similarly to
other discrete probability mass functions in R, such as the `dbinom` and
`dpois`.

``` r
library(chainbinomial)

dchainbinom(x = 2, s0 = 3, i0 = 1, prob = 0.23)
```

    ## [1] 0.1840275

We can also compute the entire final size distribution

``` r
dchainbinom(x = 0:3, s0 = 3, i0 = 1, prob = 0.23)
```

    ## [1] 0.4565330 0.2425560 0.1840275 0.1168835

Suppose instead that 2 of the 4 household members were infected
simultaneously outside of the household. The `i0` would then be 2. We
can again compute the final size distribution. Note that the number of
initial susceptible household members `s0` is now 2.

``` r
dchainbinom(x = 0:2, s0 = 2, i0 = 2, prob = 0.23)
```

    ## [1] 0.3515304 0.3717092 0.2767604

Now suppose that we don’t have observed the entire outbreak, but only a
single generation. The entire probability distribution after 1
generation can be computed using the `generations` argument. By default
the `generations` argument is Inf, meaning that the outbreak is assumed
to be completely observed.

``` r
dchainbinom(x = 0:3, s0 = 3, i0 = 1, prob = 0.23, generations = 1)
```

    ## [1] 0.456533 0.409101 0.122199 0.012167

### Simulating data

The `rchainbinom` function can be used to simulate data. Suppose we want
to simulate data on the number of infected household from 10 households,
with prob = 0.2 for the first 5, and prob = 0.4 for the last 5, all with
4 susceptible and one initial infected person. This can be done like
this:

``` r
set.seed(1)

rchainbinom(n = 10, prob = rep(c(0.2, 0.4), each = 5), s0 = 4, i0 = 1, generations = Inf)
```

    ##  [1] 0 0 3 4 2 1 4 4 4 3

The `rchainbinom` works similarly to the `rbinom` and `rpois` functions.

## Estimating the transmission probability

Suppose we have data on how many become infected in 20 households and we
want to estimate the transmission probability. The households may be of
different sizes and have different number of initial infectees. Lets
simulate some data with a know SAR = 0.3:

``` r
set.seed(123)

my_simulated_data <- data.frame(s0 = c(2,3,4,2,1,5,4,4,4,1,1,3,4,1,1,2,3,1,3,6),
                                i0 = c(1,1,1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1))

my_simulated_data$infected <- rchainbinom(n = nrow(my_simulated_data),
                                          prob = 0.3,
                                          s0 = my_simulated_data$s0,
                                          i0 = my_simulated_data$i0,
                                          generations = Inf)


my_simulated_data
```

    ##    s0 i0 infected
    ## 1   2  1        0
    ## 2   3  1        3
    ## 3   4  1        2
    ## 4   2  1        0
    ## 5   1  2        0
    ## 6   5  1        4
    ## 7   4  1        4
    ## 8   4  1        3
    ## 9   4  1        4
    ## 10  1  1        0
    ## 11  1  2        1
    ## 12  3  1        2
    ## 13  4  1        4
    ## 14  1  1        0
    ## 15  1  1        1
    ## 16  2  1        1
    ## 17  3  1        0
    ## 18  1  1        0
    ## 19  3  1        2
    ## 20  6  1        3

Now lets estimate the transmission probability using the `estimate_sar`
function

``` r
my_tp_estimate <- estimate_sar(infected = my_simulated_data$infected, 
                                 s0 = my_simulated_data$s0, 
                                 i0 = my_simulated_data$i0,
                                 generations = Inf)

my_tp_estimate$sar_hat
```

    ## NULL

We can also compute 95% confidence intervals

``` r
confint(my_tp_estimate)
```

    ##     2.5 %    97.5 % 
    ## 0.2341291 0.4493579

## Predictors of transmission, association analysis

Suppose the households differ in some systematic way and we want to see
if there are some factors that are associated with a larger of smaller
transmission probability. We can let the probability depend on a set of
predictors, similar to a Generalized Linear Model (GLM). The predictors
in this model would operate on the household level, not on the level of
individuals. One example of a predictor would be the strain or variant
of the infectious agent.

Lets simulate some data again, with a simple binary predictor called
strain_type. The transmission probability for strain_type = 0 is 0.2 and
for strain_type = 1 it is 0.5.

``` r
set.seed(10266)

# Same s0 and i0 as before.
my_simulated_data <- data.frame(s0 = c(2,3,4,2,1,5,4,4,4,1,1,3,4,1,1,2,3,1,3,6),
                                i0 = c(1,1,1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1),
                                strain_type = rep(c(0,1), each = 10))

my_simulated_data$infected <- rchainbinom(n = nrow(my_simulated_data),
                                          prob = 0.2 + my_simulated_data$strain_type*0.3,
                                          s0 = my_simulated_data$s0,
                                          i0 = my_simulated_data$i0)
```

Lets fit a model with strain_type as predictor using the `cbmod`
function. To use the formula interface that is common in R to specify
models, we can use the `model.matrix` function to make the X matrix that
is then passed on to `cbmod.`

Note that `s0`, `i0`, and `generations` should not be thought of as
predictor variables and should not be included in the X matrix. They are
also not response variables (the number of infected is the response),
but could perhaps be thought of as ‘nuisance data’.

``` r
xmat <- model.matrix(~ strain_type, data = my_simulated_data)

cbmod_res <- cbmod(y = my_simulated_data$infected, 
                   s0 = my_simulated_data$s0, 
                   i0 = my_simulated_data$i0,
                   generations = Inf,
                   x = xmat, 
                   link = 'identity')

summary(cbmod_res)
```

    ## Chain Binomial model with identity link.
    ## Model successfully fitted in 0.05 seconds
    ## 
    ## Model log-likelihood:             -19.6
    ## Null log-likelihood:              -24.0
    ## Chisq (df = 1):                   8.916
    ## p-value:                          0.003
    ## 
    ## Coefficients:
    ##             Estimate Std. Error  P-value
    ## (Intercept)    0.204     0.058     0.000
    ## strain_type    0.328     0.113     0.004

``` r
confint(cbmod_res)
```

    ##                  2.5 %    97.5 %
    ## (Intercept) 0.09064079 0.3178483
    ## strain_type 0.10600845 0.5498695

Here we used the identity link function, which is the default. This
gives the easiest interpretation of the coefficients, but will often not
work in more complicated models with more than one predictor or when the
predictor(s) are numerical rather than categorical. In that case you
should use `link='logit'`.

`tidy` and `glance` methods are also available for `cbmod` objects.

## Unobserved individuals

*This methodology is still under development so the function interface,
implementations, and underlying methodology might change.*

Sometimes there are individuals in the household whose infection status
is not known. This could be because they did not get tested, did not
consent to participate in the study or were excluded for some other
reason. These individuals will still contribute to the outbreak dynamics
within the household and their presence ought to be modelled and not
ignored, even if their infection status is unknown.

One way to deal with this is to assume an underlying chain binomial
model of the outbreak, and have a hypergeometric observational model on
top of that. The probability of observing x infected in a household of 5
initial susceptible individuals, where only 4 of them are observed
(i.e. 1 individual is not observed) can be calculated with the
`dcbhyper` function, using the `s0_obs` argument.

``` r
dcbhyper(x=0:4, s0 = 5, prob = 0.25, s0_obs = 4)
```

    ## [1] 0.2623329 0.1470408 0.1563287 0.2069753 0.2273222

There is also a function `ecbhyper` for calculating the expected value
and `estimate_sar_cbhyper` for estimating the sar parameter.

## Litterature

- Ludwig, D. (1975) Final Size Distributions for Epidemics. Mathematical
  Biosciences, 23, 33-46. <https://doi.org/10.1016/0025-5564(75)90119-4>

- Sharker Y, Kenah E (2021) Estimating and interpreting secondary attack
  risk: Binomial considered biased. PLoS Comput Biol 17(1): e1008601.
  <https://doi.org/10.1371/journal.pcbi.1008601>

- Lindstrøm JC, et al (2024) Estimating the household secondary attack
  rate with the Incomplete Chain Binomial model. In preparation.
  <https://doi.org/10.48550/arXiv.2403.03948>
