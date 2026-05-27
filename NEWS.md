
# chainbinomial 0.2 (0.1.10)

There are a few changes that breaks backwards compatibility:
* The 'sar' argument has been renamed to 'prob' in all distribution functions. Using the 
  sar argument still works, but it will give a warning, and will be removed in 
  future versions. The reason is to conform to the more precise terminology
  of calling the parameter "transmission probability"" rather than "secondary attack rate".
* In the output from estimate_sar() the point estimate is now prob_hat, not sar_hat. 

New functionality: The Longini-Koopman model.
* Distributon function: dchainbinom2
* Estimate parameters with estimate_sar() function, set model = 'lk'.





# chainbinomial 0.1.5

* Initial CRAN submission.
