## Example: prep_priors_stan()
##
## Return the prior hyperparameters for the Kronecker correlated Stan model.

priors <- prep_priors_stan(model = "model_2")

str(priors)
