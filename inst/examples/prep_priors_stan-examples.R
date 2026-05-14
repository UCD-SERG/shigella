## Example: prep_priors_stan()
##
## Return the prior hyperparameters for the Chapter 2 Stan model.

## Default priors for the Kronecker correlated model (Chapter 2)
priors <- prep_priors_stan(model = "model_2")

str(priors)
