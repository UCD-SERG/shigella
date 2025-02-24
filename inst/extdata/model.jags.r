model {
 for (subj in  1:nsubj) {
  for (cur_antigen_iso in 1:n_antigen_isos) { # `n_antigen_isos` is the number of biomarkers being modeled
   # 
     # beta is called `mu` in Teunis et al Epidemics 2016; 
     # it is the antibody growth rate during the active infection
     # this expression corresponds to equation 17 in that paper
     beta[subj, cur_antigen_iso] <- 
       log(
         y1[subj,cur_antigen_iso] / y0[subj,cur_antigen_iso]
         ) / 
       t1[subj,cur_antigen_iso]
     
    # `nsmpl` is the number of observations per `subject` 
   for(obs in 1:nsmpl[subj]) {
      
     # this is `log(y(t))` in the paper, before Gaussian noise is added
     mu.logy[subj, obs, cur_antigen_iso] <- ifelse(
        
        # `step(x)` returns 1 if x >= 0;
        # here we are determining which phase of infection we are in; 
        # active or recovery;
        # `smpl.t` is the time when the blood sample was collected, 
        # relative to estimated start of infection;
        # so we are determining whether the current observation is after `t1` 
        # the time when the active infection ended.
        step(t1[subj,cur_antigen_iso] - smpl.t[subj,obs]), 
        
        ## active infection period:
        # this is equation 15, case t <= t_1, but on a logarithmic scale
        log(y0[subj,cur_antigen_iso]) + (beta[subj,cur_antigen_iso] * smpl.t[subj,obs]),
        
        ## recovery period:
        # this is equation 15, case t > t_1
        1 / (1 - shape[subj,cur_antigen_iso]) *
           log(
              # this is `log{y_1^(1-r)}`; 
              # the exponent cancels out with the factor outside the log
              y1[subj, cur_antigen_iso]^(1 - shape[subj, cur_antigen_iso]) - 
                 
               # this is (1-r); not sure why switched from paper  
              (1 - shape[subj,cur_antigen_iso]) *
                
                  # (there's no missing y1^(r-1) term here; the math checks out)
                 
                 # alpha is `nu` in Teunis 2016; the "decay rate" parameter
                alpha[subj,cur_antigen_iso] *
                 
                 # this is `t - t_1`
                 (smpl.t[subj,obs] - t1[subj,cur_antigen_iso])))
     
     # we are fitting a loglinear model: log(Y) ~ N(mu, sigma^2)
     # this is the likelihood
     logy[subj,obs,cur_antigen_iso] ~ dnorm(mu.logy[subj,obs,cur_antigen_iso], prec.logy[cur_antigen_iso])
   }
    
  # these are random effects 
   y0[subj,cur_antigen_iso]    <- exp(par[subj,cur_antigen_iso,1])
   y1[subj,cur_antigen_iso]    <- y0[subj,cur_antigen_iso] + exp(par[subj,cur_antigen_iso,2]) # par[,,2] must be log(y1-y0)
   t1[subj,cur_antigen_iso]    <- exp(par[subj,cur_antigen_iso,3])
   alpha[subj,cur_antigen_iso] <- exp(par[subj,cur_antigen_iso,4]) # `nu` in the paper
   shape[subj,cur_antigen_iso] <- exp(par[subj,cur_antigen_iso,5]) + 1 # `r` in the paper
   
   # `n_params` is the number of model parameters; y0, y1, t1, alpha (aka nu), and r
   # this is the prior distribution
   par[subj, cur_antigen_iso, 1:n_params] ~ dmnorm(mu.par[cur_antigen_iso,], prec.par[cur_antigen_iso,,])
  }
 }
   
 # hyperpriors   
 for(cur_antigen_iso in 1:n_antigen_isos) {
    
  mu.par[cur_antigen_iso, 1:n_params] ~ dmnorm(mu.hyp[cur_antigen_iso,], prec.hyp[cur_antigen_iso,,])
  prec.par[cur_antigen_iso, 1:n_params, 1:n_params] ~ dwish(omega[cur_antigen_iso,,], wishdf[cur_antigen_iso])
  prec.logy[cur_antigen_iso] ~ dgamma(prec.logy.hyp[cur_antigen_iso,1], prec.logy.hyp[cur_antigen_iso,2])
 }
}
