# the report wording is stable for each model family

    Code
      print(report_methods(fit))
    Output
      The analysis sample comprised # cases (of # supplied; # case had no
      observed indicator and was removed). Classes were defined by # binary
      indicators.
      
      A #-class latent class model was estimated by maximum likelihood with
      the EM algorithm (Dempster, Laird, & Rubin, #), with weak priors
      that keep estimates off the boundary of the parameter space (Galindo
      Garre & Vermunt, #): # pseudo-observation for the class proportions
      and # pseudo-observation for the item-response probabilities. The
      search used # random starts; the best log-likelihood was reached by #
      of them. The algorithm converged in # iterations.
      
      # of # indicator values (#.#%) were missing, in # indicator. Missing
      values were handled by full-information maximum likelihood, which
      assumes they are missing at random.
      
      The model's log-likelihood was -#.# with # free parameters (AIC =
      #.#, BIC = #.#, sample-size-adjusted BIC = #.#).
      
      Relative entropy was #.#. The smallest class held #.#% of the sample
      (about # cases). The average posterior probability of the assigned
      class ranged from #.# to #.# across classes. These describe how
      cleanly cases are classified; they are not criteria for the number of
      classes (Masyn, #).
      
      Models were estimated with the R package mixtureEM (version #.#.#.#)
      in R #.#.#. Random starts were drawn with random_state = #. The model
      is refitted by: fit_mixture(indicators = X, n_classes = #, measurement
      = "binary", n_init = #, random_state = #)
      
      References
      Dempster, A. P., Laird, N. M., & Rubin, D. B. (#). Maximum
        likelihood from incomplete data via the EM algorithm. Journal of the
        Royal Statistical Society: Series B, #(#), #-#.
      Galindo Garre, F., & Vermunt, J. K. (#). Avoiding boundary estimates
        in latent class analysis by Bayesian posterior mode estimation.
        Behaviormetrika, #(#), #-#. https://doi.org/#.#/bhmk.#.#
      Masyn, K. E. (#). Latent class analysis and finite mixture modeling.
        In T. D. Little (Ed.), The Oxford Handbook of Quantitative Methods
        (Vol. #, pp. #-#). Oxford University Press.

---

    Code
      print(report_methods(lcga))
    Output
      The analysis sample comprised # cases. The outcome was measured on #
      occasions (time scores #, #, #, #).
      
      A #-class latent class growth analysis (linear growth, no within-class
      growth-factor variance) was estimated by maximum likelihood with the EM
      algorithm (Dempster, Laird, & Rubin, #), with weak priors that keep
      estimates off the boundary of the parameter space (Galindo Garre &
      Vermunt, #; Chen, Tan, & Zhang, #): # pseudo-observation for the
      class proportions and # pseudo-observation for the class-specific
      variances. The search used # random starts; the best log-likelihood was
      reached by # of them. The algorithm converged in # iterations.
      
      The model's log-likelihood was -#.# with # free parameters (AIC =
      #.#, BIC = #.#, sample-size-adjusted BIC = #.#).
      
      Relative entropy was #.#. The smallest class held #.#% of the sample
      (about # cases). The average posterior probability of the assigned
      class ranged from #.# to #.# across classes. These describe how
      cleanly cases are classified; they are not criteria for the number of
      classes (Masyn, #).
      
      Models were estimated with the R package mixtureEM (version #.#.#.#)
      in R #.#.#. Random starts were drawn with random_state = #. The model
      is refitted by: fit_lcga(indicator = y, n_classes = #, times = #,
      family = "gaussian", n_init = #, random_state = #)
      
      References
      Chen, J., Tan, X., & Zhang, R. (#). Inference for normal mixtures in
        mean and variance. Statistica Sinica, #(#), #-#.
      Dempster, A. P., Laird, N. M., & Rubin, D. B. (#). Maximum
        likelihood from incomplete data via the EM algorithm. Journal of the
        Royal Statistical Society: Series B, #(#), #-#.
      Galindo Garre, F., & Vermunt, J. K. (#). Avoiding boundary estimates
        in latent class analysis by Bayesian posterior mode estimation.
        Behaviormetrika, #(#), #-#. https://doi.org/#.#/bhmk.#.#
      Masyn, K. E. (#). Latent class analysis and finite mixture modeling.
        In T. D. Little (Ed.), The Oxford Handbook of Quantitative Methods
        (Vol. #, pp. #-#). Oxford University Press.

---

    Code
      print(grolts_report(lcga))
    Output
      GRoLTS checklist (van de Schoot et al., #)
      
      Answered from the fit
      #. Is the metric of time used in the statistical model reported?
          Time scores #, #, #, # for occasions T#, T#, T#, T#.
      
      #c. Is a description provided of how missing data in the analyses were
          dealt with?
          No values were missing on the # occasions.
      
      #. Is the software mentioned?
          R package mixtureEM #.#.#.#, R #.#.#.
      
      #. Is information reported about the number of random start values and
          final iterations included?
          The search used # random starts; the best log-likelihood was
          reached by # of them. The algorithm converged in # iterations.
      
      #. Are the number of cases per class reported for each model?
          Final model only: #.#%, #.#%. Pass the comparison for every
          model.
      
      #. If classification of cases in a trajectory is the goal, is entropy
          reported?
          Relative entropy #.#.
      
      #. Are characteristics of the final class solution numerically
          described?
          Growth-factor means by class are in summary(fit); class shares
          #.#%, #.#%.
      
      #. Are the syntax files available?
          Models were estimated with the R package mixtureEM (version
          #.#.#.#) in R #.#.#. Random starts were drawn with random_state
          = #. The model is refitted by: fit_lcga(indicator = y, n_classes =
          #, times = #, family = "gaussian", n_init = #, random_state = #)
      
      
      You must supply
      #. Is information presented about the mean and variance of time within
          a wave?
          The model uses one fixed time score per occasion; report how much
          actual measurement times varied within each wave.
      
      #a. Is the missing data mechanism reported?
          The estimator assumes values are missing at random; argue whether
          that is plausible. mcar_test() tests the stronger MCAR assumption.
      
      #b. Is a description provided of what variables are related to
          attrition/missing data?
          Name the variables related to attrition or missingness.
      
      #. Is information about the distribution of the observed variables
          included?
          Describe the observed continuous outcome at each occasion (means,
          SDs, skew).
      
      #a. Are alternative specifications of within-class heterogeneity
          considered and clearly documented?
          This fit is a #-class latent class growth analysis (linear growth,
          no within-class growth-factor variance). Report which other
          within-class specifications (e.g., LCGA against GMM) were fitted
          and why this one was kept.
      
      #b. Are alternative specifications of the between-class differences in
          variance-covariance matrix structure considered and clearly
          documented?
          This fit has no growth-factor covariance; its residual variances
          are free across classes. Report which other structures were
          considered.
      
      #. Are alternative shape/functional forms of the trajectories
          described?
          This fit uses linear growth. Report the other shapes fitted, if
          any.
      
      #. Are the model comparison (and selection) tools described from a
          statistical perspective?
          Pass the compare_longitudinal() result to answer this from the
          fits.
      
      #. Are the total number of fitted models reported, including a
          one-class solution?
          Pass the compare_longitudinal() result to answer this from the
          fits.
      
      #a. Is a plot included with the estimated mean trajectories of the
          final solution?
          plot(fit) draws the estimated mean trajectories; include it.
      
      #b. Are plots included with the estimated mean trajectories for each
          model?
          Plot each compared model's mean trajectories, not only the final
          one's.
      
      #c. Is a plot included of the estimated means of the final model with
          the observed individual trajectories, by class?
          Overlay the observed individual trajectories on each class's
          estimated mean.
      
      
      Not applicable to this model
      #. If covariates have been used, can analyses still be replicated?
          No covariates were used.
      

---

    Code
      print(report_methods(lta))
    Output
      The analysis sample comprised # cases. Each occasion was measured by
      # binary indicators, over # occasions.
      
      A latent transition model with # statuses over # occasions was
      estimated by maximum likelihood with the EM algorithm (Dempster, Laird,
      & Rubin, #), with weak priors that keep estimates off the boundary
      of the parameter space (Galindo Garre & Vermunt, #): #
      pseudo-observation for the initial-status and transition probabilities
      and # pseudo-observation for the item-response probabilities. The
      search used # random starts; the best log-likelihood was reached by #
      of them. The algorithm converged in # iterations.
      
      The model's log-likelihood was -#.# with # free parameters (AIC =
      #.#, BIC = #.#, sample-size-adjusted BIC = #.#).
      
      Item-response probabilities were held equal across occasions, and the
      transition probabilities were estimated for each pair of occasions.
      
      Relative entropy was #.#. The smallest status held #.#% of the
      sample at the first occasion (about # cases). These describe how
      cleanly cases are classified; they are not criteria for the number of
      statuses (Masyn, #).
      
      Models were estimated with the R package mixtureEM (version #.#.#.#)
      in R #.#.#. Random starts were drawn with random_state = #. The model
      is refitted by: fit_lta(indicators = Xl, n_statuses = #, times = #,
      measurement = "binary", n_init = #, random_state = #, standard_errors =
      FALSE)
      
      References
      Dempster, A. P., Laird, N. M., & Rubin, D. B. (#). Maximum
        likelihood from incomplete data via the EM algorithm. Journal of the
        Royal Statistical Society: Series B, #(#), #-#.
      Galindo Garre, F., & Vermunt, J. K. (#). Avoiding boundary estimates
        in latent class analysis by Bayesian posterior mode estimation.
        Behaviormetrika, #(#), #-#. https://doi.org/#.#/bhmk.#.#
      Masyn, K. E. (#). Latent class analysis and finite mixture modeling.
        In T. D. Little (Ed.), The Oxford Handbook of Quantitative Methods
        (Vol. #, pp. #-#). Oxford University Press.

