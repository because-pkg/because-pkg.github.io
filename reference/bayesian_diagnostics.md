# Assumption Diagnostics for because Models

A suite of diagnostic functions for testing regression and DAG/SEM
assumptions in `because` models.

Diagnostics evaluate posterior parameter covariance geometry and
posterior predictive discrepancies.

Included diagnostics:

- [`check_collinearity`](https://because-pkg.github.io/because/reference/check_collinearity.md)
  /
  [`bvif`](https://because-pkg.github.io/because/reference/check_collinearity.md):
  Leamer's Bayesian Variance Inflation Factor and posterior parameter
  correlation matrix.

- [`check_normality`](https://because-pkg.github.io/because/reference/check_normality.md)
  /
  [`plot_qq`](https://because-pkg.github.io/because/reference/plot_qq.md):
  Bayesian Q-Q check with empirical posterior predictive simulation
  envelopes and skewness/kurtosis discrepancies.

- [`check_homoscedasticity`](https://because-pkg.github.io/because/reference/check_homoscedasticity.md):
  Binned posterior residual standard deviations across fitted values and
  variance ratio discrepancies.

- [`check_residual_structure`](https://because-pkg.github.io/because/reference/check_residual_structure.md):
  Posterior distribution of the Variance Partition Coefficient
  (\\\lambda\_{\text{Bayes}}\\) for phylogenetic or spatial random
  effects.

- [`check_model`](https://because-pkg.github.io/because/reference/check_model.md):
  Unified multi-panel diagnostic dashboard.

## References

de Villemereuil, P., & Nakagawa, S. (2014). General Quantitative Genetic
Methods for Comparative Biology. In: Garamszegi, L. (eds) *Modern
Phylogenetic Comparative Methods and Their Application in Evolutionary
Biology*. Springer, Berlin, Heidelberg

Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and graphical statistics*, 5(3), 236-244.

Leamer, E. E. (1973). Multicollinearity: A Bayesian interpretation.
*Review of Economics and Statistics*, 55(3), 371-380.

Leamer, E. E. (1978). *Specification Searches: Ad Hoc Inference with
Nonexperimental Data*. John Wiley & Sons.

McElreath, R. (2020). *Statistical Rethinking: A Bayesian Course with
Examples in R and Stan* (2nd ed.). CRC Press.

Gelman, A., Meng, X.-L., & Stern, H. (1996). Posterior predictive
assessment of model fitness via realized discrepancies. *Statistica
Sinica*, 6(4), 733-760.

Landwehr, J. M., Pregibon, D., & Shoemaker, A. C. (1984). Graphical
methods for assessing logistic regression models. *Journal of the
American Statistical Association*, 79(385), 61-71.

Hadfield, J. D. (2010). MCMC methods for multi-response generalized
linear mixed models: The MCMCglmm R package. *Journal of Statistical
Software*, 33(2), 1-22.
