# rdirdirgamma 0.2.5 (unreleased)

* Removed unused function.
* Add Remotes in DESCRIPTION.
* Updated `renv` to 0.14.0.
* Updated `renv` lockfile.

# rdirdirgamma 0.2.4

* Add URL to Description.
* Add pkgdown documentation to Rbuildignore.

# rdirdirgamma 0.2.3

* Now summary statistics returned by `compute_ABC_cpp` are documented in the correct order.

# rdirdirgamma 0.2.2

* Interface change: number of ABC samples is always `reps`.
* Add pkgdown, better documentation.

# rdirdirgamma 0.2.1

* Add `compute_ABC_cpp` to directly compute the acceptance ratio and return the sampled distances.

# rdirdirgamma 0.2

* Adding more C++ functions to compute distances between two datasets.
* Gains RcppArmadillo dependency.
* Adding functions:
   - `colsd`
   - `colskewness`
   - `colkurtosis`
   - `norm_minkowski`
   - `compute_distances_gen_obs_cpp`
* Drafting tests for ABC, misc functions.
* Now ABC generates datasets with the same number of rows as the observed dataset.   
  Sampling and distance calculation functions are unchanged, however.
* Documented RNG seed for R and GSL.


* Removed `sample_ABC_rdirdirgamma_cpp`.
* Split code to cpp/headers.

* Add generation of samples that would be accepted by ABC.

* Summary statistics are now generic.

  Two sets:
  - standard: mean, sd
  - optimized: mean, sd, skewness, kurtosis (in development)

# rdirdirgamma 0.1.0

* Added a `NEWS.md` file to track changes to the package.
