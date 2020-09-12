# rdirdirgamma 0.2.1

* Add `compute_ABC_cpp` to directly compute the acceptance ratio.

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
