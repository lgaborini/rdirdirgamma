sample_ABC_rdirdirgamma_cpp_wrapper <- function(
   alpha_0, beta_0, nu_0,
   mtx_obs,
   method = 'euclidean', p = 2,
   reps = 1,
   n_sample, m_sample,
   summarize_ABC = c('mean', 'none')
) {

   stopifnot(!is.null(alpha_0))
   stopifnot(!is.null(beta_0))
   stopifnot(!is.null(nu_0))

   summarize_ABC <- match.arg(summarize_ABC)

   n_obs <- nrow(mtx_obs)
   p_var <- ncol(mtx_obs)

   # Precompute observed summary statistics
   mu_obs <- colMeans(mtx_obs)
   sd_obs <- apply(mtx_obs, sd, MARGIN = 1)

   # Generate one dataset, compute 2 summary statistics
}
