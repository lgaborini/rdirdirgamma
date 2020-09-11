devtools::load_all('.')

library(magrittr)

# RNG ---------------------------------------------------------------------



seed <- sample.int(20, 1)

p <- 3
conc <- 20
nu_0 <- c(0.4, 0.6, 0.5)
nu_0 <- nu_0 / sum(nu_0) * conc

set.seed(seed)
mtx_beta <- rdirichlet_beta_cpp(1000, nu_0);
mtx_gamma <- purrr::map(seq(1000), ~ rdirichlet_cpp(nu_0, seed = sample.int(2000, 1))) %>% { do.call(rbind, .) }

plot(mtx_beta[,1], mtx_beta[,2], xlim = c(0, 1), ylim = c(0, 1))
points(mtx_gamma[,1], mtx_gamma[,2], col = 'red')


# Population --------------------------------------------------------------


n <- 100
m <- 5
alpha_0 <- 1e1
beta_0 <- 1e-1
conc <- 5
nu_0 <- rep(1, p)
nu_0 <- nu_0 / sum(nu_0) * conc

# set.seed(seed)

mtx_obs <- rdirdirgamma_cpp(
   n = n,
   m = m,
   alpha_0 = alpha_0,
   beta_0 = beta_0,
   nu_0 = nu_0
)

# set.seed(seed)

mtx_obs_beta <- rdirdirgamma_beta_cpp(
   n = n,
   m = m,
   alpha_0 = alpha_0,
   beta_0 = beta_0,
   nu_0 = nu_0
)

plot(mtx_obs[,1], mtx_obs[,2], xlim = c(0, 1), ylim = c(0, 1))
points(mtx_obs_beta[,1], mtx_obs_beta[,2], col = 'red')


# ABC / Rcpp --------------------------------------------------------------



mtx_ABC <- sample_ABC_rdirdirgamma_beta_cpp(
   n_sample = n,
   m_sample = m,
   alpha_0 = alpha_0,
   beta_0 = beta_0,
   nu_0 = nu_0,
   mtx_obs = mtx_obs,
   reps = 1000,
   p_norm = Inf
)

mtx_ABC

# Crashing
#
# sample_ABC_rdirdirgamma_cpp(
#    n_sample = n,
#    m_sample = m,
#    alpha_0 = alpha_0,
#    beta_0 = beta_0,
#    nu_0 = nu_0,
#    mtx_obs = mtx_obs,
#    reps = 100,
#    p_norm = 2
# )



# Summary statistics and norms -------------------------------------------------------------------


mtx_obs <- rdirdirgamma_beta_cpp(
   n = n,
   m = m,
   alpha_0 = alpha_0,
   beta_0 = beta_0,
   nu_0 = nu_0
)


mtx_gen <- rdirdirgamma_beta_cpp(
   n = n,
   m = m,
   alpha_0 = alpha_0,
   beta_0 = beta_0,
   nu_0 = nu_0
)



mtx_obs %>%
   rdirdirgamma::get_optimized_summary_statistics_cpp()


compute_distances_gen_obs_cpp(mtx_gen, mtx_obs, p_norm = 2, use_optimized_summary = FALSE)



# Optimized summary statistics ------------------------------------------------------


mtx_obs %>%
   rdirdirgamma::get_optimized_summary_statistics_cpp()


compute_distances_gen_obs_cpp(mtx_gen, mtx_obs, p_norm = 1, use_optimized_summary = FALSE)
compute_distances_gen_obs_cpp(mtx_gen, mtx_obs, p_norm = 1, use_optimized_summary = TRUE)
