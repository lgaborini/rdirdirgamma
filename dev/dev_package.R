devtools::load_all('.')

n <- 10
m <- 5
alpha_0 <- 1e1
beta_0 <- 1e1
nu_0 <- rep(1, 5)/5

mtx_obs <- rdirdirgamma_cpp(
   n = n,
   m = m,
   alpha_0 = alpha_0,
   beta_0 = beta_0,
   nu_0 = nu_0
)


sample_ABC_rdirdirgamma_cpp_internal(
   n = n,
   m = m,
   alpha_0 = alpha_0,
   beta_0 = beta_0,
   nu_0 = nu_0,
   mtx_obs = mtx_obs,
   reps = 100,
   p_norm = 2
)
