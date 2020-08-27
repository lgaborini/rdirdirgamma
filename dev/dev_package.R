devtools::load_all('.')

library(magrittr)

# RNG ---------------------------------------------------------------------



seed <- sample.int(20, 1)

n <- 10
m <- 5
alpha_0 <- 1e1
beta_0 <- 1e1
p <- 3
nu_0 <- rep(1, p)

set.seed(seed)
mtx_beta <- rdirichlet_beta_cpp(1000, nu_0);
mtx_gamma <- purrr::map(seq(1000), ~ rdirichlet_cpp(nu_0, seed = sample.int(2000, 1))) %>% { do.call(rbind, .) }

# plot(mtx_beta[,1], mtx_beta[,2])
# points(mtx_gamma[,1], mtx_gamma[,2], col = 'red')



# Population --------------------------------------------------------------



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
