
# Dirichlet ---------------------------------------------------------------



test_that("rdirichlet_cpp works", {

   r <- expect_silent(rdirichlet_cpp(rep(1, 3)))
   expect_length(r, 3)
   expect_type(r, 'double')
   expect_is(r, 'numeric')

   # RNG works without seed
   r1 <- rdirichlet_cpp(rep(1, 3))
   Sys.sleep(0.001)
   r2 <- rdirichlet_cpp(rep(1, 3))
   expect_true(!identical(r1, r2))

   # RNG works with seed
   r1 <- rdirichlet_cpp(rep(1, 3), seed = 123)
   r2 <- rdirichlet_cpp(rep(1, 3), seed = 123)
   expect_identical(r1, r2)
})



test_that("rdirichlet_beta_cpp works", {

   r <- expect_silent(rdirichlet_beta_cpp(n = 1, rep(1, 3)))
   expect_length(r, 3)
   expect_type(r, 'double')
   expect_is(r, 'matrix')
   expect_equal(nrow(r), 1)
   expect_equal(ncol(r), 3)

   r <- expect_silent(rdirichlet_beta_cpp(n = 10, rep(1, 3)))
   expect_type(r, 'double')
   expect_is(r, 'matrix')
   expect_equal(nrow(r), 10)
   expect_equal(ncol(r), 3)

})







# Dirichlet-Dirichlet-Gamma -----------------------------------------------


n <- 100
m <- 5
p <- 3
alpha_0 <- 1e1
beta_0 <- 1e-1
conc <- 5
nu_0 <- rep(1, p)
nu_0 <- nu_0 / sum(nu_0) * conc


test_that('rdirdirgamma_cpp works', {


   mtx_obs <- expect_silent(rdirdirgamma_cpp(
      n = n,
      m = m,
      alpha_0 = alpha_0,
      beta_0 = beta_0,
      nu_0 = nu_0
   ))

   expect_type(mtx_obs, 'double')
   expect_is(mtx_obs, 'matrix')
   expect_equal(nrow(mtx_obs), n*m)
   expect_equal(ncol(mtx_obs), p)
   expect_true(all(!is.na(mtx_obs)))
   expect_true(all(!is.nan(mtx_obs)))

})


test_that('rdirdirgamma_cpp works', {

   mtx_obs <- expect_silent(rdirdirgamma_cpp(
      n = n,
      m = m,
      alpha_0 = alpha_0,
      beta_0 = beta_0,
      nu_0 = nu_0
   ))

   expect_type(mtx_obs, 'double')
   expect_is(mtx_obs, 'matrix')
   expect_equal(nrow(mtx_obs), n*m)
   expect_equal(ncol(mtx_obs), p)
   expect_true(all(!is.na(mtx_obs)))
   expect_true(all(!is.nan(mtx_obs)))
})

