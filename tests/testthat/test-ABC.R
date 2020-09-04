test_that("sample_ABC_rdirdirgamma_beta_cpp works", {

   n <- 3*8
   p <- 3
   mtx_obs <- matrix(rnorm(n), ncol = p)

   mtx_dist <- expect_silent(sample_ABC_rdirdirgamma_beta_cpp(
         n_sample = 10,
         m_sample = 10,
         alpha_0 = 1,
         beta_0 = 1,
         nu_0 = rep(1, p),
         mtx_obs,
         reps = 100,
         p_norm = 2
   ))

   expect_is(mtx_dist, 'matrix')
   expect_equal(nrow(mtx_dist), 100)
   expect_equal(ncol(mtx_dist), 2)
   expect_true(all(mtx_dist >= 0))
})

# test_that("sample_ABC_rdirdirgamma_cpp works", {
#
#    n <- 3*8
#    p <- 3
#    mtx_obs <- matrix(rnorm(n), ncol = p)
#
#    mtx_dist <- sample_ABC_rdirdirgamma_cpp(
#          n_sample = 10,
#          m_sample = 10,
#          alpha_0 = 1,
#          beta_0 = 1,
#          nu_0 = rep(1, p),
#          mtx_obs,
#          reps = 100,
#          p_norm = 2
#    )
#
#    expect_is(mtx_dist, 'matrix')
#    expect_equal(nrow(mtx_dist), 100)
#    expect_equal(ncol(mtx_dist), 2)
#    expect_true(all(mtx_dist >= 0))
# })
