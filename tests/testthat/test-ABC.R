
# Summary statistics ------------------------------------------------------

n <- 3*8
p <- 3
mtx_obs <- matrix(rnorm(n), ncol = p)


test_that('get_summary_statistics_cpp (standard) works', {

   n_summary <- expect_silent(get_number_summary_statistics(FALSE))

   mtx_summary <- expect_silent(get_summary_statistics_cpp(mtx_obs, FALSE))
   expect_is(mtx_summary, 'matrix')
   expect_equal(nrow(mtx_summary), n_summary)
   expect_equal(ncol(mtx_summary), p)
})

test_that('get_summary_statistics_cpp (optimized) works', {

   n_summary <- expect_silent(get_number_summary_statistics(TRUE))

   mtx_summary <- expect_silent(get_summary_statistics_cpp(mtx_obs, TRUE))
   expect_is(mtx_summary, 'matrix')
   expect_equal(nrow(mtx_summary), n_summary)
   expect_equal(ncol(mtx_summary), p)
})

# Test ABC ----------------------------------------------------------------


test_that("sample_ABC_rdirdirgamma_beta_cpp works", {

   n_summary <- expect_silent(get_number_summary_statistics(FALSE))

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
   expect_equal(ncol(mtx_dist), n_summary)
   expect_true(all(mtx_dist >= 0))
})


test_that("sample_ABC_rdirdirgamma_beta_cpp (optimized) works", {

   n_summary <- expect_silent(get_number_summary_statistics(TRUE))

   mtx_dist <- expect_silent(sample_ABC_rdirdirgamma_beta_cpp(
         n_sample = 10,
         m_sample = 10,
         alpha_0 = 1,
         beta_0 = 1,
         nu_0 = rep(1, p),
         mtx_obs,
         reps = 100,
         p_norm = 2,
         use_optimized_summary = TRUE
   ))

   expect_is(mtx_dist, 'matrix')
   expect_equal(nrow(mtx_dist), 100)
   expect_equal(ncol(mtx_dist), n_summary)
   expect_true(all(mtx_dist >= 0 | is.nan(mtx_dist)))
})


test_that("sample_ABC_rdirdirgamma_beta_cpp works: not enough rows", {

   mtx_dist <- expect_error(sample_ABC_rdirdirgamma_beta_cpp(
         n_sample = 2,
         m_sample = 2,
         alpha_0 = 1,
         beta_0 = 1,
         nu_0 = rep(1, p),
         mtx_obs,
         reps = 100,
         p_norm = 2
   ), regexp = 'cannot generate enough observations')

})

test_that("sample_ABC_rdirdirgamma_beta_cpp works: wrong columns", {

   mtx_dist <- expect_error(sample_ABC_rdirdirgamma_beta_cpp(
         n_sample = 10,
         m_sample = 10,
         alpha_0 = 1,
         beta_0 = 1,
         nu_0 = rep(1, p + 1),
         mtx_obs,
         reps = 100,
         p_norm = 2
   ), regexp = 'different number of columns')

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
