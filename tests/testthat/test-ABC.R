
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



# Full ABC computation ----------------------------------------------------



test_that("compute_ABC_cpp works: standard, no distances", {

   n_summary <- expect_silent(get_number_summary_statistics(FALSE))
   vec_eps <- rep(1, n_summary)

   list_ABC <- expect_silent(compute_ABC_cpp(
         n_sample = 10,
         m_sample = 10,
         alpha_0 = 1,
         beta_0 = 1,
         nu_0 = rep(1, p),
         mtx_obs,
         summarize_eps = vec_eps,
         reps = 100,
         p_norm = 2,
         use_optimized_summary = FALSE,
         return_distances = FALSE
   ))

   expect_is(list_ABC, 'list')
   expect_named(list_ABC, c('n_accepted', 'accept_ratio'))
   expect_gte(list_ABC$n_accepted, 0)
   expect_lte(list_ABC$n_accepted, 100)
   expect_gte(list_ABC$accept_ratio, 0)
   expect_lte(list_ABC$accept_ratio, 1)
})


test_that("compute_ABC_cpp works: optimized, no distances", {

   n_summary <- expect_silent(get_number_summary_statistics(TRUE))
   vec_eps <- rep(1, n_summary)

   list_ABC <- expect_silent(compute_ABC_cpp(
         n_sample = 10,
         m_sample = 10,
         alpha_0 = 1,
         beta_0 = 1,
         nu_0 = rep(1, p),
         mtx_obs,
         summarize_eps = vec_eps,
         reps = 100,
         p_norm = 2,
         use_optimized_summary = TRUE,
         return_distances = FALSE
   ))

   expect_is(list_ABC, 'list')
   expect_named(list_ABC, c('n_accepted', 'accept_ratio'))
   expect_gte(list_ABC$n_accepted, 0)
   expect_lte(list_ABC$n_accepted, 100)
   expect_gte(list_ABC$accept_ratio, 0)
   expect_lte(list_ABC$accept_ratio, 1)
})

test_that("compute_ABC_cpp works: standard, with distances", {

   n_summary <- expect_silent(get_number_summary_statistics(TRUE))
   vec_eps <- rep(1, n_summary)

   list_ABC <- expect_silent(compute_ABC_cpp(
         n_sample = 10,
         m_sample = 10,
         alpha_0 = 1,
         beta_0 = 1,
         nu_0 = rep(1, p),
         mtx_obs,
         summarize_eps = vec_eps,
         reps = 100,
         p_norm = 2,
         use_optimized_summary = TRUE,
         return_distances = TRUE
   ))

   expect_is(list_ABC, 'list')
   expect_named(list_ABC, c('n_accepted', 'accept_ratio', 'd_ABC'))
   expect_gte(list_ABC$n_accepted, 0)
   expect_lte(list_ABC$n_accepted, 100)
   expect_gte(list_ABC$accept_ratio, 0)
   expect_lte(list_ABC$accept_ratio, 1)
   expect_is(list_ABC$d_ABC, 'matrix')
   expect_equal(nrow(list_ABC$d_ABC), 100)
   expect_equal(ncol(list_ABC$d_ABC), n_summary)
   expect_true(all(list_ABC$d_ABC >= 0  | is.nan(list_ABC$d_ABC)))
})


test_that("compute_ABC_cpp works: optimized, with distances", {

   n_summary <- expect_silent(get_number_summary_statistics(TRUE))
   vec_eps <- rep(1, n_summary)

   list_ABC <- expect_silent(compute_ABC_cpp(
         n_sample = 10,
         m_sample = 10,
         alpha_0 = 1,
         beta_0 = 1,
         nu_0 = rep(1, p),
         mtx_obs,
         summarize_eps = vec_eps,
         reps = 100,
         p_norm = 2,
         use_optimized_summary = TRUE,
         return_distances = TRUE
   ))

   expect_is(list_ABC, 'list')
   expect_named(list_ABC, c('n_accepted', 'accept_ratio', 'd_ABC'))
   expect_gte(list_ABC$n_accepted, 0)
   expect_lte(list_ABC$n_accepted, 100)
   expect_gte(list_ABC$accept_ratio, 0)
   expect_lte(list_ABC$accept_ratio, 1)
   expect_is(list_ABC$d_ABC, 'matrix')
   expect_equal(nrow(list_ABC$d_ABC), 100)
   expect_equal(ncol(list_ABC$d_ABC), n_summary)
   expect_true(all(list_ABC$d_ABC >= 0 | is.nan(list_ABC$d_ABC)))
})
