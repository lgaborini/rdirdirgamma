

test_that("column-wise sd", {

   mtx <- matrix(rnorm(12), nrow = 4, ncol = 3)

   mtx_sd_target <- apply(mtx, 2, sd)

   mtx_sd_target_cpp <- expect_silent(colsd(mtx))

   expect_equal(mtx_sd_target, mtx_sd_target_cpp)
})


# Norms -------------------------------------------------------------------

test_that('norm_minkowski is correct', {

   v <- c(1, 2, 0)

   expect_equal(norm_minkowski(v, 1), 3)
   expect_equal(norm_minkowski(v, 2), 5)
   expect_equal(norm_minkowski(v, Inf), 2)
   expect_equal(norm_minkowski(v, 0), 3)

})




# ABC distances -----------------------------------------------------------

test_that('compute_distances_gen_obs_cpp is correct', {

   set.seed(123)
   mtx_gen <- matrix(rnorm(12), nrow = 4, ncol = 3)
   set.seed(124)
   mtx_obs <- matrix(rnorm(12), nrow = 4, ncol = 3)

   v <- expect_silent(compute_distances_gen_obs_cpp(mtx_gen, mtx_gen))
   expect_true(all.equal(v, rep(0, 3)))

   v <- expect_silent(compute_distances_gen_obs_cpp(mtx_gen, mtx_obs))
   expect_length(v, 3)
   expect_true(all(v > 0))
})