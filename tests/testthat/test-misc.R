
# Colwise functions -------------------------------------------------------


test_that("column-wise sd", {

   mtx <- matrix(rnorm(12), nrow = 4, ncol = 3)

   mtx_sd_target <- apply(mtx, 2, sd)

   mtx_sd_target_cpp <- expect_silent(colsd(mtx))

   expect_equal(mtx_sd_target, mtx_sd_target_cpp)
})

test_that("column-wise kurtosis", {

   mtx <- matrix(rnorm(12), nrow = 4, ncol = 3)

   mtx_k_target <- matrix(moments::kurtosis(mtx), nrow = 1)
   mtx_k_target_cpp <- expect_silent(colkurtosis(mtx))

   expect_equal(mtx_k_target, mtx_k_target_cpp)
})

test_that("column-wise skewness", {

   mtx <- matrix(rnorm(12), nrow = 4, ncol = 3)

   mtx_s_target <- matrix(moments::skewness(mtx), nrow = 1)
   mtx_s_target_cpp <- expect_silent(colskewness(mtx))

   expect_equal(mtx_s_target, mtx_s_target_cpp)
})


# Norms -------------------------------------------------------------------

dist_minkowski <- function(x, p_norm) {
   p <- length(x)
   dist(rbind(x, rep(0, p)), p = p_norm, method = 'minkowski')[1]
}

test_that('norm_minkowski is correct', {

   v <- c(1, 2, 0)

   expect_equal(norm_minkowski(v, 1), dist_minkowski(v, 1))
   expect_equal(norm_minkowski(v, 2), dist_minkowski(v, 2))
   expect_equal(norm_minkowski(v, Inf), max(abs(v)))
   expect_error(norm_minkowski(v, 0))

   v <- c(-1, 2, 0)

   expect_equal(norm_minkowski(v, 1), dist_minkowski(v, 1))
   expect_equal(norm_minkowski(v, 2), dist_minkowski(v, 2))
   expect_equal(norm_minkowski(v, Inf), max(abs(v)))
   expect_error(norm_minkowski(v, 0))

})




# ABC distances -----------------------------------------------------------

test_that('compute_distances_gen_obs_cpp is correct', {

   set.seed(123)
   mtx_gen <- matrix(rnorm(12), nrow = 4, ncol = 3)
   set.seed(124)
   mtx_obs <- matrix(rnorm(12), nrow = 4, ncol = 3)

   v <- expect_silent(compute_distances_gen_obs_cpp(mtx_gen, mtx_gen))
   expect_true(all.equal(v, rep(0, 2)))

   v <- expect_silent(compute_distances_gen_obs_cpp(mtx_gen, mtx_obs))
   expect_length(v, 2)
   expect_true(all(v > 0))
})
