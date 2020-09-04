

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
