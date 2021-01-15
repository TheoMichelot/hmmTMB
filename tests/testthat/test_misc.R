
context("Miscellaneous")

test_that("make_formulas returns correct output", {
    # Define state-specific formulas
    f1 <- ~ state1(x1) + x2
    f2 <- ~ state1(x1) + state2(s(x2, bs = "cs"))
    input_forms <- list(step = list(shape = f1), count = list(lambda = f2))
    
    f <- make_formulas(input_forms = input_forms,
                       var_names = c("step", "count"),
                       par_names = list(step = c("shape", "scale"), count = "lambda"),
                       n_states = 2)
    
    # Check format
    expect_equal(length(f), length(input_forms))
    expect_equal(unname(sapply(f, length)), c(2, 1))
    
    # Check formulas in output
    f <- unlist(f)
    expect_equal(f[[1]], ~1 + x1 + x2)
    expect_equal(f[[2]], ~1 + x2)
    expect_equal(f[[3]], ~1)
    expect_equal(f[[4]], ~1)
    expect_equal(f[[5]], ~1 + x1)
    expect_equal(f[[6]], ~1 + s(x2, bs = "cs"))    
})
