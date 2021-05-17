
context("Dist class")

test_that("Distribution can be accessed and correct information is returned", {
    obj <- dist_norm
    expect_equal(obj$name(), "norm")
    expect_equal(obj$pdf(), dnorm)    
    expect_equal(obj$rng(), rnorm)
    expect_equal(obj$link(), list(mean = identity, sd = log))
    expect_equal(obj$invlink(), list(mean = identity, sd = exp))
    expect_equal(obj$npar(), 2)
    expect_equal(obj$parnames(), c("mean", "sd"))
    expect_equal(obj$parapprox(), function(x) {return(c(mean(x), sd(x)))})
    expect_equal(obj$fixed(), c(FALSE, FALSE))
    expect_equal(obj$code(), 1)
})

test_that("n2w and w2n return right format (dist_norm)", {
    par <- list(mean = c(-5, 5), sd = c(2, 10))
    wpar <- dist_norm$n2w(par = par)
    expect_is(wpar, "numeric")
    expect_equal(length(wpar), length(unlist(par)))

    par2 <- dist_norm$w2n(wpar = wpar)
    expect_is(par2, "list")
    expect_equal(length(par2), length(par))
    
    par3 <- dist_norm$w2n(wpar = wpar, as_matrix = TRUE)
    expect_is(par3, "matrix")
    expect_equal(dim(par3), c(length(par$mean), length(par)))
    expect_equal(colnames(par3), names(par))
})

test_that("pdf_apply returns the right output (dist_norm)", {
    x <- pi
    par <- c(mean = 1, sd = 5)
    expect_equal(dist_norm$pdf_apply(x = x, par = par, log = FALSE),
                 dnorm(x, mean = par[1], sd = par[2], log = FALSE))
    expect_equal(dist_norm$pdf_apply(x = x, par = par, log = TRUE),
                 dnorm(x, mean = par[1], sd = par[2], log = TRUE))
})
