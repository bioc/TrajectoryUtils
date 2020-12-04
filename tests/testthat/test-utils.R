# Various utilities work as expected.
# library(testthat); library(TrajectoryUtils); source("setup.R"); source("test-utils.R")

set.seed(100)
x <- matrix(runif(1000), ncol = 4)
group <- sample(1:10, 250, TRUE)
y <- Matrix::rsparsematrix(nrow(x), ncol(x), 0.1)

test_that('rowmean works as expected', {
    out <- rowmean(x, group)
    expect_identical(out[1,], colMeans(x[group==1,])) 
    expect_identical(out[10,], colMeans(x[group==10,])) 

    # Works for sparse matrices.
    out <- rowmedian(y, group)
    expect_identical(out[1,], colMedians(y[group==1,])) 
    expect_identical(out[10,], colMedians(y[group==10,])) 
})

test_that('weighted rowmean works as expected', {
    ref <- rowmean(x, group)
    sgroups <- factor2matrix(group)
    alt <- rowmean(x, sgroups)
    expect_equal(alt, ref)

    # Mean-centering is done correctly.
    mat <- matrix(runif(10 * nrow(x)), nrow(x))
    colnames(mat) <- letters[1:10]

    out <- rowmean(x, mat)
    expect_identical(rownames(out), colnames(mat))
    out2 <- rowmean(x, mat/rowMeans(mat))
    expect_equal(out, out2)
})

library(DelayedMatrixStats)
test_that('rowmedian works as expected', {
    out <- rowmedian(x, group)
    expect_identical(out[1,], colMedians(x[group==1,])) 
    expect_identical(out[10,], colMedians(x[group==10,])) 

    # Works for sparse matrices.
    out <- rowmedian(y, group)
    expect_identical(out[1,], colMedians(y[group==1,])) 
    expect_identical(out[10,], colMedians(y[group==10,])) 
})

test_that('weighted rowmedian works as expected', {
    ref <- rowmedian(x, group)
    sgroups <- factor2matrix(group)
    alt <- rowmedian(x, sgroups)
    expect_equal(alt, ref)

    # Mean-centering is done correctly.
    mat <- matrix(runif(10 * nrow(x)), nrow(x))
    colnames(mat) <- letters[1:10]

    out <- rowmedian(x, mat)
    expect_identical(rownames(out), colnames(mat))
    out2 <- rowmedian(x, mat/rowMeans(mat))
    expect_equal(out, out2)
})
