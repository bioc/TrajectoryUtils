# Various utilities work as expected"
# library(testthat); library(TrajectoryUtils); source("test-utils.R")

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
