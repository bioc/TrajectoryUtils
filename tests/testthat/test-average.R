# Tests the averagePseudotime function.

test_that("averagePseudotime works", {
    pseudo <- matrix(rnorm(1000), ncol=4)
    pseudo[sample(length(pseudo), 200)] <- NA
    expect_equal(averagePseudotime(pseudo), rowMeans(pseudo, na.rm=TRUE))

    out <- PseudotimeOrdering(pseudo)
    expect_equal(averagePseudotime(out), rowMeans(pseudo, na.rm=TRUE))
})
