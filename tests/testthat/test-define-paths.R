# This tests the defineMSTPaths function.
# library(testthat); library(TrajectoryUtils); source("test-define-paths.R")

library(igraph)
g <- make_graph(c("A", "B", "B", "C", "B", "D"), directed=FALSE)

test_that("defineMSTPaths works with preset roots", {
    paths <- defineMSTPaths(g, roots="A")
    expect_identical(paths, list(c("A", "B", "C"), c("A", "B", "D")))

    paths <- defineMSTPaths(g, roots="B")
    expect_identical(paths, list(c("B", "A"), c("B", "C"), c("B", "D")))

    # Testing for multiple components:
    g2 <- make_graph(c("A", "B", "B", "C", "B", "D", "E", "F"), directed=FALSE)
    expect_error(defineMSTPaths(g2, roots="B"), "one node")

    paths <- defineMSTPaths(g2, roots=c("D", "E"))
    expect_identical(paths, list(c("D", "B", "A"), c("D", "B", "C"), c("E", "F")))
})

test_that("defineMSTPaths works with per-node times", {
    paths <- defineMSTPaths(g, times=c(A=0, B=1, C=2, D=3))
    expect_identical(paths, list(c("A", "B", "C"), c("A", "B", "D")))

    paths <- defineMSTPaths(g, times=c(A=0, B=-1, C=2, D=3))
    expect_identical(paths, list(c("B", "A"), c("B", "C"), c("B", "D")))

    # Multiple convergent paths generated:
    paths <- defineMSTPaths(g, times=c(A=0, B=5, C=2, D=3))
    expect_identical(paths, list(c("A", "B"), c("C", "B"), c("D", "B")))

    paths <- defineMSTPaths(g, times=c(A=0, B=5, C=2, D=10))
    expect_identical(paths, list(c("A", "B", "D"), c("C", "B", "D")))

    # Works for multiple components.
    g2 <- make_graph(c("A", "B", "B", "C", "B", "D", "E", "F"), directed=FALSE)
    paths <- defineMSTPaths(g2, times=c(A=0, B=1, C=2, D=3, E=10, F=0))
    expect_identical(paths, list(c("A", "B", "C"), c("A", "B", "D"), c("F", "E")))
})

set.seed(1000)
test_that("per-node timings can be computed from clusters", {
    clusters <- sample(LETTERS[1:4], 100, replace=TRUE)
    timings <- runif(100)

    ref <- defineMSTPaths(g, times=vapply(split(timings, clusters), mean, 0))
    out <- defineMSTPaths(g, times=timings, cluster=clusters)
    expect_identical(ref, out)

    ref <- defineMSTPaths(g, times=vapply(split(timings, clusters), median, 0))
    out <- defineMSTPaths(g, times=timings, cluster=clusters, use.median=TRUE)
    expect_identical(ref, out)
})
