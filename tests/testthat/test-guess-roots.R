# This tests that the root-guessing machinery works as expected.
# library(testthat); library(TrajectoryUtils); source('test-guess-roots.R')

library(igraph)
edges <- c("A", "B", "B", "C", "C", "D", "C", "E")
g <- make_graph(edges, directed=FALSE)

# Works with multiple components.
edges2 <- c(edges, "F", "G", "G", "H")
g2 <- make_graph(edges2, directed=FALSE)

empty <- make_graph(character(0))

test_that("guessMSTRoots works with method='degree1'", {
    expect_true(guessMSTRoots(g) %in% c("A", "D", "E"))
    
    multi <- guessMSTRoots(g2)
    expect_identical(length(multi), 2L)
    expect_true(multi[1] %in% c("A", "D", "E"))
    expect_true(multi[2] %in% c("F", "H"))

    expect_identical(guessMSTRoots(empty), character(0))
})

test_that("guessMSTRoots works with method='maxstep'", {
    expect_identical(guessMSTRoots(g, method="maxstep"), "A")
    expect_identical(guessMSTRoots(g2, method="maxstep"), c("A", "F"))

    # Mixing it up to make it use another endpoint.
    edges <- c("A", "B", "B", "C", "C", "D", "C", "E", "E", "F", "F", "G", "G", "H")
    g <- make_graph(edges, directed=FALSE)
    expect_identical(guessMSTRoots(g, method="maxstep"), "H")
})

test_that("guessMSTRoots works with method='maxlen'", {
    # Without weights. 
    expect_identical(guessMSTRoots(g, method="maxlen"), "A")
    expect_identical(guessMSTRoots(g2, method="maxlen"), c("A", "F"))

    # Adding some weights.
    E(g)$weight <- c(1,1,1,10)
    expect_identical(guessMSTRoots(g, method="maxlen"), "E")
})

test_that("guessMSTRoots works with method='minstep'", {
    expect_identical(guessMSTRoots(g, method="minstep"), "C")
    expect_identical(guessMSTRoots(g2, method="minstep"), c("C", "G"))
})

test_that("guessMSTRoots works with method='minlen'", {
    edges <- c("A", "C", "B", "C", "C", "D", "D", "E", "E", "F", "E", "G", "E", "H")
    g <- make_graph(edges, directed=FALSE)
    expect_identical(guessMSTRoots(g, method="minlen"), "E")

    # Throwing on some weights, but it doesn't really matter,
    # as the number of paths is the determining factor here.
    E(g)$weight <- c(1,1,0.1,10,1,1,1)
    expect_identical(guessMSTRoots(g, method="minlen"), "E")
})
