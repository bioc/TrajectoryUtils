# This tests the splitByBranches function.
# library(testthat); library(TrajectoryUtils); source("test-split-branches.R")

library(igraph)

test_that("splitByBranches works correctly", {
    G <- make_graph(c("A", "B", "B", "C", "C", "D", "D", "E", "D", "F", "F", "G"), directed=FALSE)
    out <- splitByBranches(G)
    expect_identical(out, list(LETTERS[4:1], LETTERS[5:4], c("G", "F", "D"))) 

    # Works with multiple components.
    G <- make_graph(c("A", "B", "B", "C", "C", "D", "D", "E", "F", "G"), directed=FALSE)
    out <- splitByBranches(G)
    expect_identical(out, list(LETTERS[5:1], LETTERS[7:6]))

    # Works with no-edge components.
    G <- make_graph(c("A", "B", "B", "C", "C", "D", "D", "E", "F", "F"), directed=FALSE)
    out <- splitByBranches(simplify(G))
    expect_identical(out, list(LETTERS[5:1], "F"))

    # Works well enough with an empty graph.
    expect_identical(splitByBranches(make_graph(character(0))), list())
})
