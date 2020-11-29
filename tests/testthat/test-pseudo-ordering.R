# This tests the capabilities of the PseudotimeOrdering class.
# library(testthat); library(TrajectoryUtils); source("test-pseudo-ordering.R")

ncells <- 200
npaths <- 5
orderings <- matrix(rnorm(1000), ncells, npaths)

test_that("PseudotimeOrdering constructors work as expected", {
    pto <- PseudotimeOrdering(orderings)
    expect_identical(pathStat(pto), orderings)

    pto <- PseudotimeOrdering(list(ordering=orderings))
    expect_identical(pathStatNames(pto), "ordering")

    pto <- PseudotimeOrdering(list(ordering=orderings), 
        pathData=DataFrame(row.names=LETTERS[1:npaths], X=1:5))
    expect_identical(rownames(pathData(pto)), LETTERS[1:npaths])
    expect_identical(colnames(pathStat(pto)), LETTERS[1:npaths])
 
    pto <- PseudotimeOrdering(list(ordering=orderings), 
        cellData=DataFrame(row.names=1:ncells))
    expect_identical(rownames(cellData(pto)), as.character(1:ncells))
    expect_identical(rownames(pathStat(pto)), as.character(1:ncells))
})

test_that("general dim utilities work", {
    pto <- PseudotimeOrdering(orderings)
    expect_equal(npaths(pto), npaths)
    expect_equal(ncells(pto), ncells)

    expect_null(cellnames(pto))
    replacements <- paste0("CELL_", seq_len(ncells))
    cellnames(pto) <- replacements
    expect_identical(cellnames(pto), replacements)

    expect_null(pathnames(pto))
    replacements <- paste0("CELL_", seq_len(npaths))
    pathnames(pto) <- replacements
    expect_identical(pathnames(pto), replacements)
})

test_that("pathStats getters and setters work", {
    pto <- PseudotimeOrdering(orderings)
    pathStat(pto, 2) <- orderings + 1
    expect_identical(pathStat(pto, 2), orderings+1)
    pathStat(pto, "YAY") <- orderings * 2
    expect_identical(pathStat(pto, "YAY"), orderings * 2)

    pathStatNames(pto) <- LETTERS[1:3]
    expect_identical(pathStatNames(pto), LETTERS[1:3])

    pathStats(pto) <- list(alternative=orderings)
    expect_identical(pathStats(pto), List(alternative=orderings))
    expect_identical(pathStatNames(pto), "alternative")
    expect_identical(pathStat(pto, 1), orderings)
})

test_that("annotation getters and setters work", {
    pto <- PseudotimeOrdering(orderings)
    pathData(pto)$X <- letters[1:npaths]
    expect_identical(pathData(pto)$X, letters[1:npaths])

    cellData(pto)$X <- 1:ncells
    expect_identical(cellData(pto)$X, 1:ncells)

    pto$X <- ncells:1
    expect_identical(cellData(pto)$X, ncells:1)
    expect_identical(pto$X, ncells:1)
})

