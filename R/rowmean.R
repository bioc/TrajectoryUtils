#' Compute column means based on a grouping variable
#'
#' Computes the column mean or median for each group of rows in a matrix.
#'
#' @param x A numeric matrix or matrix-like object.
#' @param group A vector or factor specifying the group assignment for each row of \code{x}.
#' @param ... Further arguments to pass to \code{\link{rowMeans}} or \code{\link{rowMedians}}.
#'
#' @details
#' The naming scheme here is somewhat inspired by the \code{rowsum} function from the \pkg{DelayedArray} package,
#' which is in turn a transposed counterpart to \code{\link{colsum}}.
#' Admittedly, it is rather confusing when \code{\link{rowMeans}} computes the mean for a row across all columns
#' while \code{rowmean} computes the mean for a column across a subset of rows, but there you have it.
#'
#' @return A numeric matrix with one row per level of \code{group},
#' where the value for each column contains the mean or median across the subset of rows corresponding that level.
#'
#' @author Aaron Lun
#' 
#' @examples
#' x <- matrix(runif(100), ncol = 5)
#' group <- sample(1:8, 20, TRUE)
#' (xmean <- rowmean(x, group))
#' (xmeds <- rowmedian(x, group))
#'
#' @export
#' @importFrom Matrix colMeans
rowmean <- function(x, group, ...) {
    .rowstats(x, group, FUN=colMeans)
}

#' @export
#' @rdname rowmean
#' @importFrom MatrixGenerics colMedians
rowmedian <- function(x, group, ...) {
    .rowstats(x, group, FUN=colMedians, ...)
}

.rowstats <- function(x, group, FUN, ...) {
    by.group <- split(seq_len(nrow(x)), group)
    output <- matrix(0, length(by.group), ncol(x), dimnames=list(names(by.group), colnames(x)))    
    for (i in seq_along(by.group)) {
        output[i,] <- FUN(x[by.group[[i]],,drop=FALSE], ...)
    }
    output
}
