#' Compute the average pseudotime
#'
#' Compute the average pseudotime for each cell across all paths in which it is involved.
#' 
#' @param x A numeric matrix-like object containing pseudotime orderings.
#' Alternatively, a \linkS4class{PseudotimeOrdering} object containing such a matrix.
#' @param i Integer scalar or string specifying the entry of \code{\link{pathStats}(x)} containing the pseudotime matrix,
#' if \code{x} is a PseudotimeOrdering object.
#'
#' @return A numeric vector containing the average pseudotime for each cell.
#'
#' @details
#' Averaging the pseudotime is a convenient way to consolidate multiple paths into a single ordering for, e.g., visualization.
#' It is permissible as cells involved in multiple paths should generally have similar pseudotimes in each path,
#' under assumption that the cell is involved in the part of the trajectory that is common to those paths.
#' In such cases, the average is just a way of compressing those pseudotimes into a single value.
#' Conversely, for cells that are unique to a single path, the average collapses to that path's pseudotime (assuming that all other values are \code{NA}).
#' 
#' @author Aaron Lun
#'
#' @examples
#' pseudotimes <- matrix(rnorm(200), ncol=2)
#' pseudotimes[1:40,1] <- NA
#' pseudotimes[61:100,2] <- NA
#' pseudotimes[41:60,] <- runif(20)
#' averagePseudotime(pseudotimes)
#'
#' pto <- PseudotimeOrdering(pseudotimes)
#' averagePseudotime(pto)
#' 
#' @export
#' @importFrom Matrix rowMeans
averagePseudotime <- function(x, i=1L) {
    if (is(x, "PseudotimeOrdering")) {
        x <- pathStat(x, i)
    }
    rowMeans(x, na.rm=TRUE)
}
