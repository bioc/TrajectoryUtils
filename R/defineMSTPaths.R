#' Define paths through the MST
#'
#' Define paths through the MST, either from pre-specified root nodes or based on external timing information.
#'
#' @param g A \link{graph} object containing a minimum spanning tree, e.g., from \code{\link{createClusterMST}}.
#' @param roots A character vector specifying the root node for each component in \code{g}.
#' @param times A numeric vector of length equal to the number of nodes in \code{g},
#' specifying the external time associated with each node.
#' This should be named with the name of each node.
#' Alternatively, a numeric vector of length equal to the number of cells, in which case \code{clusters} must be specified.
#' @param clusters A vector or factor specifying the assigned cluster for each cell,
#' where each cluster corresponds to a node in \code{g}.
#'
#' Alternatively, a matrix with number of rows equal to \code{nrow(x)}, 
#' containing soft assignment weights for each cluster (column).
#' All weights should be positive and sum to 1 for each row.
#' 
#' This only has an effect if \code{times} is set to a vector of length equal to the number of cells.
#' @param use.median Logical scalar indicating whether the time for each cluster is defined as the median time across its cells.
#' The mean is used by default.
#' This only has an effect if \code{clusters} is specified.
#'
#' @return
#' A list of character vectors.
#' Each vector contains the names of nodes in \code{g} and defines a path through the MST from a root to an endpoint node.
#' 
#' @author Aaron Lun
#'
#' @details
#' When \code{roots} is specified, a path is defined from the root to each endpoint node (i.e., with degree 1) in \code{g}.
#' We expect one root node to be specified for each component in \code{g}.
#'
#' When \code{times} is specified, a path is defined from each local minima in time to the nearest local maxima within each component of \code{g}.
#' Timing information can be defined from experimental metadata or with computational methods like RNA velocity.
#' 
#' @examples
#' library(igraph)
#' test.g <- make_graph(c("A", "B", "B", "C", "B", "D"), directed=FALSE)
#' defineMSTPaths(test.g, roots="A")
#' defineMSTPaths(test.g, roots="B")
#' 
#' defineMSTPaths(test.g, times=c(A=0, B=1, C=2, D=3))
#' defineMSTPaths(test.g, times=c(A=0, B=-1, C=2, D=3))
#' defineMSTPaths(test.g, times=c(A=0, B=5, C=2, D=3))
#'
#' @seealso
#' \code{\link{guessMSTRoots}}, to obtain \code{roots} without any prior information.
#'
#' \code{\link{splitByBranches}}, for a root-free way of obtaining paths.
#' 
#' @export
#' @importFrom igraph V degree graph_from_adjacency_matrix decompose shortest_paths
defineMSTPaths <- function(g, roots, times=NULL, clusters=NULL, use.median=FALSE) {
    forest <- decompose(g)
    output <- vector("list", length(forest))

    if (!missing(roots)) {
        for (i in seq_along(forest)) {
            tree <- forest[[i]]
            named <- names(V(tree))
            if (length(named)==1L) {
                output[[i]] <- list(named)
                next
            }

            is.root <- roots %in% named
            if (sum(is.root)!=1) {
                stop("'roots' should contain one node per component of 'g'")
            }
            cur.root <- roots[is.root]

            deg <- degree(tree)
            deg1 <- named[deg == 1L]
            paths <- shortest_paths(tree, from=cur.root, to=setdiff(deg1, cur.root))$vpath
            output[[i]] <- lapply(paths, names)
        }
    } else {
        if (is.null(times)) {
            stop("'times' must be specified if 'roots' is missing")
        }
        if (!is.null(clusters)) {
            FUN <- if (use.median) rowmedian else rowmean
            times <- FUN(cbind(times), clusters)[,1]
        }

        for (i in seq_along(forest)) {
            tree <- forest[[i]]
            mat <- tree[]
            keep <- outer(times[rownames(mat)], times[colnames(mat)], "<")
            mat <- mat * keep
            g2 <- graph_from_adjacency_matrix(mat, mode="directed")

            all.nodes <- names(V(g2))
            all.starts <- all.nodes[degree(g2, mode="in") == 0]
            all.ends <- all.nodes[degree(g2, mode="out") == 0]

            collected <- vector("list", length(all.starts))
            for (s in seq_along(collected)) {
                # If it's not a neighboring local maxima, it won't be reachable,
                # because the directions of the edges won't allow it.
                reachable <- subcomponent(g2, all.starts[s], mode="out")
                cur.ends <- intersect(all.nodes[reachable], all.ends)
                paths <- shortest_paths(g2, from=all.starts[s], to=cur.ends, output="vpath")$vpath
                collected[[s]] <- lapply(paths, names)
            }

            output[[i]] <- unlist(collected, recursive=FALSE)
        }
    }

    if (length(output)) {
        unlist(output, recursive=FALSE)
    } else {
        output
    }
}
