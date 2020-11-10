#' Guess the roots of a MST
#'
#' Pick nodes to use as the root(s) of an MST using a variety of ad hoc methods.
#'
#' @param g A \link{graph} object containing a MST.
#' All nodes should be named.
#' @param method String specifying the method to use to pick the root.
#'
#' @return A character vector containing the identity of the root for each component in \code{g}.
#'
#' @details
#' When \code{method="degree1"}, an arbitrary node of degree 1 is chosen as the root for each component.
#' This aims to reduce the number of branch events by starting from an already-terminal node.
#'
#' When \code{method="maxstep"}, we pick the node of degree 1 that has the highest average number of steps to reach all other nodes of degree 1.
#' This aims to maximize the number of shared clusters between different paths when traversing \code{g} from the root to the other terminal nodes,
#' under the philosophy that branch events should occur as late as possible.
#' When \code{method="maxlen"}, we instead pick the node of degree 1 that has the highest average distance to reach all other nodes of degree 1.
#' This also considers the distance spanned by each cluster.
#'
#' When \code{method="minstep"}, we pick the node that has the lowest average number of steps to reach all nodes of degree 1.
#' This aims to minimize the number of shared clusters between different paths under the philosophy that branch events should occur as early as possible.
#' When \code{method="minlen"}, we instead pick the node that has the highest average distance to reach all nodes of degree 1.
#'
#' @author Aaron Lun,
#' based on code by Kelly Street
#'
#' @examples
#' library(igraph)
#' edges <- c("A", "B", "B", "C", "C", "D", "C", "E")
#' g <- make_graph(edges, directed=FALSE)
#'
#' guessMSTRoots(g)
#' guessMSTRoots(g, method="maxstep")
#' guessMSTRoots(g, method="minstep")
#'
#' # Works with multiple components.
#' edges2 <- c(edges, "F", "G", "G", "H")
#' g2 <- make_graph(edges2, directed=FALSE)
#' guessMSTRoots(g2)
#'
#' @export
#' @importFrom igraph shortest_paths distances decompose degree V
guessMSTRoots <- function(g, method=c("degree1", "maxstep", "maxlen", "minstep", "minlen")) {
    method <- match.arg(method)

    forest <- decompose(g)
    roots <- vector("list", length(forest))

    for (i in seq_along(forest)) {
        tree <- forest[[i]]
        deg <- degree(tree)
        named <- names(V(tree))
        deg1 <- named[deg==1L]

        if (method %in% c("degree1", "maxstep", "maxlen")) {
            candidates <- deg1
        } else {
            candidates <- named
        }

        if (method=="degree1") {
            roots[[i]] <- candidates[1]
        } else {
            averages <- numeric(length(candidates))
            names(averages) <- candidates
            for (j in candidates) {
                ends <- setdiff(deg1, j)

                if (method %in% c("minstep", "maxstep")) {
                    paths <- shortest_paths(g, from = j, to = ends, output = 'vpath')$vpath
                    averages[j] <- mean(lengths(paths))
                } else {
                    dists <- distances(g, v = j, to = ends)
                    averages[j] <- mean(dists)
                }
            }

            if (method %in% c("maxstep", "maxlen")) {
                FUN <- which.max
            } else {
                FUN <- which.min
            }
            roots[[i]] <- candidates[FUN(averages)]
        }
    }

    unlist(roots)
}
