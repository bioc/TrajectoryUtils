#' Split a graph into branch-free paths
#'
#' Split a graph (usually a MST, at least a DAG) into branch-free components, i.e., paths between branch points or terminal nodes.
#'
#' @param g A \link{graph} object such as that produced by \code{\link{createClusterMST}}.
#'
#' @return A list of character vectors containins unbranched paths between branch points or terminal nodes.
#'
#' @details
#' This function implements another strategy to define paths through the MST,
#' by simply considering all linear sections (i.e., connected by nodes of degree 2) between branch points or terminii.
#' The idea is to enable characterisation of the continuum without external information or assumptions statements about where the root is positioned,
#' which would otherwise be required in functions such as \code{\link{defineMSTPaths}}.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(igraph)
#' test.g <- make_graph(c("A", "B", "B", "C", "C", "D", 
#'     "D", "E", "D", "F", "F", "G"), directed=FALSE)
#' splitByBranches(test.g)
#'
#' @seealso
#' \code{\link{defineMSTPaths}}, for the root-based method of defining paths.
#'
#' @export
#' @importFrom igraph subcomponent shortest_paths V degree
splitByBranches <- function(g) {
    all.nodes <- names(V(g))
    deg <- degree(g) 
    of.interest <- all.nodes[deg != 2 & deg != 0]

    collected <- vector("list", length(of.interest))
    for (i in seq_along(of.interest)) {
        reachable <- subcomponent(g, of.interest[i], mode="out")
        cur.ends <- intersect(all.nodes[reachable], of.interest)
        paths <- shortest_paths(g, from=of.interest[i], to=cur.ends, output="vpath")$vpath

        for (p in seq_along(paths)) {
            curpath <- names(paths[[p]])
            middle <- curpath[-c(1L, length(curpath))]
            if (any(of.interest %in% middle)) {
                paths[p] <- list(NULL)
                next
            }

            last <- curpath[length(curpath)]
            if (which(of.interest==last) >= i) {
                paths[p] <- list(NULL)
                next
            }

            paths[[p]] <- curpath
        }

        collected[[i]] <- paths 
    }

    output <- unlist(collected, recursive=FALSE)
    output <- output[!vapply(output, is.null, FALSE)]

    # Sticking back the empties.
    c(output, as.list(all.nodes[deg==0]))
}
