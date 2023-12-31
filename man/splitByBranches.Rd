% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitByBranches.R
\name{splitByBranches}
\alias{splitByBranches}
\title{Split a graph into branch-free paths}
\usage{
splitByBranches(g)
}
\arguments{
\item{g}{A \link{graph} object such as that produced by \code{\link{createClusterMST}}.}
}
\value{
A list of character vectors containins unbranched paths between branch points or terminal nodes.
}
\description{
Split a graph (usually a MST, at least a DAG) into branch-free components, i.e., paths between branch points or terminal nodes.
}
\details{
This function implements another strategy to define paths through the MST,
by simply considering all linear sections (i.e., connected by nodes of degree 2) between branch points or terminii.
The idea is to enable characterisation of the continuum without external information or assumptions statements about where the root is positioned,
which would otherwise be required in functions such as \code{\link{defineMSTPaths}}.
}
\examples{
library(igraph)
test.g <- make_graph(c("A", "B", "B", "C", "C", "D", 
    "D", "E", "D", "F", "F", "G"), directed=FALSE)
splitByBranches(test.g)

}
\seealso{
\code{\link{defineMSTPaths}}, for the root-based method of defining paths.
}
\author{
Aaron Lun
}
