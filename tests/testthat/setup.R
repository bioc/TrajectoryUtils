factor2matrix <- function(g) {
    ugroup <- sort(unique(g))
    sgroups <- matrix(0, length(g), length(ugroup))
    colnames(sgroups) <- ugroup
    sgroups[cbind(seq_along(g), match(g, ugroup))] <- 1
    sgroups
}
