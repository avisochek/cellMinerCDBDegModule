#' Combine GMT files for gene set enrichment analysis
#' 
#' @param set1 first gene set
#' @param set2 second gene set
#' 
#' @@examples
#' 
#' 
#' @author Nana A. Kusi
#' 
#' 
#' 
#' 
combineGeneSets <- function(set1, set2) {
  # Merge the two lists by combining gene sets with the same key and unioning genes in overlapping sets
  combinedSets <- lapply(union(names(set1), names(set2)), function(x) {
    unique(c(set1[[x]], set2[[x]]))
  })
  names(combinedSets) <- union(names(set1), names(set2))
  return(combinedSets)
}