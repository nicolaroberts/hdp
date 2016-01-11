#' Extract major components from the raw clusters
#' @param x hdpSampleChain or hdpSampleMulti object
#' @param cos.merge Merge components with cosine similarity above this threshold (default 0.90)
#' @param prop.ex The proportion of data items explained by the components (default 0.97).
#'  Only used for hdpSampleChain as this can be calculated automatically for hdpMultiChain objects.
#'@param redo Logical. If true - constituent chains with previously calculated components
#'  will be re-calculated. Only used for hdpSampleMulti.
#' @return A hdpSampleChain or hdpSampleMulti object updated with component information
#' @aliases hdp_extract_components
#' @import clue
#' @include hdp_extract_comp_single.R hdp_extract_comp_multi.R
#' @export
#' @examples
#' hdp_extract_components(mut_example_chain)

setGeneric("hdp_extract_components",
           function(x, cos.merge=0.9, prop.ex=0.97, redo=TRUE)
             standardGeneric("hdp_extract_components"))

#' @describeIn hdp_extract_components Extract components for hdpSampleChain
setMethod("hdp_extract_components",
          signature = "hdpSampleChain",
          definition = function(x, cos.merge, prop.ex) {
            ans <- hdp_extract_comp_single(x, prop.ex, cos.merge)
            return(ans)
          })

#' @describeIn hdp_extract_components Extract components for hdpSampleMulti
setMethod("hdp_extract_components",
          signature = "hdpSampleMulti",
          definition = function(x, cos.merge, redo) {
            ans <- hdp_extract_comp_multi(x, cos.merge, redo)
            return(ans)
          })
