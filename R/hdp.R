#' @import methods
NULL

#' Fake categorical count data
#'
#' Fake categorical count data with 10 samples and 6 categories.
#' Generated from two underlying categorical data distributions with a
#' different average mixture ratio in the first five samples from the last
#' five samples.
#'
#' @format A numeric count matrix with 10 rows and 6 columns
"example_data_hdp"


#' Posterior sampling chain from TCGA somatic mutation data
#'
#' HDP posterior sampling chain with data from SomaticCancerAlterations package.
#' Categories are the 96 base substitution types defined by local trinucleotide
#' content, and the samples include 291 glioblastoma
#' multiforme (gbm) samples, and 293 kidney chromophobe cancer (kirc) samples.
#' Data was derived from exome-sequencing studies.
#' Each sample was assigned to a unique child DP node, with one parent DP node for
#' gbm, one parent node for kirc samples, and one grandparent DP node at the top level.
#' Chain initialised with 8 clusters, then run through 3000 burn-in iterations before
#' collecting 50 posterior samples with 50 iterations between each.
#'
#' @format A hdpSampleChain object with 50 posterior samples
"tcga_example_chain"

#' Multiple posterior sampling chains with TCGA somatic mutation data
#'
#' Five independent HDP posterior sampling chains with data from SomaticCancerAlterations package.
#' Categories are the 96 base substitution types defined by local trinucleotide
#' content, and the samples include 291 glioblastoma
#' multiforme (gbm) samples, and 293 kidney chromophobe cancer (kirc) samples.
#' Data was derived from exome-sequencing studies.
#' Each sample was assigned to a unique child DP node, with one parent DP node for
#' gbm, one parent node for kirc samples, and one grandparent DP node at the top level.
#' Each of the five chains were initialised with 8 clusters, then run through 3000 burn-in iterations before
#' collecting 50 posterior samples with 50 iterations between each.
#'
#' @format A hdpSampleMulti object with 250 posterior samples, 50 from each chain
"tcga_example_multi"
