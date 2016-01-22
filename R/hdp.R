#' @import methods
NULL

.onUnload <- function(libpath){
  library.dynam.unload("hdp", libpath)
}

#' Fake categorical count data
#'
#' Fake categorical count data with 10 samples and 6 categories.
#' Generated from two underlying categorical data distributions with a
#' different average mixture ratio in the first five samples from the last
#' five samples.
#'
#' @format A numeric count matrix with 10 rows and 6 columns
"example_data_hdp"

#' Cancer mutation count data
#'
#' Mutation count data from SomaticCancerAlterations package.
#' Categories are the 96 base substitution types defined by local trinucleotide
#' content, and the samples include 100 lung adenocarcinomas, 100 ovarian serous
#' carcinomas, and 100 skin cutaneous melanomas.
#' Data is derived from TCGA exome-sequencing studies.
#'
#' @format A matrix of mutation counts with 300 rows (one per cancer sample)
#'  and 96 columns (one per mutation category)
"mut_count"


#' Lung squamous cell carcinoma mutation data
#'
#' Mutation count data from SomaticCancerAlterations package.
#' Categories are the 96 base substitution types defined by local trinucleotide
#' content, and the samples are 100 lung squamous cell carcinomas.
#' Data is derived from TCGA exome-sequencing studies.
#'
#' @format A matrix of mutation counts with 100 rows (one per cancer sample)
#'  and 96 columns (one per mutation category)
"lusc_count"



#' Posterior sampling chain with cancer mutation data
#'
#' HDP posterior sampling chain with data from SomaticCancerAlterations package.
#' Categories are the 96 base substitution types defined by local trinucleotide
#' content, and the samples include 100 lung adenocarcinomas, 100 ovarian serous
#' carcinomas, and 100 skin cutaneous melanomas.
#' Data is derived from TCGA exome-sequencing studies.
#' Each sample was assigned to a unique child DP node, with one parent DP node per
#' cancer type, and one grandparent DP node at the top level.
#' Chain initialised with 10 clusters, then run through 4000 burn-in iterations
#' before collecting 50 posterior samples with 50 iterations between each.
#'
#' @format A hdpSampleChain object with 50 posterior samples
"mut_example_chain"

#' Multiple posterior sampling chains with cancer mutation data
#'
#' Four independent HDP sampling chains with data from
#' SomaticCancerAlterations package.
#' Categories are the 96 base substitution types defined by local trinucleotide
#' content, and the samples include 100 lung adenocarcinomas, 100 ovarian serous
#' carcinomas, and 100 skin cutaneous melanomas.
#' Data is derived from TCGA exome-sequencing studies.
#' Each sample was assigned to a unique child DP node, with one parent DP node per
#' cancer type, and one grandparent DP node at the top level.
#' Each chain initialised with 10 clusters, then run 4000 burn-in iterations
#' before collecting 50 posterior samples with 50 iterations between each.
#'
#' @format A hdpSampleMulti object with 200 posterior samples, 50 from each chain
"mut_example_multi"


#' Posterior sampling chains for LUSC data, conditional on previous HDP
#'
#' Four independent HDP sampling chains for lung squamous cell carcinoma
#' data, conditioned on previous HDP output (mut_example_multi) for three other
#' cancer types: lung adenocarcinoma, ovarian cancer, and melanoma.
#' Categories are the 96 base substitution types defined by local trinucleotide
#' content.
#' Nodes 2:304 are frozen (old data from other cancer types),
#' node 305 is the parent node for the LUSC samples,
#' and nodes 306:405 hold the data for the 100 new LUSC samples.
#'
#' @format A hdpSampleMulti object with 200 posterior samples, 50 from each chain
"lusc_multi"
