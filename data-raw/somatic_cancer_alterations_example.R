# example hdpSampleChain and hdpSampleMulti object from real dataset of
# somatic mutation counts in TCGA cancers.

library(GenomicRanges)
library(Biostrings)
library(hdp)

# load human reference genome hg19/GRCh37
genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

# load alterations from
# glioblastoma multiforme (gbm), and kidney chromophobe cancer (kirc).
data(gbm_tcga, package='SomaticCancerAlterations')
data(kirc_tcga, package='SomaticCancerAlterations')


# only keep SNP type, add cancer type as metadata column, and only keep
# necessary metadata. Then concatenate and sort.
for (cancer_name in c("gbm", "kirc")){
  raw <- get(paste0(cancer_name, "_tcga"))
  snv <- raw[which(raw$Variant_Type == "SNP")]
  mcols(snv) <- data.frame(patient_ID=paste(cancer_name, snv$Patient_ID, sep='_'),
                           reference=snv$Reference_Allele,
                           alternate=snv$Tumor_Seq_Allele2)
  assign(cancer_name, snv)
}
variants <- sort(c(gbm, kirc))
remove(cancer_name, kirc, kirc_tcga, gbm, gbm_tcga, raw, snv)

# only keep standard chromosomes, and ditch MT
variants <- keepStandardChromosomes(variants)
variants <- variants[seqnames(variants)!='MT']
variants <- dropSeqlevels(variants, 'MT')

# chromosomes with 'chr' start
chrlevels <- paste0('chr', seqlevels(variants))
names(chrlevels) <- seqlevels(variants)
variants <- renameSeqlevels(variants, chrlevels)
remove(chrlevels)

# find the trinuc context of each mutation
bases <- c('A', 'C', 'G', 'T')
trinuc_levels <- paste0(rep(bases, each=16), rep(rep(bases, each=4), 4), rep(bases, 16))

get_trinuc <- function(seqname){
  pos <- start(variants[seqnames(variants)==seqname])
  view <- Views(genome[[seqname]], start=pos-1, end=pos+1)
  ans <- factor(as.character(view), levels=trinuc_levels, labels=1:64)
  return(as.numeric(ans))
}

trinuc <- sapply(seqlevels(variants), get_trinuc)
variants$trinuc <- factor(unlist(trinuc, use.name=FALSE), levels=1:64, labels=trinuc_levels)
remove(trinuc)

# are there any places where the reference doesn't match up?
bad <- which(substr(variants$trinuc, 2, 2) != variants$reference)
# variants <- variants[-bad]

# convert mutations wrt A or G ref base into reverse complement
variants$ref <- variants$reference
variants$alt <- variants$alternate
variants$context <- variants$trinuc
torc <- which(variants$reference %in% c('A', 'G'))
variants$ref[torc] <- as.character(reverseComplement(DNAStringSet(variants$ref[torc])))
variants$alt[torc] <- as.character(reverseComplement(DNAStringSet(variants$alt[torc])))
variants$context[torc] <- as.character(reverseComplement(DNAStringSet(variants$context[torc])))


# define classes
variants$class <- paste(variants$ref, variants$alt, 'in', variants$context, sep='.')
variants$class <- factor(variants$class)

# tally mutation classes in each tumour
tcga_snv_count_data <- table(variants$patient_ID, variants$class)
class(tcga_snv_count_data) <- "matrix"

# how many of each cancer type?
tcga_names <- rownames(tcga_snv_count_data)
cancer_names <- sub('_', '', regmatches(tcga_names, regexpr('^.*_', tcga_names)))
table(cancer_names)

# initialise HDP
ppindex <- c(0, rep(1, 2), rep(2, 291), rep(3, 293))
cpindex <- c(1, rep(2, 2), rep(3, 291), rep(4, 293))
hdp <- hdp_init(ppindex, cpindex, hh=rep(1, 96), alphaa=rep(1, 4), alphab=rep(1, 4))

# add data
hdp <- hdp_setdata(hdp, 4:numdp(hdp), tcga_snv_count_data)

# activate DPs, 8 initial components
hdp <- dp_activate(hdp, 1:numdp(hdp), 8, seed=1)

# posterior sampling
tcga_example_chain <- hdp_posterior(hdp, 3000, 50, 50, 5, seed=10)

# save to data/
devtools::use_data(tcga_example_chain)

# multiple indep chains

tcga_example_chain_v2 <- hdp_posterior(hdp, 3000, 50, 50, 5, seed=20)
tcga_example_chain_v3 <- hdp_posterior(hdp, 3000, 50, 50, 5, seed=30)
tcga_example_chain_v4 <- hdp_posterior(hdp, 3000, 50, 50, 5, seed=40)
tcga_example_chain_v5 <- hdp_posterior(hdp, 3000, 50, 50, 5, seed=50)

tcga_example_multi <- hdp_multi_chain(list(
  tcga_example_chain,
  tcga_example_chain_v2,
  tcga_example_chain_v3,
  tcga_example_chain_v4,
  tcga_example_chain_v5))

# save to data/
devtools::use_data(tcga_example_multi)
