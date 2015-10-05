# example hdpSampleChain and hdpSampleMulti object from real dataset of
# somatic mutation counts in TCGA cancers.

library(hdp)

# load alterations from
# glioblastoma multiforme (gbm), and kidney chromophobe cancer (kirc).
data(gbm_tcga, package='SomaticCancerAlterations')
data(kirc_tcga, package='SomaticCancerAlterations')


# only keep SNP type, add cancer type to sample name, and only keep
# necessary metadata. Then concatenate and sort.
for (cancer_name in c("gbm", "kirc")){
  raw <- get(paste0(cancer_name, "_tcga"))
  snv <- raw[which(raw$Variant_Type == "SNP")]
  mcols(snv) <- data.frame(sampleID=paste(cancer_name, snv$Patient_ID, sep='_'),
                           ref=snv$Reference_Allele,
                           alt=snv$Tumor_Seq_Allele2)
  assign(cancer_name, snv)
}
variants <- sort(c(gbm, kirc))
remove(cancer_name, kirc, kirc_tcga, gbm, gbm_tcga, raw, snv)


# tally mutations per sample in 96 base subsitution classes by trinuc context
tcga_snv_count_data <- nrmisc::tally_mutations_96(variants)


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
devtools::use_data(tcga_example_chain, overwrite = TRUE)

# multiple indep chains

tcga_example_chain_v2 <- hdp_posterior(hdp, 3000, 25, 50, 5, seed=20)
tcga_example_chain_v3 <- hdp_posterior(hdp, 3000, 25, 50, 5, seed=30)
tcga_example_chain_v4 <- hdp_posterior(hdp, 3000, 25, 50, 5, seed=40)
tcga_example_chain_v5 <- hdp_posterior(hdp, 3000, 25, 50, 5, seed=50)

tcga_example_multi <- hdp_multi_chain(list(
  tcga_example_chain,
  tcga_example_chain_v2,
  tcga_example_chain_v3,
  tcga_example_chain_v4,
  tcga_example_chain_v5))

# save to data/
devtools::use_data(tcga_example_multi, overwrite = TRUE)
