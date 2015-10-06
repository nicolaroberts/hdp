# example hdpSampleChain and hdpSampleMulti object from real dataset of
# somatic mutation counts in TCGA cancers.

library(hdp)

# load alterations from
# lung adenocarcinoma, ovarian serous cystadenocarcinoma, skin cutaneous melanoma
data(luad_tcga, package='SomaticCancerAlterations')
data(ov_tcga, package='SomaticCancerAlterations')
data(skcm_tcga, package='SomaticCancerAlterations')


# only keep SNP type, add cancer type to sample name, and only keep
# necessary metadata. Only keep 100 samples. Then concatenate and sort.
for (cancer_name in c("luad", "ov", "skcm")){
  raw <- get(paste0(cancer_name, "_tcga"))
  snv <- raw[which(raw$Variant_Type == "SNP")]
  snv <- snv[which(snv$Patient_ID %in% levels(snv$Patient_ID)[1:100])]
  mcols(snv) <- data.frame(sampleID=paste(cancer_name, snv$Patient_ID, sep='_'),
                           ref=snv$Reference_Allele,
                           alt=snv$Tumor_Seq_Allele2)
  assign(cancer_name, snv)
}
variants <- sort(c(luad, ov, skcm))
remove(cancer_name, luad, luad_tcga, ov, ov_tcga, skcm, skcm_tcga, raw, snv)


# tally mutations per sample in 96 base subsitution classes by trinuc context
mut_count <- nrmisc::tally_mutations_96(variants)

# save to data/
devtools::use_data(mut_count, overwrite = TRUE)

# how many of each cancer type?
samp_names <- rownames(mut_count)
cancer_names <- sub('_', '', regmatches(samp_names, regexpr('^.*_', samp_names)))
num_samp <- table(cancer_names)
num_samp

# initialise HDP
ppindex <- c(0,
             rep(1, 3),
             rep(2, num_samp[1]),
             rep(3, num_samp[2]),
             rep(4, num_samp[3]))
cpindex <- c(1,
             rep(2, 3),
             rep(3, num_samp[1]),
             rep(4, num_samp[2]),
             rep(5, num_samp[3]))
hdp <- hdp_init(ppindex,
                cpindex,
                hh=rep(1, 96),
                alphaa=rep(1, 5),
                alphab=rep(1, 5))

# add data
hdp <- hdp_setdata(hdp, 5:numdp(hdp), mut_count)

# activate DPs, 10 initial components
hdp <- dp_activate(hdp, 1:numdp(hdp), 10, seed=1)

# multiple indep posterior sampling chains
chlist <- vector("list", 4)

for (i in 1:4){
  chlist[[i]] <- hdp_posterior(hdp,
                               burnin=4000,
                               n=50,
                               space=50,
                               cpiter=3,
                               seed=i*1e4)
}

# one example chain
mut_example_chain <- chlist[[1]]

# save to data/
devtools::use_data(mut_example_chain, overwrite = TRUE)

# example multi object
mut_example_multi <- hdp_multi_chain(chlist)

# save to data/
devtools::use_data(mut_example_multi, overwrite = TRUE)
