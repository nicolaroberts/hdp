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



## lung squamous cell carcinoma, condition on previous chains


data("lusc_tcga", package="SomaticCancerAlterations")

lusc <- lusc_tcga[which(lusc_tcga$Variant_Type == "SNP")]
lusc <- lusc[which(lusc$Patient_ID %in% levels(lusc$Patient_ID)[1:101])]
mcols(lusc) <- data.frame(sampleID=paste('lusc', lusc$Patient_ID, sep='_'),
                          ref=lusc$Reference_Allele,
                          alt=lusc$Tumor_Seq_Allele2)

remove(lusc_tcga)

# tally mutation counts in 96 base substitution classes defined by trinucleotide context
lusc_count <- nrmisc::tally_mutations_96(lusc)

# remove sample with the largest burden (very different to others)
lusc_count <- lusc_count[-which.max(rowSums(lusc_count)),]
head(lusc_count[,1:5])

# save to data/
devtools::use_data(lusc_count, overwrite = TRUE)

# run four chains off end of previous
hdpStatelist <- lapply(chains(mut_example_multi), final_hdpState)
chlist <- vector("list", 4)

for (i in 1:4){
  hdp <- hdpStatelist[[i]]
  hdp <- hdp_addconparam(hdp, 1, 1)
  hdp <- hdp_adddp(hdp, 101,
                   ppindex=c(1, rep(305, 100)),
                   cpindex=c(2, rep(6, 100)))
  hdp <- hdp_setdata(hdp, 306:405, lusc_count)
  hdp <- dp_freeze(hdp, 2:304)
  hdp <- dp_activate(hdp, 305:405, initcc=base(hdp)@numclass, seed=i*1e5)

  chlist[[i]] <- hdp_posterior(hdp, burnin=1500, n=50,
                               space=50, cpiter=3, seed=i*1e6)
}

lusc_multi <- hdp_multi_chain(chlist)
# save to data/
devtools::use_data(lusc_multi, overwrite = TRUE)
