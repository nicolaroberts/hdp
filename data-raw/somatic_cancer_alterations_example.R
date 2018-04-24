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

# initialise HDP
hdp_mut <- hdp_init(ppindex = c(0, rep(1, 3), rep(2:4, each=100)),
                cpindex = c(1, rep(2, 3), rep(3:5, each=100)),
                hh=rep(1, 96),
                alphaa=rep(1, 5),
                alphab=rep(1, 5))

# add data to leaf nodes (one per cancer sample, in row order of mut_count)
hdp_mut <- hdp_setdata(hdp_mut, 5:numdp(hdp_mut), mut_count)

# multiple indep posterior sampling chains
chlist <- vector("list", 4)

for (i in 1:4){

  # activate DPs, 10 initial components
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), 10, seed=i*200)

  chlist[[i]] <- hdp_posterior(hdp_activated,
                               burnin=5000,
                               n=50,
                               space=200,
                               cpiter=3,
                               seed=i*1e3)
}

# example multi object
mut_example_multi <- hdp_multi_chain(chlist)

# save to data/
devtools::use_data(mut_example_multi, overwrite = TRUE)


par(mfrow=c(2,2))
lapply(chains(mut_example_multi), plot_lik, start=1000)
lapply(chains(mut_example_multi), plot_numcluster)
lapply(chains(mut_example_multi), plot_data_assigned)


mut_example_multi <- hdp_extract_components(mut_example_multi)

par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(mut_example_multi, bty="L")

bases <- c("A", "C", "G", "T")
trinuc_context <- paste0(rep(rep(bases, times=6), each=4),
                         rep(c("C", "T"), each=48),
                         rep(bases, times=24))
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
for (i in 0:(length(comp_categ_counts(mut_example_multi)) -1)) {
  print(plot_comp_distn(mut_example_multi, comp=i, cat_names=trinuc_context,
                        grouping=group_factor, col=RColorBrewer::brewer.pal(6, "Set2"),
                        col_nonsig="grey80", show_group_labels=TRUE))
}

plot_dp_comp_exposure(mut_example_multi, main_text="Lung adenocarcinoma",
                      dpindices=4+(1:100),
                      col=RColorBrewer::brewer.pal(12, "Set3"),
                      incl_nonsig=FALSE, incl_comp0 = FALSE)

plot_dp_comp_exposure(mut_example_multi, main_text="Ovarian cancer",
                      dpindices=104+(1:100),
                      col=RColorBrewer::brewer.pal(12, "Set3"),
                      incl_nonsig=FALSE, incl_comp0 = F)

plot_dp_comp_exposure(mut_example_multi, main_text="Melanoma",
                      dpindices=204+(1:100),
                      col=RColorBrewer::brewer.pal(12, "Set3"),
                      incl_nonsig=FALSE, incl_comp0 = F)

plot_dp_comp_exposure(mut_example_multi,
                      dpindices=2:4, incl_numdata_plot=FALSE,
                      col=RColorBrewer::brewer.pal(12, "Set3"),
                      incl_nonsig=T,
                      dpnames=c("Lung Adeno", "Ovarian", "Melanoma"))
