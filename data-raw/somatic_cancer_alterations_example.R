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

# check diagnostics
par(mfrow=c(2,2))
lapply(chains(mut_example_multi), plot_lik, start=1000)
lapply(chains(mut_example_multi), plot_numcluster)
lapply(chains(mut_example_multi), plot_data_assigned)

# save to data/
devtools::use_data(mut_example_multi, overwrite = TRUE)



##### lung, conditioned on prior signatures

cosmic.sigs <- read.table('http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt', header=TRUE, sep='\t')
#  sort by Substitution Type and Trinucleotide
cosmic.sigs <- cosmic.sigs[order(cosmic.sigs$Substitution.Type, cosmic.sigs$Trinucleotide),]
sigs <- as.matrix(cosmic.sigs[,grep('Signature', colnames(cosmic.sigs))])
# number of prior signatures to condition on (30)
nps <- ncol(sigs)

luad_prior <- hdp_prior_init(prior_distn = sigs,
                             prior_pseudoc = rep(1000, nps),
                             hh=rep(1, 96),
                             alphaa=c(1, 1),
                             alphab=c(1, 1))

luad_prior <- hdp_addconparam(luad_prior,
                              alphaa = c(1,1),
                              alphab = c(1,1))

luad_prior <- hdp_adddp(luad_prior,
                        numdp = 101,
                        ppindex = c(1, rep(1+nps+1, 100)),
                        cpindex = c(3, rep(4, 100)))

luad_prior <- hdp_setdata(luad_prior,
                          dpindex = (1+nps+1)+1:100,
                          mut_count[1:100,])

chlist <- vector("list", 4)
for (i in 1:4){
  luad_activated <- dp_activate(luad_prior,
                                dpindex = (1+nps+1)+0:100,
                                initcc = nps+5,
                                seed = i*1000)

  chlist[[i]] <- hdp_posterior(luad_activated,
                               burnin = 4000,
                               n = 50,
                               space = 100,
                               cpiter = 3,
                               seed = i*1e6)
}

luad_multi <- hdp_multi_chain(chlist)

par(mfrow=c(2,2))
p1 <- lapply(chains(luad_multi), plot_lik, bty='L', start=1000)
p2 <- lapply(chains(luad_multi), plot_numcluster, bty='L')
p3 <- lapply(chains(luad_multi), plot_data_assigned, bty='L')

# save to data/
devtools::use_data(luad_multi, overwrite = TRUE)

