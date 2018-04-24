# small example dataset
# 10 samples from two groups
# 2 signatures across 6 possible values
# different signature ratios in the two groups

set.seed(12)

# signatures
sigA <- c(1, 1, 2, 5, 2, 1)
sigB <- c(5, 1, 1, 1, 1, 2)
sigs <- rbind(sigA/sum(sigA), sigB/sum(sigB))

# per-sample probability of exposure to each sig
sig_exposure_per_sample <- rbind(gtools::rdirichlet(5, c(3, 1)), gtools::rdirichlet(5, c(1, 3)))

# vector of probs per sample (exposures*sigs)
data_distn_persamp <- sig_exposure_per_sample %*% sigs
# store as list for easy iteration through
data_distn_persamp_list <- plyr::alply(data_distn_persamp, 1)

# number of data items per sample (10 samples)
num_data_items <- floor(10^rgamma(10, shape=20, rate=10))

# simulate data
example_data_hdp <- t(mapply(rmultinom, n=1, size=num_data_items, prob=data_distn_persamp_list))
example_data_hdp

# save to data/ dir
devtools::use_data(example_data_hdp, overwrite=TRUE)

rm(example_data_hdp, data_distn_persamp, sig_exposure_per_sample, sigs,
   data_distn_persamp_list, num_data_items, sigA, sigB)

# larger example dataset with some known priors

# 100 samples from one group
# 4 signatures across 10 possible values
# 2 kept as known priors

set.seed(25)

# signatures
sigA <- sample(8, 10, replace=TRUE)
sigB <- sample(10, 10, replace=TRUE)
sigC <- sample(15, 10, replace=TRUE)
sigD <- sample(7, 10, replace=TRUE)
sigs <- cbind(sigA/sum(sigA),
              sigB/sum(sigB),
              sigC/sum(sigC),
              sigD/sum(sigD))

par(mfrow=c(2,2))
apply(sigs, 2, barplot)

# per-sample probability of exposure to each sig
sig_exposure_per_sample <- gtools::rdirichlet(100, c(0.7, 0.7, 0.7, 0.7))

# vector of probs per sample (exposures*sigs)
data_distn_persamp <- sig_exposure_per_sample %*% t(sigs)
# store as list for easy iteration through
data_distn_persamp_list <- plyr::alply(data_distn_persamp, 1)

# number of data items per sample (100 samples)
num_data_items <- floor(10^rgamma(100, shape=20, rate=10))

# simulate data
example_data_hdp_prior <- t(mapply(rmultinom, n=1, size=num_data_items, prob=data_distn_persamp_list))

# known prior sigs
example_known_priors <- sigs[,1:2]

# save to data/ dir
devtools::use_data(example_data_hdp_prior, overwrite=TRUE)
devtools::use_data(example_known_priors, overwrite=TRUE)


# extract components
hdp_p <- hdp_prior_init(example_known_priors, rep(1000, 2), hh=rep(1, 10),
             alphaa=c(1,1), alphab=c(1,1))
hdp_p <- hdp_addconparam(hdp_p, alphaa=c(1,1), alphab=c(1,1))
hdp_p <- hdp_adddp(hdp_p, 101, c(1, rep(4, 100)), c(3, rep(4, 100)))
hdp_p <- hdp_setdata(hdp_p, 5:104, example_data_hdp_prior)

chlist <- vector("list", 4)

for (i in 1:4){

  activated_hdp <- dp_activate(hdp_p, dpindex=4:104, initcc=4, seed=i*1e3)

  chlist[[i]] <- hdp_posterior(activated_hdp,
                               burnin=2500,
                               n=50,
                               space=100,
                               cpiter=2,
                               seed=i*1e5)
}

chlist

hdp_p <- hdp_multi_chain(chlist)
devtools::use_data(hdp_p, overwrite=TRUE)
