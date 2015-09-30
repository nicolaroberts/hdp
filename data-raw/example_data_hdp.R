# small example dataset
# 10 samples from two groups
# 2 signatures across 6 possible values
# different signature ratios in the two groups

set.seed(10)

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

# save to data/ dir
devtools::use_data(example_data_hdp)

# library(hdp)
# my_hdp <- hdp_init(ppindex=0, cpindex=1, hh=rep(1, 6), alphaa=rep(1, 3), alphab=rep(2, 3))
# my_hdp <- hdp_adddp(my_hdp, 2, 1, 2)
# my_hdp <- hdp_adddp(my_hdp, 10, c(rep(2, 5), rep(3, 5)), 3)
# my_hdp <- hdp_setdata(my_hdp, 4:13, example_data_hdp)
# my_hdp <- dp_activate(my_hdp, 1:13, 2)
#
# input_list <- list(my_hdp, my_hdp, my_hdp, my_hdp, my_hdp)
#
# output <- parallel::mclapply(input_list, hdp_posterior, burnin=1000, n=50, space=50, cpiter=5)
#
# example_multi <- hdp_multi_chain(output)
#
# # save to data/ dir
# devtools::use_data(example_multi, overwrite=TRUE)


