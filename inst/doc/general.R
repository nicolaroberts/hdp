## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

knitr::opts_chunk$set(fig.width=6, fig.height=5)

## ----load----------------------------------------------------------------
library(hdp)

## ----toydata-------------------------------------------------------------
example_data_hdp

## ----quick_init----------------------------------------------------------
set.seed(10)
quick_hdp <- hdp_quick_init(example_data_hdp)
quick_hdp

## ----quickchain----------------------------------------------------------
quick_chain <- hdp_posterior(quick_hdp, burnin=100, n=100, space=10, seed=1234)
quick_chain

## ----plotchain-----------------------------------------------------------
plot_lik(quick_chain, bty="L")

plot_numcluster(quick_chain, bty="L")

plot_data_assigned(quick_chain, bty="L")

## ----cullpost------------------------------------------------------------
quick_chain_v2 <- cull_posterior_samples(quick_chain, ncull=20)
quick_chain_v2

## ----extract_single------------------------------------------------------
quick_chain_v2 <- hdp_extract_components(quick_chain_v2)
quick_chain_v2

## ----plotcompsize--------------------------------------------------------
plot_comp_size(quick_chain_v2, bty="L", lab=c(3, 5, 7))

## ----plotcompdistn-------------------------------------------------------
par(mfrow=c(2,2), mar=c(3, 2, 2, 1))
plot_comp_distn(quick_chain_v2, cat_names=paste0("Categ", 1:6), col="skyblue3")

## ----plotdpexposure------------------------------------------------------
plot_dp_comp_exposure(quick_chain_v2, dpindices=2:6, main_text="First five samples",
                      col=RColorBrewer::brewer.pal(3, "Set1"))
plot_dp_comp_exposure(quick_chain_v2, dpindices=7:11, main_text="Last five samples",
                      col=RColorBrewer::brewer.pal(3, "Set1"))

## ----hdpinit-------------------------------------------------------------
my_hdp <- hdp_init(ppindex=c(0, 1, 1), 
                   cpindex=c(1, 2, 2), 
                   hh=rep(1, 6), 
                   alphaa=rep(2, 2), 
                   alphab=rep(0.5, 2))
my_hdp


## ----add-----------------------------------------------------------------
my_hdp <- hdp_addconparam(my_hdp, 
                          alphaa=rep(2, 2), 
                          alphab=rep(0.5, 2))
my_hdp <- hdp_adddp(my_hdp, 
                    numdp=10, 
                    ppindex=rep(2:3, each=5), 
                    cpindex=rep(3:4, each=5))
my_hdp

## ----setdata-------------------------------------------------------------
my_hdp <- hdp_setdata(my_hdp,
                      dpindex=4:13,
                      data=example_data_hdp)
my_hdp

## ----activate------------------------------------------------------------
my_hdp <- dp_activate(my_hdp,
                      dpindex=1:13,
                      initcc=2,
                      seed=5678)
my_hdp

## ----runmultichain-------------------------------------------------------
chlist <- vector("list", 4)

for (i in 1:4){
  chlist[[i]] <- hdp_posterior(my_hdp, 
                               burnin=1000,
                               n=50,
                               space=50,
                               cpiter=3, 
                               seed=i*1000)
}

chlist[1:2]

multi <- hdp_multi_chain(chlist)
multi

## ----plotchainmulti, fig.width=7, fig.height=5---------------------------
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(multi), plot_lik, bty="L")
p2 <- lapply(chains(multi), plot_numcluster, bty="L")
p3 <- lapply(chains(multi), plot_data_assigned, bty="L")

## ----extractmulti--------------------------------------------------------
multi <- hdp_extract_components(multi)
multi

plot_comp_size(multi, bty="L", lab=c(3, 5, 7))

par(mfrow=c(2,2), mar=c(3, 2, 2, 1))
plot_comp_distn(multi, cat_names=paste0("Categ", 1:6), col="skyblue3")

par(mfrow=c(1,1))
plot_dp_comp_exposure(multi, dpindices=4:8, main_text="First five samples",
                      col=RColorBrewer::brewer.pal(3, "Set1"))
plot_dp_comp_exposure(multi, dpindices=9:13, main_text="Last five samples",
                      col=RColorBrewer::brewer.pal(3, "Set1"))

## ----hdpstate------------------------------------------------------------

# number of data categories
numcateg(quick_hdp)

# number of DP nodes
numdp(quick_hdp)

# number of concentration parameters
numconparam(quick_hdp)

# 'base' distribution above the top DP node
base(quick_hdp)

# concentration parameter details
conparam(quick_hdp)

# DP node detail
dp(quick_hdp)[1:5]

# DP 'state' for each node
# 2 is activated (included in posterior sampling)
# 1 is frozen (conditioned on during posterior sampling)
# 0 is heldout (ignored during postering sampling)
dpstate(quick_hdp)

# index of parent node for each node
ppindex(quick_hdp)

# index of concentration parameter for each node
cpindex(quick_hdp)


## ----hdpsamplechain------------------------------------------------------

# random seed
hdp_seed(quick_chain)

# settings of the posterior sampling chain
hdp_settings(quick_chain)

# instance of the hdpState object at the end of the chain
final_hdpState(quick_chain)

# data likelihood (given model) after every iteration
lik(quick_chain)[1:10]

# number of raw clusters in each posterior sample
numcluster(quick_chain)[1:10]

# concentration parameter values at each posterior sample
cp_values(quick_chain)[1:10]

# List of matrices (one from each posterior sample) counting the category-cluster 
# data assignment across all DP nodes. Number of rows is the number of 
# categories (constant), and number of columns is the number of clusters 
# in that posterior sample (variable).
clust_categ_counts(quick_chain)[1:3]

# List of matrices (one from each posterior sample) counting within-DP cluster 
# assignment (aggregating across data categories). Number of rows is the number of
# DPs (constant), and number of columns is the number of clusters in that posterior sample (variable).
clust_dp_counts(quick_chain)[1:3]


## ----chainsaccess--------------------------------------------------------
chlist <- chains(multi)
chlist[3:4]

## ----comp_access---------------------------------------------------------
# List of matrices (one for each component) counting the sample-category data 
# assignment across all DP nodes. Number of rows is the number of posterior samples, 
# and number of columns is the number of data categories.
comp_categ_counts(multi)[[2]][1:10,]

# List of matrices (one for each DP) counting sample-component assignment 
# (aggregating across data categories). Number of rows is the number of posterior 
# samples, and number of columns is the number of components.
comp_dp_counts(multi)[[2]][1:10,]

# List with elements "mean" and "cred.int", containing matrices with the mean 
# (and lower/upper 95% credibility interval) over data categories for each component. 
# Number of rows is the number of components, and number of columns is the 
# number of data categories.
comp_categ_distn(multi)

# List with elements "mean" and "cred.int", containing matrices with the mean 
# (and lower/upper 95% credibility interval) distribution over components for each DP. 
# Number of rows is the number of DPs, and number of columns is the number of components.
comp_dp_distn(multi)$mean
comp_dp_distn(multi)$cred.int[4:7]


## ----sessionInfo---------------------------------------------------------
devtools::session_info()

