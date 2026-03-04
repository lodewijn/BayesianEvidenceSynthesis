# clear console and env
rm(list = ls(all.names = TRUE))
cat("\014")

# reset results folder
unlink("results", recursive = TRUE)
dir.create("results", recursive = TRUE)
dir.create("sim", recursive = TRUE)

# load packages
dependencies <- c(
  "MASS", "bain", "metafor", "lme4",
  "mvtnorm", "doSNOW", "data.table"
)
lapply(dependencies, library, character.only = TRUE)

# Create hypergrid with simulation parameters and save it as .RData file extension
hyper_parameters <- list(
  ndataset = 1:1000,
  hyp_type = c("original", "sym0", "pos08", "neg08"),
  scale = c("cor", "cov"),
  errorsd = c(0.81, 0.5, 0),
  n = c(20, 80, 200, 500),
  k = c(2, 3, 10)
)

summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)

# derive es and hyp_val from hyp_type
summarydata$es <- NA_real_
summarydata$hyp_val <- NA_real_

summarydata$es[summarydata$hyp_type == "original"] <- 0.1
summarydata$hyp_val[summarydata$hyp_type == "original"] <- 0.1

summarydata$es[summarydata$hyp_type == "sym0"] <- 0
summarydata$hyp_val[summarydata$hyp_type == "sym0"] <- 0

summarydata$es[summarydata$hyp_type == "pos08"] <- 0.8
summarydata$hyp_val[summarydata$hyp_type == "pos08"] <- 0.8

summarydata$es[summarydata$hyp_type == "neg08"] <- -0.8
summarydata$hyp_val[summarydata$hyp_type == "neg08"] <- -0.8

# seeds
set.seed(6164900)
summarydata$seed <- sample(1:.Machine$integer.max, nrow(summarydata))

saveRDS(summarydata, "./sim/summarydata.RData")

# prepare parallel processing
library(doSNOW)
nclust <- parallel::detectCores()
cl <- makeCluster(nclust)
registerDoSNOW(cl)

# add progression bar
pb <- txtProgressBar(min = 0, max = nrow(summarydata), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

# run simulation
tab <- foreach(
  rownum = seq_len(nrow(summarydata)),
  .options.snow = opts,
  .packages = c("bain", "mvtnorm", "lme4", "metafor"),
  .combine = rbind
) %dopar% {
  
# set seed
  row <- summarydata[rownum, ]
  set.seed(row$seed)
  
  n <- row$n
  k <- row$k
  es <- row$es
  hyp_val <- row$hyp_val
  errorsd <- row$errorsd
  scale <- row$scale
  
  Sigma <- matrix(c(1, es, es, 1), 2, 2)
  
  # generate datasets
  dfs <- lapply(seq_len(k), function(i) {
    df <- rmvnorm(n, sigma = Sigma)
    noise <- matrix(rnorm(2 * n, sd = errorsd), ncol = 2)
    df + noise
  })
  
  # IPD meta-analysis
  df_join <- scale(do.call(rbind, dfs))
  df_join <- data.frame(
    X1 = df_join[, 1],
    X2 = df_join[, 2],
    X3 = rep(seq_len(k), each = n)
  )
  
  res_ipd <- lmer(X1 ~ X2 + (1 | X3), data = df_join)
  ci_ipd <- confint(res_ipd, level = 0.90)
  
  # obtain estimates for effect size and their standard errors for every dataset
  res <- sapply(dfs, function(x) {
    x <- as.matrix(x)
    # correlation
    if (scale == "cor") {
      est <- cor(x)[2, 1]
      se <- sqrt((1 - est^2) / (n - 2))
      
      # covariance
    } else {
      est <- cov(x)[2, 1]
      se <- sqrt((var(x[,1]) * var(x[,2]) + est^2) / (n - 1))
    }
    c(est, se)
  })
  
  df_meta <- data.frame(t(res))
  res_met <- tryCatch(
    rma(yi = df_meta[,1], sei = df_meta[,2], level = 90),
    error = function(e) NULL
  )
  
  # classical vote count
  classic_votecount <- round(mean(apply(res, 2, function(x) {
    (x[1] - hyp_val) / x[2] > 1.644854
  })))
  
  # Bayes factors
  colnames(res) <- paste0("r", seq_len(k))
  sig <- lapply(res[2, ], function(x) matrix(x^2))
  
  #run bf_individual to extract product bf and geometric product bf
  bf_individual <- lapply(
    paste0(colnames(res), ">", hyp_val), # for every group, we hypothesize that r_k > hyp_val
    bain,                 # call bain
    x = res[1,],          # estimates
    Sigma = sig,          # all rho's are assumed to be independent
    n = rep(n, k),        # pass the named vector of sample sizes per group to bain
    group_parameters = 1, # every group k has 1 parameter which is rho_xy
    joint_parameters = 0) # they do not share parameters 

  # extract BF_ic and BF_iu for the parameter of every group
  BFs <- t(sapply(bf_individual, function(x)
    c(x$fit$BF.c[1], x$fit$BF.u[1])
  ))

  # obtain geometric product and regular product
  gp_and_prod <- apply(BFs, 2, function(x){
    prod_bf <- prod(x)         #obtain product bf
    c(prod_bf^(1/k), prod_bf)  #concatenate geometric product and regular product
  })
  
  # create bf_together object to obtain the BFs for the (group1, group2, group3) > hyp_val
  bf_together <- bain(
    res[1, ],
    hypothesis = paste0("(", paste(colnames(res), collapse = ", "), ") > ", hyp_val),
    n = sum(rep(n, k)), # n = sum of sample sizes over all groups.
    Sigma = diag(res[2, ]^2) # assume independence between groups, square standard errors
  )
  
  # returns in order: gpbf_ic, gpbf_iu, prodbf_ic, prodbf_iu, tbf_ic, tbf_iu
  c(
    rownum,
    c(0, 10)[(ci_ipd[4,1] > hyp_val) + 1], #return 10 if IPD is significantly greater than hyp, 0 if false
    if (is.null(res_met)) {
      NA # return NA if RMA cannot be calculated (to make sure simulation doesn't crash)
    } else {
      c(0,10)[(res_met$ci.lb > hyp_val) + 1] #return 10 if RMA is significantly greater than hyp, 0 if false
    }
    ,
    c(0, 10)[classic_votecount + 1], #return 10 if classic_votecount = T, 0 if false
    gp_and_prod[1, ],
    gp_and_prod[2, ],
    bf_together$fit$BF.c[1],
    bf_together$fit$BF.u[1]
  )
}

#Close cluster
stopCluster(cl)

# End of simulation -------------------------------------------------------
stop("End of simulation")

# Merge files -------------------------------------------------------------
library(data.table)

# read in the simulation conditions 
res <- readRDS("./sim/summarydata.RData")
setDT(res)

# algorithms (geometric product Bayes Factor, product Bayes Factor, together Bayes Factor)
algorithms <- c("gpbf", "prodbf", "tbf")
hyps <- c("_ic", "_iu") # using both complementary and unconstrained
alg_names <- c(
  "ipd", "rma", "votecount",
  paste0(rep(algorithms, each = 2), hyps)
)

# give appropriate names to the simulation results and omit identification variable 'V1'
tab <- as.data.table(tab)
setnames(tab, c("V1", alg_names))
tab[, V1 := NULL]

# cbind conditions and results
res <- cbind(res, tab)

# write results to .RData and .csv extension and delete .txt files in the results folder.
fwrite(res, file.path("sim", paste0("sim_results_1000it_", Sys.Date(), ".csv")))
saveRDS(res, file.path("sim", paste0("sim_results_1000it_", Sys.Date(), ".RData")))

# END OF FILE
tabres <- res[
  , lapply(.SD, function(x) mean(x > 3)),
  .SDcols = alg_names,
  by = c("hyp_type", "scale", "es", "hyp_val", "errorsd", "n", "k")
]

fwrite(tabres, "tabres.csv")
saveRDS(tabres, "tabres.RData")
#colMeans(res[, .SD > 3, .SDcols = alg_names, .gro])