# clear console and env
rm(list = ls(all.names = TRUE))
cat("\014")

# Delete folder where previous results are in and create new 'Results' folder
unlink("results", recursive = TRUE)
dir.create("results")

# load necessary packages
dependencies <- c("MASS", "bain", "metafor", "lme4", "mvtnorm")
lapply(dependencies, function(x) library(x, character.only = TRUE))

# Check package versions
versions <- c(
  compareVersion(as.character(packageVersion("bain")), "0.2.8"),
  compareVersion(as.character(packageVersion("MASS")), "7.3.58")
)
if (!all(versions == 0)) stop("Using the incorrect version of one or more packages.")

# set simulation conditions
hyper_parameters <- list(
  ndataset     = 1:1000,
  rho_type     = c("ordered", "equal"),
  errorsd_type = c("equal_0", "equal_05", "unequal"),
  n            = c(20, 80, 200, 500),
  k            = c(2, 3, 10)
)

# Create hypergrid with simulation parameters
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)

# store errorsd values per variable (X, Y, Z)
summarydata$errorsd_X <- ifelse(summarydata$errorsd_type == "equal_0",  0,
                                ifelse(summarydata$errorsd_type == "equal_05", 0.5, 0.81))
summarydata$errorsd_Y <- ifelse(summarydata$errorsd_type == "equal_0",  0,
                                ifelse(summarydata$errorsd_type == "equal_05", 0.5, 0))
summarydata$errorsd_Z <- ifelse(summarydata$errorsd_type == "equal_0",  0,
                                ifelse(summarydata$errorsd_type == "equal_05", 0.5, 0.5))

# seeds
set.seed(6164900)
summarydata$seed <- sample(1:.Machine$integer.max, nrow(summarydata))
saveRDS(summarydata, file = "./sim/summarydata.RData")

# prepare parallel processing
library(doSNOW)
nclust <- parallel::detectCores()
cl <- makeCluster(nclust)
registerDoSNOW(cl)

# add progress bar
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
  
  n        <- row$n
  k        <- row$k
  errorsds <- c(row$errorsd_X, row$errorsd_Y, row$errorsd_Z)
  
  # rho1 = cor(X,Y), rho2 = cor(X,Z), rho3 = cor(Y,Z)
  rhos <- switch(row$rho_type,
                 ordered = c(rho1 = -0.1, rho2 = 0.0, rho3 = 0.1),
                 equal   = c(rho1 =  0.1, rho2 = 0.1, rho3 = 0.1)
  )
  
  # 3x3 correlation matrix
  Sigma <- matrix(c(
    1,         rhos[1],    rhos[2],
    rhos[1],   1,          rhos[3],
    rhos[2],   rhos[3],    1
  ), 3, 3)
  
  # generate k datasets with 3 columns (X, Y, Z)
  datasets <- lapply(seq_len(k), function(i) {
    
    
    df <- mvtnorm::rmvnorm(n, sigma = Sigma)
    
    # add measurement error (different errorsd per variable)
    noise <- sapply(errorsds, function(sd) rnorm(n, sd = sd))
    df <- df + noise
    
    # extract correlations
    R <- cor(df)
    r1 <- R[1, 2]  # cor(X, Y) = rho1
    r2 <- R[1, 3]  # cor(X, Z) = rho2
    r3 <- R[2, 3]  # cor(Y, Z) = rho3
    
    # and SEs
    se1 <- sqrt((1 - r1^2) / (n - 2))
    se2 <- sqrt((1 - r2^2) / (n - 2))
    se3 <- sqrt((1 - r3^2) / (n - 2))
    
    list(
      est = c(rho1 = r1, rho2 = r2, rho3 = r3),
      se  = c(rho1 = se1, rho2 = se2, rho3 = se3)
    )
  })
  
  # test ordered hypothesis once per dataset
  BFs_ordered <- t(sapply(datasets, function(d) {
    
    est <- d$est
    sig <- lapply(d$se, function(se) matrix(se^2))
    names(sig) <- names(est)
    
    # Hi: rho3 > rho2 > rho1
    bf <- bain(
      x                = est,
      hypothesis       = "rho3 > rho2 > rho1",
      Sigma            = sig,
      n                = rep(n, 3),
      group_parameters = 1,
      joint_parameters = 0
    )
    
    c(bf$fit$BF.c[1], bf$fit$BF.u[1])
  }))
  
  # product BF across k datasets
  prod_bf_ic <- prod(BFs_ordered[, 1])
  prod_bf_iu <- prod(BFs_ordered[, 2])
  
  # geometric mean product BF
  gp_bf_ic <- prod_bf_ic^(1 / k)
  gp_bf_iu <- prod_bf_iu^(1 / k)
  
  fits <- c(rownum, gp_bf_ic, gp_bf_iu, prod_bf_ic, prod_bf_iu)
  
  write.table(
    x         = t(fits),
    file      = sprintf("results/results_%d.txt", Sys.getpid()),
    sep       = "\t",
    append    = TRUE,
    row.names = FALSE,
    col.names = FALSE
  )
  NULL
}

stopCluster(cl)

# End of simulation -------------------------------------------------------
stop("End of simulation")

# Merge files -------------------------------------------------------------
library(data.table)

res <- readRDS("./sim/summarydata_ordered.RData")
setDT(res)

f   <- list.files("results/", full.names = TRUE)
tab <- rbindlist(lapply(f, fread))
setorderv(tab, cols = "V1", order = 1L, na.last = FALSE)

alg_names <- c("gpbf_ordered_ic", "gpbf_ordered_iu",
               "prodbf_ordered_ic", "prodbf_ordered_iu")

if (!(tab[1, 1] == 1 & tab[nrow(tab), 1] == nrow(res) & length(unique(tab$V1)) == nrow(res))) {
  stop("Results not the same length as number of simulation iterations")
}

names(tab) <- c("V1", alg_names)
tab[, "V1" := NULL]

res <- cbind(res, tab)
rm(tab)

fwrite(res, file.path("sim", paste0("sim_results_ordered_", Sys.Date(), ".csv")))
saveRDS(res, file.path("sim", paste0("sim_results_ordered_", Sys.Date(), ".RData")))

# summary: proportion of iterations where BF > 3
tabres <- res[, lapply(.SD, function(x) mean(x > 3)),
              .SDcols = alg_names,
              by = c("errorsd_type", "n", "k")]
write.csv(tabres, "tabres_ordered.csv", row.names = FALSE)
saveRDS(tabres, "tabres_ordered.RData")