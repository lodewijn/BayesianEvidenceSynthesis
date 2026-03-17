# Read in original data set and replication (including es = 0)
dat1 <- readRDS(file.path("Sim", "sim_results_2022-05-08.RData"))
dat2 <- readRDS(file.path("Sim", "sim_results_2026-03-10.RData"))

# Check if seeds are the same now
which(dat1$seed != dat2$seed) # all seeds align

# Check if they have the same number of columns and rows 
nrow(dat1) == nrow(dat2) # TRUE
ncol(dat1) == ncol(dat2) # TRUE

# Make sure they have the same column names
which(names(dat1) != names(dat2)) # Column 10 has a different name

# Rename columns of both data sets to make sure they align
names(dat1) <- c("ndataset", "es", "reliability", "n", "k", "hyp_val", "seed", "IPD", 
                 "RMA", "VC", "gpbf_ic", "gpbf_iu", "PBF", "prodbf_iu", 
                 "tbf_ic", "tbf_iu")

names(dat2) <- names(dat1)

# Check if data sets are the same
dat1 == dat2 # They are not in all columns

# Especially the prodbf_iu, and gpbfs differ, but not prodbf_ic

# Test if these differences are actually meaningful 
mean(abs(dat1$prodbf_iu - dat2$prodbf_iu), na.rm = TRUE) # 3.182259e-14
mean(abs(dat1$gpbf_ic - dat2$gpbf_ic), na.rm = TRUE) # 1.986921e-13
mean(abs(dat1$gpbf_iu - dat2$gpbf_iu), na.rm = TRUE) # 1.136935e-15
mean(abs(dat1$tbf_ic - dat2$tbf_ic), na.rm = TRUE) # 2.600731e-13
mean(abs(dat1$tbf_iu - dat2$tbf_iu), na.rm = TRUE) # 3.170125e-14

# It stands out that all mean differences between the columns that are not 
# exactly the same are very small (< 10^-13). This means they are neglectible.

# Test if colums are (nearly) equal
all.equal(dat1$prodbf_iu, dat2$prodbf_iu, tolerance = 1e-6) # TRUE
all.equal(dat1$gpbf_ic, dat2$gpbf_ic, tolerance = 1e-6) # TRUE
all.equal(dat1$gpbf_iu, dat2$gpbf_iu, tolerance = 1e-6) # TRUE
all.equal(dat1$tbf_ic, dat2$tbf_ic, tolerance = 1e-6) # TRUE
all.equal(dat1$tbf_iu, dat2$tbf_iu, tolerance = 1e-6) # TRUE

# They are equal, so this is likely not the reason why the confusion matrices differ.
# This may be due to the fact that MASS 7.3.55 was not available for my version 
# of R (and I was not able to install it any more), so I used MASS 7.3.58.

# Run confusion matrices manually to see where they differ
res1 <- dat1[!es < .1, sapply(.SD, function(var) {
  table(ordered(es > .1, levels = c("FALSE", "TRUE")), 
        ordered(var > 3, levels = c("FALSE", "TRUE")))
}), .SDcols = varsout]
rownames(res1) <- c("TN", "FN", "FP", "TP")

res2 <- dat2[!es < .1, sapply(.SD, function(var) {
  table(ordered(es > .1, levels = c("FALSE", "TRUE")), 
        ordered(var > 3, levels = c("FALSE", "TRUE")))
}), .SDcols = varsout]
rownames(res2) <- c("TN", "FN", "FP", "TP")

res1
res2

# Check if the distribution of PBF values > 3 is equal across the two data sets
table(dat1$PBF > 3)
table(dat2$PBF > 3)

all(dat1$PBF == dat2$PBF) # FALSE

# See if these differences are meaningful
all.equal(dat1$PBF, dat2$PBF, tolerance = 1e-6) # PBF in dat2 appears to be categorical!

# Convert PBF in dat2 to numeric
dat2$PBF <- as.numeric(dat2$PBF)

# Check tables again
table(dat1$PBF > 3) == table(dat2$PBF > 3) # They match perfectly now

all.equal(dat1$PBF, dat2$PBF, tolerance = 1e-6) # TRUE - so columns are the same

# Check if the other variables are in the correct form
sapply(dat1[, varsout, with = FALSE], class)
sapply(dat2[, varsout, with = FALSE], class)

# Change all values to numeric where needed
dat2[, c("IPD", "RMA", "VC") := lapply(.SD, as.numeric), .SDcols = c("IPD", "RMA", "VC")]

# Check results again
res1 <- dat1[!es < .1, sapply(.SD, function(var) {
  table(ordered(es > .1, levels = c("FALSE", "TRUE")), 
        ordered(var > 3, levels = c("FALSE", "TRUE")))
}), .SDcols = varsout]
rownames(res1) <- c("TN", "FN", "FP", "TP")

res2 <- dat2[!es < .1, sapply(.SD, function(var) {
  table(ordered(es > .1, levels = c("FALSE", "TRUE")), 
        ordered(var > 3, levels = c("FALSE", "TRUE")))
}), .SDcols = varsout]
rownames(res2) <- c("TN", "FN", "FP", "TP")

res1
res2

# The PBF results align now, but those of the other methods do not. 

# Maybe this is still due to the MASS difference?
all.equal(dat1$IPD, dat2$IPD, tolerance = 1e-6) # mean diff 10
all.equal(dat1$RMA, dat2$RMA, tolerance = 1e-6) # mean diff 10
all.equal(dat1$VC, dat2$VC, tolerance = 1e-6) # mean diff 10

# Since the only two possible values are 0 or 10, this means that some rows 
# may have shifted between significant/non-significant between data sets.

# Check which rows are differing between data sets
ipd_diff <- which(dat1$IPD != dat2$IPD)
rma_diff <- which(dat1$RMA != dat2$RMA)
vc_diff  <- which(dat1$VC  != dat2$VC)

# What percentage of rows differs between data sets
sum(dat1$IPD != dat2$IPD)/nrow(dat1) * 100 # 2.85%
sum(dat1$RMA != dat2$RMA)/nrow(dat1) * 100 # 2.76 %
sum(dat1$VC != dat2$VC)/nrow(dat1) * 100 # 2.50%

# How many rows differ in the same three variables simultaneously?
length(intersect(ipd_diff, intersect(rma_diff, vc_diff))) # 7
# Only 7 rows differ in all three columns simultaneously.

# This means that it is unlikely that there are structural differences in the
# data used to calculate these estimates (due to MASS versions), 
# but it is possible that the methods for IPD, RMA, and VC are more 
# sensitive to numerical fluctuations.
