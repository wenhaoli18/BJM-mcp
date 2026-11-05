# Fit the pbc3 demo BJM (same spec as BJM-attribution/code/fit.R).
suppressMessages(library(BJM)); data(pbc3)
markers <- c("serBilir","albumin","prothrombin","platelets")
surv_dat <- pbc3[!duplicated(pbc3$id), ]
survival_fit_all <- survivalSub(surv_dat, form_marginal_surv = Surv(years, status3) ~ age + sex, form_conditional_cr = NULL)
fx <- setNames(lapply(markers, \(m) as.formula(paste(m, "~ year + age + sex + (years) + (years) * year"))), markers)
rx <- setNames(rep(list(~ year | id), length(markers)), markers)
fitdat <- rep(list(pbc3[pbc3$status3 == 1, ]), length(markers))
t0 <- Sys.time()
long_fit_all <- longitudinalSub(fitdat, fx, rx)
cat("fit time", format(Sys.time()-t0), "\n")
out <- file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "..", "models", "fits_pbc3.rds")
saveRDS(list(s=survival_fit_all, l=long_fit_all, markers=markers), out)
print(long_fit_all)
