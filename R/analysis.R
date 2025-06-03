######################.
##### Print info #####
######################.

chk(1, "Data analysis: START")
.t_start <- Sys.time()
print("")
print("CONFIG")
print("------------")
print(paste("model_version:", cfg$model_version))
print(paste("sample size:", cfg$samp_size))
print(paste("sex:", cfg$model_sex))
print("------------")



######################################################.
##### Create global model objects (sex-specific) #####
######################################################.

if (cfg$model_sex=="Female") {
  par_y <- par_y_F; terms_y <- terms_y_F; terms_y2 <- terms_y2_F;
  par_x <- par_x_F; terms_x <- terms_x_F; terms_x2 <- terms_x2_F;
  par_s <- par_s_F; terms_s <- terms_s_F; terms_s2 <- terms_s2_F;
} else {
  par_y <- par_y_M; terms_y <- terms_y_M; terms_y2 <- terms_y2_M;
  par_x <- par_x_M; terms_x <- terms_x_M; terms_x2 <- terms_x2_M;
  par_s <- par_s_M; terms_s <- terms_s_M; terms_s2 <- terms_s2_M;
}



#############################################.
##### Construct log likelihood function #####
#############################################.

chk(3, "construct_negloglik: START")
n <- attr(dat_prc, "n")
print(paste0("Using ", cfg$sim_n_cores, " cores."))
if (cfg$parallelize) {
  n_batches <- cfg$sim_n_cores
} else {
  n_batches <- 2
}
folds <- cut(c(1:n), breaks=n_batches, labels=FALSE)
batches <- lapply(c(1:n_batches), function(batch) { c(1:n)[folds==batch] })
dat_objs_wrapper <- lapply(batches, function(i) { dat_objs[i] })

if (cfg$parallelize) {
  cl <- parallel::makeCluster(cfg$sim_n_cores)
  objs_to_export <- c(
    "f_x", "f_y", "icll", "lik_fn", "batches", "uncompress", "par_init",
    "par_y", "terms_y", "terms_y2", "par_x", "terms_x", "terms_x2", "par_s",
    "terms_s", "terms_s2"
  )
  parallel::clusterExport(cl, objs_to_export, envir=.GlobalEnv)
  negloglik <- construct_negloglik(parallelize=T, use_counter=T)
} else {
  negloglik <- construct_negloglik(parallelize=F, use_counter=T)
}
chk(3, "construct_negloglik: END")



############################.
##### Run optimization #####
############################.

chk(4, "optim: START")
counter <- 0
st <- system.time({ nll_init <- negloglik(par_init) })
print(st)
print(paste("negloglik(par_init):", nll_init))
opt <- stats::optim(
  par = par_init,
  fn = negloglik,
  method = "Nelder-Mead",
  control = list(maxit=cfg$opt_maxit, reltol=cfg$opt_reltol)
)
print("optim() finished.")
print(opt)

chk(4, "optim: END")



###########################.
##### Compute Hessian #####
###########################.

chk(5, "hessian: START")
parallel::clusterExport(cl, "opt", envir=.GlobalEnv)
hessian_est <- numDeriv::hessian(
  func = negloglik,
  x = opt$par,
  method = "Richardson", # "Richardson" "complex"
  method.args = list(
    eps = 0.0001,
    d = 0.0001, # d gives the fraction of x to use for the initial numerical approximation
    zero.tol = sqrt(.Machine$double.eps/7e-7),
    r = cfg$opt_r, # r gives the number of Richardson improvement iterations (repetitions with successly smaller d
    v = 2 # v gives the reduction factor
  )
)

tryCatch(
  expr = {
    hessian_inv <- solve(hessian_est)
  },
  error = function(e) {
    print("Hessian not invertible.")
    print(hessian_est)
  }
)

chk(5, "hessian: END")
parallel::stopCluster(cl)



##################################.
##### Save and print results #####
##################################.

saveRDS(
  list(opt=opt, hessian_inv=hessian_inv),
  paste0("ests_", cfg$model_version,
         ifelse(cfg$samp_size==0, "_full_", "_partial_"),
         substr(cfg$model_sex,1,1), "_", format(Sys.time(), "%Y%m%d"),
         ".rds")
)

# Parse results
res <- data.frame(
  "par" = character(),
  "par_init" = double(),
  "est" = double(),
  "se" = double(),
  "ci_lo" = double(),
  "ci_up" = double()
)
for (i in c(1:length(par_init))) {
  est <- as.numeric(opt$par[i])
  se <- sqrt(diag(hessian_inv))[i]
  res[i,] <- c(
    par = names(par_init)[i],
    par_init = round(par_init[i],4),
    est = round(est,4),
    se = round(se,4),
    ci_lo = round(est-1.96*se,4),
    ci_up = round(est + 1.96*se,4)
  )
}
print(res)

chk(1, "Data analysis: END")
print("Runtime (analysis):")
print(Sys.time()-.t_start)
