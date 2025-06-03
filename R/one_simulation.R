#' Run a single simulation replicate
#'
#' @return A list of simulation results

one_simulation <- function() {
  
  if (cfg$model_version!=L$model_version) {
    stop("cfg$model_version must equal L$model_version.")
  }
  
  chk(0, "START")
  
  # Generate and transform dataset
  dat <- generate_data(
    n = L$n,
    max_time = L$max_time,
    params = L$par
  )
  
  dat_objs <- transform_dataset(dat=dat, window_start=1, window_end=9999)
  dat_i_names <- names(dat_objs[[1]]$dat_i)
  n <- attr(dat, "n")
  n_batches <- 2
  folds <- cut(c(1:n), breaks=n_batches, labels=FALSE)
  batches <- lapply(c(1:n_batches), function(batch) { c(1:n)[folds==batch] })
  dat_objs_wrapper <- lapply(batches, function(i) { dat_objs[i] })

  chk(1, "Data generated")
  
  # Set initial parameter vector
  if (L$model_version==1) {
    par_init <- log(c(
      a_x=0.055, g_x1=1.2, g_x2=1.1, t_x=1,
      a_s=0.22, g_s1=1.8, g_s2=1.2, t_s=1, 
      beta_x=1.4,
      a_y=0.025, g_y1=1.3, g_y2=1, t_y=1
    ))
  }

  chk(2, "construct_negloglik: START")
  negloglik <- construct_negloglik(parallelize=F, temp=dat_objs_wrapper)
  chk(2, "construct_negloglik: END")
  chk(2, "optim: START")
  opt <- stats::optim(
    par = par_init,
    fn = negloglik,
    method = "Nelder-Mead",
    control = list(maxit=cfg$opt_maxit, reltol=cfg$opt_reltol)
  )
  chk(2, "optim: END")
  chk(2, "hessian: START")
  hessian_est <- numDeriv::hessian(
    func = negloglik,
    x = opt$par,
    method = "Richardson",
    method.args = list(
      eps = 0.0001,
      d = 0.0001,
      zero.tol = sqrt(.Machine$double.eps/7e-7),
      r = cfg$opt_r,
      v = 2
    )
  )
  hessian_inv <- solve(hessian_est)
  chk(2, "hessian: END")
  
  res <- list(conv=opt$convergence, iter=as.numeric(opt$counts[1]))
  pn <- names(par_init)
  for (i in c(1:length(par_init))) {
    res[[paste0("lik_M_",pn[i],"_est")]] <- as.numeric(opt$par[i])
    res[[paste0("lik_M_",pn[i],"_se")]] <- sqrt(diag(hessian_inv))[i]
  }
  
  chk(3, "END")
  
  return(res)
  
}
