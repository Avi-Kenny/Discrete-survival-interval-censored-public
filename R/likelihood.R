#' Transform dataset
#' @param dat A dataset returned by generate_dataset()
transform_dataset <- function(dat, window_start, window_end) {
  
  # This procedure assumes the dataset is sorted by id
  row_map <- new.env()
  id_prev <- -1
  i <- 1
  row_start <- 1
  for (id in dat$id) {
    if (id!=id_prev) {
      row_map[[as.character(id_prev)]] <- c(row_start:(i-1))
      row_start <- i
    }
    id_prev <- id
    i <- i + 1
  }
  row_map[[as.character(id_prev)]] <- c(row_start:(i-1))
  
  # Extract and remove attributes
  n <- attr(dat, "n")
  s_i <- attr(dat, "s_i")
  t_i <- attr(dat, "t_i")
  attr(dat, "n") <- NULL
  attr(dat, "s_i") <- NULL
  attr(dat, "t_i") <- NULL
  
  # Construct a dataframe specific to each individual
  dat_objs <- lapply(c(1:n), function(i) {
    
    d <- list()
    d$dat_i <- dat[row_map[[as.character(i)]],]
    
    # Start and end times
    d$s_i <- s_i[i]
    d$t_i <- t_i[i]
    
    # Create new variables
    d$dat_i$init_visit <- 0
    d$dat_i$j <- d$dat_i$t_end
    d$dat_i[1,"init_visit"] <- 1

    # Apply spline bases to dataframe
    for (s in spl) {
      for (k in c(1:s$df)) {
        d$dat_i[[paste0(s$name,"_",k)]] <- signif(sapply(
          d$dat_i[[s$var]],
          function(x) { do.call(s$name, list(x,k)) }
        ),4)
      }
    }
    
    # Calculate the set X_i to sum over
    d$X_i_set <- list()
    for (j in c(1:(d$t_i-d$s_i+2))) {
      
      x_ <- c(rep(0,d$t_i-d$s_i-j+2), rep(1,j-1))
      T_pm <- T_plusminus(s_i=d$s_i, t_i=d$t_i, x=x_, v=d$dat_i$v)
      case_i_ <- case(T_pm$T_minus, T_pm$T_plus)
      
      if (case_i_!=9 &&
          all(d$dat_i$u==x_*d$dat_i$delta) &&
          all(d$dat_i$delta==g_delta(case=case_i_, s_i=d$s_i, t_i=d$t_i,
                                     T_minus=T_pm$T_minus, T_plus=T_pm$T_plus))
      ) {

        if (length(x_)==1) {
          x_prev <- 0
        } else {
          x_prev <- c(0,x_[1:(length(x_)-1)])
        }
        
        d$X_i_set <- c(d$X_i_set, list(list(
          x = sum(x_),
          x_prev = max(sum(x_)-1,0)
        )))
        
      }
      
    }
    
    return(d)
    
  })
  
  rm(dat)
  
  return(dat_objs)
  
}

#' Negative log-likelihood across individuals and time
#' @param parallelize Whether or not to parallelize computation
#' @return Numeric likelihood
#' @notes This corresponds to the missing data structure
construct_negloglik <- function(
  parallelize=FALSE, use_counter=F, temp=FALSE
) {
  
  if (!identical(temp, FALSE)) { dat_objs_wrapper <- temp }
  
  negloglik <- function(par) {
    
    if (use_counter) {
      counter <<- counter + 1
      if (counter%%10==0) {
        print(paste0(counter, " negloglik calls processed: ",Sys.time()))
      }
    }

    # Compute the negative likelihood across individuals
    if (parallelize) {
      
      nll <- -1 * sum(log(unlist(
        parallel::parLapply(cl, dat_objs_wrapper, function(d) {
          lapply(d, function(d2) { lik_fn(d2, par) })
        })
      )))
      return(nll)
      
    } else {
      
      nll <- -1 * sum(log(unlist(
        lapply(dat_objs_wrapper, function(d) {
          lapply(d, function(d2) { lik_fn(d2, par) })
        })
      )))
      return(nll)
      
    }
    
  }
  
  return(negloglik)
  
}



#' Calculate likelihood for individual i
#'
#' @param d One component of the list dat_objs
#' @param par TO DO
#' @return Numeric likelihood
#' @note dat_objs is accessed globally
lik_fn <- function(d, par) {
  
  {
    
    # Precompute values
    f_xy_vals <- apply(X=d$dat_i, MARGIN=1, FUN = function(r) {
      
      f_x_00 <- f_x(x=0, x_prev=0, s=r[["init_visit"]], par_x=par[par_x],
                    terms_x=terms_x(r=r), par_s=par[par_s], terms_s=terms_s(r=r))
      f_x_01 <- f_x(x=1, x_prev=0, s=r[["init_visit"]], par_x=par[par_x],
                    terms_x=terms_x(r=r), par_s=par[par_s], terms_s=terms_s(r=r))
      f_x_11 <- f_x(x=1, x_prev=1, s=r[["init_visit"]], par_x=par[par_x],
                    terms_x=terms_x(r=r), par_s=par[par_s], terms_s=terms_s(r=r))
      f_y_0 <- f_y(y=r[["y"]], par_y=par[par_y], terms_y=terms_y(r=r, x=0))
      f_y_1 <- f_y(y=r[["y"]], par_y=par[par_y], terms_y=terms_y(r=r, x=1))
      
      return(c(f_x_00*f_y_0, f_x_01*f_y_1, f_x_11*f_y_1))
      
    })
    
    len_obs <- round(d$t_i - d$s_i + 1)
    
    f2_fnc <- function(x_obj) {
      
      x <- uncompress(len_obs, x_obj$x)
      x_prev <- uncompress(len_obs, x_obj$x_prev)
      
      return(prod(sapply(X=c(1:length(x)), FUN = function(j) {
        x_prev_j <- x_prev[j]
        x_j <- x[j]
        if (x_prev_j==0 && x_j==0) {
          return(f_xy_vals[1,j])
        } else if (x_prev_j==0 && x_j==1) {
          return(f_xy_vals[2,j])
        } else if (x_prev_j==1 && x_j==1) {
          return(f_xy_vals[3,j])
        }
      })))
      
    }
    
    f2 <- sum(unlist(lapply(d$X_i_set, f2_fnc)))
      
  }
  
  if (f2<=0) { f2 <- 1e-10 }
  
  return(f2)
  
}



#' Calculate likelihood component f_x
#'
#' @param x Seroconversion indicator (time j)
#' @param x_prev Seroconversion indicator (time j-1)
#' @param s Indicator that this is the first observation for an individual
#' @param par_x Parameter list corresponding to seroconversion model
#' @param terms_x Data values corresponding to par_x
#' @param par_s Parameter list corresponding to initial status model
#' @param terms_s Data values corresponding to par_s
#' @return Numeric likelihood
f_x <- function(x, x_prev, s, par_x, terms_x, par_s, terms_s) {
  if (s==0) {
    if (x==1 && x_prev==1) {
      return(1)
    } else {
      prob <- icll(as.numeric(par_x%*%terms_x))
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  } else {
    prob <- icll(as.numeric(par_s%*%terms_s))
    if (x==1) { return(prob) } else { return(1-prob) }
  }
}



#' Calculate likelihood component f_y
#'
#' @param y Event indicator (time j)
#' @param par_y Parameter list corresponding to outcome model
#' @param terms_y Data values corresponding to par_y
#' @return Numeric likelihood
f_y <- function(y, par_y, terms_y) {
  prob <- icll(as.numeric(par_y%*%terms_y))
  if (y==1) { return(prob) } else { return(1-prob) }
}
