#' Generate dataset
#'
#' @param n Number of patients in cohort
#' @param max_time Number of time points in the study
#' @param par List of parameters governing the distribution; see levels.R
#' @return A dataframe, one row per patient
generate_data <- function(n, max_time, par) {
  
  # Generate dataframe to hold results
  dat <- data.frame(
    "id" = integer(),
    "t_end" = integer(),
    "w_1" = integer(),
    "w_2" = double(),
    "x" = integer(),
    "y" = integer(),
    "v" = integer(),
    "d" = integer(),
    "u" = integer()
  )
  
  # Generate baseline covariates
  id <- c(1:n)
  w_1 <- sample(c(0,1), size=n, replace=T)
  w_2 <- sample(c(13:60), size=n, replace=T)
  
  # Scale age variable
  w_2 <- scale_age(w_2)
  
  # Sample start times
  s_i <- sample(c(1:max_time), size=n, replace=T)
  
  # Loop through individuals/time to generate events
  p <- par
  vec_T_minus <- vec_T_plus <- vec_case <- vec_s_i <- vec_t_i <- c()
  for (i in c(1:n)) {
    
    # Initial values
    x <- y <- v <- w_2_vec <- c()
    event <- 0
    known_pos <- known_pos_prev <- 0
    j <- 1
    w_1_ <- w_1[i]
    w_2_ <- w_2[i]
    s_i_ <- s_i[i]
    cal_time <- s_i_
    
    while (!event && j<=max_time) {
      
      # Sample serostatus (x)
      if (j==1) {
        # Sample baseline serostatus
        p_sero <- icll(
          p$a_s + p$g_s[1]*w_1_ + p$g_s[2]*w_2_ + p$t_s*cal_time
        )
      } else {
        p_sero <- x_prev + (1-x_prev) * icll(
          p$a_x + p$g_x[1]*w_1_ + p$g_x[2]*w_2_ + p$t_x*cal_time
        )
      }
      x[j] <- x_prev <- rbinom(n=1, size=1, prob=p_sero)
      
      # Sample testing
      p_test <- (1-known_pos_prev) * icll(
        p$a_v + p$g_v[1]*w_1_ + p$g_v[2]*w_2_
      )
      v[j] <- rbinom(n=1, size=1, prob=p_test)

      # Calculate "known positive" variable
      known_pos[j] <- known_pos_prev + (1-known_pos_prev)*v[j]*x[j]
      known_pos_prev <- known_pos[j]
      
      # Sample events
      p_event <- icll(
        p$a_y + p$g_y[1]*w_1_ + p$g_y[2]*w_2_ + p$t_y*cal_time +
          p$beta_x*x[j]
      )
      
      event <- rbinom(n=1, size=1, prob=p_event)
      y[j] <- event
      
      # Increment time-varying variables
      j <- round(j+1)
      cal_time <- cal_time + 1
      w_2_vec <- c(w_2_vec, w_2_)
      w_2_ <- w_2_ + 0.01
      
    }
    
    # Need to subtract one from j since it is incremented at the end of the loop
    j <- j-1
    
    # End time variable (t_i)
    t_i <- cal_time - 1
    
    # Calculate time(s) of most recent negative test and/or positive test
    T_pm <- T_plusminus(s_i=s_i_, t_i=t_i, x=x, v=v)
    
    # Calculate case indicators
    case_i <- case(T_pm$T_minus, T_pm$T_plus)

    # Calculate Delta
    d <- g_delta(case=case_i, s_i=s_i_, t_i=t_i, T_minus=T_pm$T_minus,
                 T_plus=T_pm$T_plus)
    
    # Add results to dataframe
    dat <- rbind(dat, list(
      id=rep(i,j), t_end=c(s_i_:t_i), w_1=rep(w_1_,j), w_2=w_2_vec,
      y=y, v=v, delta=d, u=x*d # x=x, 
    ))
    
    # Store additional vectors
    vec_T_minus[i] <- T_pm$T_minus
    vec_T_plus[i] <- T_pm$T_plus
    vec_case[i] <- case_i
    vec_s_i[i] <- s_i_
    vec_t_i[i] <- t_i
    
  }
  
  attr(dat, "n") <- n
  attr(dat, "max_time") <- max_time
  attr(dat, "par") <- par
  attr(dat, "T_minus") <- vec_T_minus
  attr(dat, "T_plus") <- vec_T_plus
  attr(dat, "case") <- vec_case
  attr(dat, "s_i") <- vec_s_i
  attr(dat, "t_i") <- vec_t_i
  
  return (dat)
  
}
