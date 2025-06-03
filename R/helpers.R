#' Alias for as.integer
#' 
In <- as.integer

#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {1/(1+exp(-x))}



#' Inverse of complementary log-log link function
#' 
#' @param x Numeric input
#' @return Numeric output
icll <- function(x) { 1 - exp(-exp(x)) }



#' Compute "case" indicator variables
#' 
#' @param T_minus Calendar time of first negative test (0 = no NEG tests)
#' @param T_plus Calendar time of first positive test (0 = no POS tests)
#' @return An integer (1-4) representing the case, as follows:
#'     - Case 1: No testing data
#'     - Case 2: Only NEG tests
#'     - Case 3: NEG test then POS test
#'     - Case 4: first test POS
case <- function(T_minus, T_plus) {
  
  if (T_minus==0) {
    if (T_plus==0) { return(1) } else { return(4) }
  } else {
    if (T_plus==0) { return(2) } else { return(3) }
  }
  
}



#' Calculate time(s) of most recent negative test and/or positive test
#' 
#' @param s_i Start time
#' @param t_i End time
#' @param x Vector of seroconversion indicators for one individual
#' @param v Vector of testing indicators for one individual
#' @return A list containing T_minus and T_plus (0 = no NEG/POS tests)
T_plusminus <- function(s_i, t_i, x, v) {
  
  len <- t_i-s_i+1
  if (!length(v)==len || !length(x)==len) { stop("Incorrect vector lengths.") }
  
  some_tests_neg <- In(sum(v*(1-x))>0)
  some_tests_pos <- In(sum(v*x)>0)
  
  if (some_tests_neg) {
    # Time of most recent negative test
    T_minus_index <- which.max(c(s_i:t_i)*v*(1-x))
    T_minus <- c(s_i:t_i)[T_minus_index]
  } else {
    T_minus <- 0
  }
  if (some_tests_pos) {
    # Time of positive test
    T_plus_index <- which.max(v*x)
    T_plus <- c(s_i:t_i)[T_plus_index]
  } else {
    T_plus <- 0
  }
  
  return(list(T_minus=T_minus, T_plus=T_plus))
  
}



#' Compute Delta variable
#' 
#' @param case Case 
#' @param J Maximum time
#' @param x Vector of seroconversion indicators for one individual
#' @param v Vector of testing indicators for one individual
#' @return The vector Delta
g_delta <- function(case, s_i, t_i, T_minus, T_plus) {
  
  if (case==1) {
    d <- rep(0, t_i-s_i+1)
  } else if (case==2) {
    d <- c(rep(1,T_minus-s_i+1), rep(0,t_i-T_minus))
  } else if (case==3) {
    d <- c(rep(1,T_minus-s_i+1), rep(0,T_plus-T_minus-1), rep(1,t_i-T_plus+1))
  } else if (case==4) {
    d <- c(rep(0,T_plus-s_i), rep(1,t_i-T_plus+1))
  } else {
    stop("Error in computing g_delta.")
  }
  if (length(d)!=(t_i-s_i+1)) { stop("Delta is of the wrong length.") }
  
  return(d)
  
}



#' Helper function for debugging; prints timestamps
#' 
#' @param num Number
#' @param msg Message
chk <- function(num, msg="") {
  if (msg=="") {
    str <- paste0("Check ", num, ": ", Sys.time())
  } else {
    str <- paste0("Check ", num, " (", msg, "): ", Sys.time())
  }
  print(str)
}



#' Scale age variable (scaled variable goes into the model)
#' 
#' @param age Age, in completed years
scale_age <- function(age) { age / 100 }



#' Unscale age variable (unscaled variable used for display)
#'
#' @param age Age, in completed years
unscale_age <- function(age) { round(age * 100) }



#' Scale calendar time variable (scaled variable goes into the model)
#' 
#' @param year Year
#' @param st Start year of model
#' @param unit One of c("month", "year"); currently, only "year" implemented
scale_time <- function(year, st, unit="year") { year - st + 1 }



#' Unscale calendar time variable (unscaled variable used for display)
#'
#' @param year Year
#' @param st Start year of model
#' @param unit One of c("month", "year"); currently, only "year" implemented
unscale_time <- function(year, st, unit="year") { year + st - 1 }



#' Create a vector of ones and zeros
#' 
#' @param len Length of vector to output
#' @param num_ones Number of "1"s in the tail
uncompress <- function(len, num_ones) {
  c(rep(0,len-num_ones), rep(1,num_ones))
}



#' Create a natural cubic spline basis
#' 
#' @param which Which basis to construct
#' @param window_start Start year
construct_basis <- function(which, linear=FALSE) {
  
  if (which==c("year (10,13,16,19,22)")) {
    grid <- scale_time(seq(2010,2022, length.out=500), st=2010)
    k <- scale_time(seq(2010,2022, length.out=5), st=2010)
    num_df <- 4
  } else if (which=="age (13,30,40,60)") {
    grid <- scale_age(seq(13,60, length.out=500))
    k <- scale_age(c(13,30,40,60))
    num_df <- 3
  } else if (which=="age (13,40,50,60)") {
    grid <- scale_age(seq(13,60, length.out=500))
    k <- scale_age(c(13,40,50,60))
    num_df <- 3
  } else if (which=="age (13,20,30,40,60)") {
    grid <- scale_age(seq(13,60, length.out=500))
    k <- scale_age(c(13,20,30,40,60))
    num_df <- 4
  }
  
  if (linear) {
    b <- Vectorize(function(x, i) { max(0, x-k[i]) }, vectorize.args="x")
  } else {
    b <- Vectorize(function(x, i) {
      splines::ns(
        x = x,
        knots = k[2:num_df],
        intercept = F,
        Boundary.knots = k[c(1,num_df+1)]
      )[i]
    }, vectorize.args="x")
  }
  y <- matrix(NA, nrow=length(grid), ncol=num_df)
  for (i in c(1:num_df)) { y[,i] <- sapply(grid, function(x) { b(x, i=i) }) }
  rm(b)
  
  return(function(x=NA, i=NA) {
    row <- unlist(lapply(x, function(x) { which.min(abs(x-grid)) } ))
    if (is.na(i)) { return(y[row,]) } else { return(y[row,i]) }
  })
  
}
