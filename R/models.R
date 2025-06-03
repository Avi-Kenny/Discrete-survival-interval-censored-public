# This file sets the following global variables:
#     spl: a list governing which spline bases to apply to the dataset
#     b9, b10, etc: spline bases
#     par_init: initial parameter values (for optimization)
#     par_y: a character vector naming the parameters corresponding to the outcome model; used by f_y() and prob()
#     terms_y: data values corresponding to par_y; used by f_y()
#     terms_y2: a function that generates the "data values" (for predicted probabilities); used by prob()
#     par_x, terms_x, etc.: analogous but for the seroconversion discrete hazard model
#     par_s, terms_s, etc.: analogous but for the initial status model

# Construct spline bases
if (cfg$model_version==2) {
  spl <- list(list(name="b13", var="w_1", df=4),
              list(name="b14", var="t_end", df=4),
              list(name="b15", var="w_1", df=3),
              list(name="b16", var="w_1", df=3))
  b13 <- construct_basis("age (13,20,30,40,60)", linear=F)
  b14 <- construct_basis("year (10,13,16,19,22)", linear=F)
  b15 <- construct_basis("age (13,30,40,60)", linear=F)
  b16 <- construct_basis("age (13,40,50,60)", linear=F)
}

# Set parameter initial values
if (cfg$model_version==1) {
  par_init <- c(a_x=-6.2967, g_x1=-0.1535, g_x2=0.9796, t_x=0.5343, a_s=-2.3111, g_s1=-0.5649, g_s2=0.6198, t_s=0.4245, beta_x=1.401, a_y=-5.5786, g_y1=0.3278, g_y2=4.2046, t_y=-0.7198)
} else if (cfg$model_version==2) {
  if (cfg$model_sex=="Female") { par_init <- c(a_x=-4.81, g_x1=0.13, g_x2=-0.73, g_x3=1.41, g_x4=-2.66, t_x1=-0.14, a_s=-3.11, g_s1=3.67, g_s2=1.94, g_s3=4.27, g_s4=0.37, t_s1=-0.03, beta_x1=0.81, beta_x2=3.52, beta_x3=-2.22, beta_x4=-0.06, beta_x5=-0.06, beta_x6=0.1, beta_x7=0.01, a_y=-7.27, g_y1=0.81, g_y2=0.3, g_y3=3.66, g_y4=2.76, t_y1=-0.25, t_y2=-0.2, t_y3=-0.62, t_y4=-0.08) }
  if (cfg$model_sex=="Male") { par_init <- c(a_x=-8.38, g_x1=-4.78, g_x2=6.45, g_x3=-2.58, t_x1=-0.22, a_s=-3.61, g_s1=3.85, g_s2=2.93, g_s3=3.93, g_s4=1.79, t_s1=-0.04, beta_x1=0.49, beta_x2=4.04, beta_x3=-0.09, beta_x4=-0.08, beta_x5=-0.18, beta_x6=-0.09, beta_x7=0, a_y=-6.94, g_y1=1.51, g_y2=1.23, g_y3=3.78, g_y4=2.62, t_y1=-0.3, t_y2=-0.02, t_y3=-0.33, t_y4=0.03) }
}

# Outcome model
if (cfg$model_version==1) {
  par_y_F <- c("a_y", "t_y", "g_y1", "g_y2", "beta_x")
  terms_y_F <- function(r, x) { c(1, r[["j"]], r[["w_1"]], r[["w_3"]], x) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(1, j, w_1, w_3, x) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version==2) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6", "beta_x7", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x*r[["b16_1"]], x*r[["b16_2"]], x*r[["b16_3"]], x*r[["j"]]*r[["b16_1"]], x*r[["j"]]*r[["b16_2"]], x*r[["j"]]*r[["b16_3"]], x, 1, r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]], r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x*b16(w_1,1), x*b16(w_1,2), x*b16(w_1,3), x*j*b16(w_1,1), x*j*b16(w_1,2), x*j*b16(w_1,3), x, 1, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4), b14(j,1), b14(j,2), b14(j,3), b14(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
}

# Seroconversion model
if (cfg$model_version==1) {
  par_x_F <- c("a_x", "t_x1", "g_x1", "g_x2")
  terms_x_F <- function(r) { c(1, r[["j"]], r[["w_1"]], r[["w_3"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, j, w_1, w_3) }
  par_x_M <- par_x_F; terms_x_M <- terms_x_F; terms_x2_M <- terms_x2_F;
} else if (cfg$model_version==2) {
  par_x_F <- c("a_x", "t_x1", "g_x1", "g_x2", "g_x3", "g_x4")
  terms_x_F <- function(r) { c(1, r[["j"]], r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, j, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4)) }
  par_x_M <- c("a_x", "t_x1", "g_x1", "g_x2", "g_x3")
  terms_x_M <- function(r) { c(1, r[["j"]], r[["b15_1"]], r[["b15_2"]], r[["b15_3"]]) }
  terms_x2_M <- function(j, w_1, w_2) { c(1, j, b15(w_1,1), b15(w_1,2), b15(w_1,3)) }
}

# Initial status model
if (cfg$model_version==1) {
  par_s_F <- c("a_s", "t_s1", "g_s1", "g_s2")
  terms_s_F <- function(r) { c(1, r[["j"]], r[["w_1"]], r[["w_3"]]) }
  terms_s2_F <- function(j, w_1, w_2) { c(1, j, w_1, w_3) }
  par_s_M <- par_s_F; terms_s_M <- terms_s_F; terms_s2_M <- terms_s2_F;
} else if (cfg$model_version==2) {
  par_s_F <- c("a_s", "g_s1", "g_s2", "g_s3", "g_s4", "t_s1")
  terms_s_F <- function(r) { c(1, r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]], r[["j"]]) }
  terms_s2_F <- function(j, w_1, w_2) { c(1, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4), j) }
  par_s_M <- par_s_F; terms_s_M <- terms_s_F; terms_s2_M <- terms_s2_F;
}
