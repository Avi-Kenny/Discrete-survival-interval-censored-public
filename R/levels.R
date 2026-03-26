# Set simulation levels
if (cfg$run_sims && Sys.getenv("sim_run") %in% c("first", "")) {
  
  level_sets <- list()
  
  # Simulation 1: basic
  par_10 <- list(
    a_x = -3,
    g_x = c(0.3,0.2),
    t_x1 = -0.1,
    a_s = -1.6,
    g_s = c(0.5,0.3),
    t_s1 = 0.1,
    beta_x = 0.4,
    a_y = -3.5,
    g_y = c(0.2,0.1),
    t_y = -0.1,
    a_v = -2.4,
    g_v = c(0.2,0.1)
  )
  level_sets[["level_set_1"]] <- list(
    n = 1000,
    max_time = 20,
    model_version = 1,
    par = list("10% testing"=par_10)
  )
  
  # Simulation 2: expanded level sets (24 combinations)
  par_full <- list()
  counter <- 1
  for (a_x in c(-3,-2)) {
    for (beta_x in c(0.2,0.4)) {
      for (t_y in c(-0.1,-0.05)) {
        label <- paste0("a_x=",a_x,", beta_x=",beta_x,", t_y=",t_y)
        par_new <- par_10
        par_new$a_x <- a_x
        par_new$beta_x <- beta_x
        par_new$t_y <- t_y
        par_full[[label]] <- par_new
        counter <- counter + 1
      }
    }
  }
  level_sets[["level_set_2"]] <- list(
    n = 1000,
    max_time = 20,
    model_version = 1,
    par = par_full
  )
  
  level_set <- level_sets[[cfg$sim_level_set]]
  
}
