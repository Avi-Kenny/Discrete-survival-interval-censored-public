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
  
  level_set <- level_sets[[cfg$sim_level_set]]
  
}
