# Set this configuration manually
cfg <- list(
  
  # Packages to load
  pkgs = c("dplyr", "survival", "data.table", "tidyr", "memoise", "Rsolnp",
           "numDeriv", "lubridate", "ggplot2", "readstata13", "magrittr"),
  
  # These options are for MAIN.R
  run_sims = F,
  run_analysis = T,
  run_process_results = F,

  # These options are for data analysis
  model_version = 2,
  w_start = 2010,
  w_end = 2022,
  age_end = 60,
  model_sex = Sys.getenv("model_sex"), # This is set in the Slurm `sbatch` call
  samp_size = 0, # 0 means the full dataset is used, otherwise a subsample is used
  save_dat = T,
  
  # These options are for simulations
  # model_version = 1,
  sim_level_set = "level_set_1",
  sim_run_or_update = "run",
  sim_num = 1000,
  # sim_n_cores = 300, # For parallelizing via job arrays
  sim_n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), # For parallelizing across multiple CPUs within a single task
  sim_stop_at_error = F,
  
  # These options are for both sims and analysis
  parallelize = T,
  opt_maxit = 5000,
  opt_r = 4, # 2 for speed, 4 for accuracy
  opt_reltol = 1e-5,
  
  # These options are for process_results.R and process_results2.R
  process_sims = F,
  process_analysis = T,
  w_2 = 0,
  ests_M = "objs/ests_2_full_M_20250202.rds",
  ests_F = "objs/ests_2_full_F_20250202.rds"
  
)

# Change library path
.libPaths(c("/home/akenny/R_lib", "/hpc/home/ak811/R_lib", .libPaths()))

# Set cluster config
if (Sys.getenv("HOME")=="/home/akenny") {
  # Bionic
  cluster_config <- list(
    js = "slurm",
    dir = paste0("/home/akenny/", Sys.getenv("proj"),
                 "/Code__", Sys.getenv("proj"))
  )
} else if (Sys.getenv("HOME")=="/hpc/home/ak811") {
  # DCC
  cluster_config <- list(
    js = "slurm",
    dir = paste0("/hpc/home/ak811/", Sys.getenv("proj"),
                 "/Code__", Sys.getenv("proj"))
  )
} else {
  cluster_config <- list(js="", dir="")
}

# Load packages (if running locally or not running sims)
if (Sys.getenv("RSTUDIO")=="1" || !cfg$run_sims) {
  for (pkg in cfg$pkgs) {
    suppressMessages({ do.call("library", list(pkg)) })
  }
}
