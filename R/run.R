run_on_cluster(
  
  first = {
    
    # Simulation setup
    sim <- new_sim()
    sim %<>% set_config(
      num_sim = cfg$sim_num,
      parallel = cfg$parallelize,
      n_cores = cfg$sim_n_cores,
      stop_at_error = cfg$sim_stop_at_error,
      packages = cfg$pkgs
    )
    sim <- do.call(set_levels, c(list(sim), level_set))
    
    # Simulation script
    sim %<>% set_script(one_simulation)
    
  },
  
  main = { sim %<>% run() },
  
  last = {},
  
  cluster_config = cluster_config
  
)