# Set configuration
source("R/config.R", local=T)

# Load SimEngine + functions
{
  library(SimEngine)
  source("R/helpers.R", local=T)
  source("R/one_simulation.R", local=T)
  source("R/generate_data.R", local=T)
  source("R/likelihood.R", local=T)
}

# Run model spec file
source("R/models.R", local=T)

if (cfg$run_analysis) {
  
  # Process data and run analysis
  source("R/process_data.R", local=T)
  source("R/analysis.R", local=T)
  
} else if (cfg$run_sims) {
  
  # Set level sets
  source("R/levels.R", local=T)
  
  # Run simulation
  source("R/run.R", local=T)
  
}

if (cfg$run_process_results) {
  
  # Tables and figures
  source("R/process_results.R", local=T)
  
}
