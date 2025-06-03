# Discrete-survival-interval-censored-public

This repository contains code to reproduce all tables and figures related to "A discrete-time survival model to handle interval-censored covariates, with applications to HIV cohort studies"

Steps to reproduce results:
1) Request access to the dataset "SurveillanceEpisodesHIV.dta" from the [AHRI Data Repository](https://data.ahri.org/index.php/catalog/1183), administered by the Africa Health Research Institute in South Africa.
2) In `R/config`, set configuration variables according to the stage of the workflow you wish to run. In particular, set `run_sims = T` to run the simulation code, set `run_analysis = T` to run the data analysis code, and set `run_process_results = T` to generate the tables and figures. Also, set `model_version = 1` for the simulations and `model_version = 2` for the data analysis.
3) Using a job scheduler on a cluster computing system (such as Slurm), submit file `MAIN.R` as a job, along with environment variable `model_sex=Male` or `model_sex=Female` (if running the data analysis). If running locally, set this environment variable using `Sys.setenv(model_sex="Male")`. If running on a system other than Slurm, it may be necessary to change the configuration variable `sim_n_cores` in `R/config`.
