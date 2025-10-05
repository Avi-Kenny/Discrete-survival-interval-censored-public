#################.
##### Setup #####
#################.

# Get current date
c_date <- format(Sys.time(), "%Y-%m-%d")

# Read in result objects
cfg$ests_M <- readRDS(cfg$ests_M)
cfg$ests_F <- readRDS(cfg$ests_F)



#####################################################.
##### VIZ (simulations): all params (one model) #####
#####################################################.

if (cfg$process_sims) {
  
  sim <- readRDS("SimEngine.out/sim_20240719.rds")
  
  p_set <- "10% testing"
  
  # p_names should match those used by negloglik()
  p_names <- c("a_x", "g_x1", "g_x2", "t_x1",
               "a_s", "g_s1", "g_s2", "t_s1",
               "beta_x",
               "a_y", "g_y1", "g_y2", "t_y")
  p <- sim$levels$par[[p_set]]
  true_vals <- c(p$a_x, p$g_x, p$t_x,
                 p$a_s, p$g_s, p$t_s,
                 p$beta_x,
                 p$a_y, p$g_y, p$t_y)
  names(true_vals) <- p_names
  v <- paste0("lik_M_",p_names,"_est")
  r <- dplyr::filter(sim$results, par==p_set)
  x <- unlist(lapply(v, function(col) { r[,col] }))
  df_true <- data.frame(
    which = factor(v, levels=v),
    val = true_vals
  )
  df_plot <- data.frame(
    x = x,
    y = rep(0, length(x)),
    which = rep(factor(v, levels=v), each=round(length(x)/length(v)))
  )
  
  plot <- ggplot(df_plot, aes(x=x, y=y, color=which)) +
    geom_jitter(width=0, height=1, alpha=0.3, size=3) +
    geom_vline(aes(xintercept=val), df_true, alpha=0.5, linetype="dashed") +
    facet_wrap(~which, ncol=1, strip.position="left") +
    labs(y=NULL, title=p_set) +
    ylim(-2,2) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text.y.left = element_text(angle=0),
          legend.position="none")
  
  ggsave(
    filename=paste0("../Figures + Tables/", c_date, " sim_est_scatterplot.pdf"),
    plot=plot, device="pdf", width=8, height=5
  )
  
  # Summary stats
  summ_mean <- summ_bias <- summ_mean2 <- summ_sd <- summ_cov <- list()
  for (i in c(1:length(p_names))) {
    p <- p_names[i]
    summ_mean[[i]] <- list(stat="mean",
                           name=paste0(p,"__est"),
                           x=paste0("lik_M_",p,"_est"), na.rm=T)
    summ_bias[[i]] <- list(stat="bias",
                           name=paste0(p,"__bias"),
                           estimate=paste0("lik_M_",p,"_est"),
                           truth=true_vals[[p]])
    summ_mean2[[i]] <- list(stat="mean",
                           name=paste0(p,"__sd_est"),
                           x=paste0("lik_M_",p,"_se"), na.rm=T)
    summ_sd[[i]] <- list(stat="sd", name=paste0(p,"__sd_actual"),
                         x=paste0("lik_M_",p,"_est"), na.rm=T)
    summ_cov[[i]] <- list(stat="coverage", name=paste0(p,"__cov"),
                          truth=true_vals[i],
                          estimate=paste0("lik_M_",p,"_est"),
                          se=paste0("lik_M_",p,"_se"), na.rm=T)
  }
  summ <- do.call(
    SimEngine::summarize,
    c(list(sim), summ_mean, summ_bias, summ_mean2, summ_sd, summ_cov)
  )
  l_id <- 1
  summ2 <- summ[summ$level_id==l_id]
  df_results <- data.frame(
    "par" = character(),
    "truth" = double(),
    "est" = double(),
    "bias" = double(),
    "sd_est" = double(),
    "sd_actual" = double(),
    "coverage" = double()
  )
  for (p in p_names) {
    df_results[nrow(df_results)+1,] <- c(
      p,
      true_vals[[p]],
      round(summ[[paste0(p,"__est")]], 3),
      round(summ[[paste0(p,"__bias")]], 3),
      round(summ[[paste0(p,"__sd_est")]], 3),
      round(summ[[paste0(p,"__sd_actual")]], 3),
      round(summ[[paste0(p,"__cov")]], 3)
    )
  }
  utils::write.table(
    x = df_results,
    file = paste0("../Figures + Tables/", c_date, " sims_sd_and_coverage.csv"),
    sep = ",",
    row.names = FALSE
  )
  
}



########################################################.
##### Functions for plotting modeled probabilities #####
########################################################.

#' Return modeled probability
#' @param type One of c("sero", "init", "mort (HIV-neg)", "mort (HIV-pos)")
#' @param j Calendar time, unscaled (e.g., 2010)
#' @param w_1 Age, unscaled (in completed years)
#' @param w_2 Geography (0="PIPSA South", 1="PIPSA North")
#' @param w_3 Sex (0=female, 1=male)
#' @param m An integer representing the model version number
#' @param which One of c("est", "ci_lo", "ci_up")
#' @return Numeric probability
prob <- function(type, j, w_1, w_2, w_3, year_start, which="est") {
  
  # Scale time and age
  j <- scale_time(j, st=year_start, unit="year")
  w_1 <- scale_age(w_1)
  
  # Set ests objcet (sex-specific)
  if (w_3==1) {
    ests <- cfg$ests_M
  } else {
    ests <- cfg$ests_F
  }
  
  # Set params and terms (sex-specific)
  if (type=="sero") {
    if (w_3==1) {
      p2 <- par_x_M
      A <- t(matrix(terms_x2_M(j, w_1, w_2)))
    } else {
      p2 <- par_x_F
      A <- t(matrix(terms_x2_F(j, w_1, w_2)))
    }
  } else if (type=="init") {
    if (w_3==1) {
      p2 <- par_s_M
      A <- t(matrix(terms_s2_M(j, w_1, w_2)))
    } else {
      p2 <- par_s_F
      A <- t(matrix(terms_s2_F(j, w_1, w_2)))
    }
  } else {
    if (type=="mort (HIV-neg)") { x <- 0 }
    if (type=="mort (HIV-pos)") { x <- 1 }
    if (w_3==1) {
      p2 <- par_y_M
      A <- t(matrix(terms_y2_M(x, j, w_1, w_2)))
    } else {
      p2 <- par_y_F
      A <- t(matrix(terms_y2_F(x, j, w_1, w_2)))
    }
  }
  
  indices <- as.numeric(sapply(p2, function(p) {
    which(names(ests$opt$par)==p)
  }))
  beta <- matrix(ests$opt$par[indices])
  Sigma <- ests$hessian_inv[indices,indices]
  if (which=="est") {
    fac <- 0
  } else if (which=="ci_lo") {
    fac <- -1.96
  } else if (which=="ci_up") {
    fac <- 1.96
  }
  est <- c(A %*% beta)
  se <- c(sqrt(A %*% Sigma %*% t(A)))
  prob <- icll(est+fac*se)

  return(prob)
  
}



#' Return plot of modeled mortality probabilities (HIV-neg vs. HIV-pos) with CIs
#' @param x_axis One of c("Year", "Age"); the variable to go on the X-axis
#' @param m An integer representing the model version number
#' @param w_start An integer representing the window start calendar year
#' @param w_start An integer representing the window end calendar year
#' @param y_max Maximum Y value for the plot
#' @param title Boolean; if F, title is suppressed
#' @return ggplot2 object
plot_mort3 <- function(x_axis, w_start, w_end, y_max=NA, title=T) {
  
  if (x_axis=="Age") {
    
    grid <- seq(13,60,0.1)
    breaks <- seq(20,60, length.out=5)
    color <- "Year"
    if (w_start==2000) {
      outer <- c(2000,2010,2020)
    } else if (w_start==2010) {
      outer <- c(2012,2016,2020)
    } else if (w_start==2017) {
      outer <- c(2017,2019,2022)
    }
    prob2 <- function(type, which, outer) {
      1000 * sapply(grid, function(inner) {
        prob(type=type, j=outer, w_1=inner, w_2=cfg$w_2, w_3=sex,
             year_start=w_start, which=which)
      })
    }
    
  } else if (x_axis=="Year") {
    
    grid <- seq(w_start,w_end,0.01)
    breaks <- seq(w_start,w_end, length.out=5)
    color <- "Age"
    outer <- c(20,35,50)
    prob2 <- function(type, which, outer) {
      1000 * sapply(grid, function(inner) {
        prob(type=type, j=inner, w_1=outer, w_2=cfg$w_2, w_3=sex,
             year_start=w_start, which=which)
      })
    }
    
  }
  
  init <- T
  
  for (o in outer) {
    for (sex in c(0,1)) {
      
      plot_data <- data.frame(
        x = rep(grid,2),
        Rate = c(prob2(type="mort (HIV-neg)", which="est", outer=o),
                 prob2(type="mort (HIV-pos)", which="est", outer=o)),
        ci_lo = c(prob2(type="mort (HIV-neg)", which="ci_lo", outer=o),
                  prob2(type="mort (HIV-pos)", which="ci_lo", outer=o)),
        ci_up = c(prob2(type="mort (HIV-neg)", which="ci_up", outer=o),
                  prob2(type="mort (HIV-pos)", which="ci_up", outer=o)),
        color = rep(c("HIV-negative","HIV-positive"), each=length(grid)),
        outer = o,
        sex = ifelse(sex, "Male", "Female")
      )
      
      if (x_axis=="Year") {
        plot_data$outer <- paste0("Age: ", plot_data$outer)
      }
      
      if (init) {
        plot_data2 <- plot_data
        init <- F
      } else {
        plot_data2 <- rbind(plot_data2, plot_data)
      }
      
    }
  }
  
  if (title) {
    title <- "Conditional mortality rates by HIV status, over calendar time"
  } else {
    title <- NULL
  }
  
  plot <- ggplot(
    plot_data2,
    aes(x=x, y=Rate, color=color, linetype=color)
  ) +
    geom_line() +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_up, fill=color),
      alpha = 0.2,
      linetype = "blank"
    ) +
    facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(outer)) +
    coord_cartesian(ylim=c(0,y_max)) +
    labs(
      title = title,
      x = x_axis,
      color = "HIV Status",
      linetype = "HIV Status",
      fill = "HIV Status",
      y = "Deaths per 1,000 person-years"
    ) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=breaks) +
    scale_color_manual(values=c("forestgreen", "#56B4E9")) +
    scale_linetype_manual(values=c("longdash", "solid")) +
    scale_fill_manual(values=c("forestgreen", "#56B4E9"))
  
  return(plot)
  
}



#' Return plot of modeled mortality probabilities (HIV-neg vs. HIV-pos) with CIs
#' @param m An integer representing the model version number
#' @param w_start An integer representing the window start calendar year
#' @param y_max Maximum Y value for the plot
#' @param title Boolean; if F, title is suppressed
#' @return ggplot2 object
plot_sero3 <- function(type, w_start, y_max=NA, title=T) {
  
  grid <- seq(13,60,0.1)
  breaks <- seq(20,60, length.out=5)
  color <- "Year"
  if (w_start==2000) {
    outer <- c(2000,2010,2020)
  } else if (w_start==2010) {
    outer <- c(2012,2016,2020)
  } else if (w_start==2017) {
    outer <- c(2017,2019,2022)
  }
  prob2 <- function(which, outer) {
    sapply(grid, function(inner) {
      prob(type=type, j=outer, w_1=inner, w_2=cfg$w_2, w_3=sex,
           year_start=w_start, which=which)
    })
  }
  
  init <- T
  
  for (o in outer) {
    for (sex in c(0,1)) {
      
      plot_data <- data.frame(
        x = rep(grid,2),
        Rate = prob2(which="est", outer=o),
        ci_lo = prob2(which="ci_lo", outer=o),
        ci_up = prob2(which="ci_up", outer=o),
        outer = o,
        sex = ifelse(sex, "Male", "Female")
      )
      
      if (init) {
        plot_data2 <- plot_data
        init <- F
      } else {
        plot_data2 <- rbind(plot_data2, plot_data)
      }
      
    }
  }
  
  if (type=="sero") {
    plot_title <- "Probability of seroconversion (in one year), by age"
    y <- "Probability of seroconversion (in one year)"
  } else if (type=="init") {
    plot_title <- "Probability that initial status is HIV+, by age"
    y <- "Probability that initial status is HIV+"
  }
  if (!title) { plot_title <- NULL }
  
  plot <- ggplot(
    plot_data2,
    aes(x=x, y=Rate)
  ) +
    geom_line(color="forestgreen") +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_up),
      alpha = 0.2,
      fill = "forestgreen"
    ) +
    facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(outer)) +
    coord_cartesian(ylim=c(0,y_max)) +
    labs(
      title = plot_title,
      x = "Age",
      y = y
    ) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=breaks)

  return(plot)
  
}



###############################################.
##### VIZ: Plotting modeled probabilities #####
###############################################.

if (cfg$process_analysis) {
  
  # Mortality by calendar time, with CIs
  plot_05 <- plot_mort3(x_axis="Year", w_start=cfg$w_start, w_end=cfg$w_end,
                        y_max=55, title=F)
  ggsave(
    filename = paste0("../Figures + Tables/", c_date, " p5 (mort CIs, by year)",
                      ".pdf"),
    plot = plot_05, device="pdf", width=8, height=5
  )
  
  # Mortality by age, with CIs
  plot_06 <- plot_mort3(x_axis="Age", w_start=cfg$w_start, w_end=cfg$w_end,
                        y_max=75, title=F)
  ggsave(
    filename = paste0("../Figures + Tables/", c_date, " p6 (mort CIs, by age) ",
                      ".pdf"),
    plot = plot_06, device="pdf", width=8, height=5
  )
  
  # Seroconversion by age, with CIs
  plot_07 <- plot_sero3(type="sero", w_start=cfg$w_start, y_max=0.05, title=F)
  ggsave(
    filename = paste0("../Figures + Tables/", c_date, " p7 (sero CIs, by age) ",
                      ".pdf"),
    plot = plot_07, device="pdf", width=8, height=5
  )
  
  # Initial status by age, with CIs
  plot_08 <- plot_sero3(type="init", w_start=cfg$w_start, y_max=0.75, title=F)
  ggsave(
    filename = paste0("../Figures + Tables/", c_date, " p8 (init CIs, by age) ",
                      ".pdf"),
    plot = plot_08, device="pdf", width=8, height=5
  )
  
}



#############################################.
##### VIZ: HRs as functions of age/time #####
#############################################.

if (cfg$process_analysis) {
  
  # Functions to return spline basis as matrix
  A_b9 <- function(w_1) {
    t(matrix(c(b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4))))
  }
  A_b10 <- function(j) {
    t(matrix(c(b10(j,1), b10(j,2), b10(j,3), b10(j,4))))
  }
  A_b12 <- function(j) {
    t(matrix(c(b12(j,1), b12(j,2), b12(j,3), b12(j,4))))
  }
  A_b13 <- function(w_1) {
    t(matrix(c(b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4))))
  }
  A_b14 <- function(j) {
    t(matrix(c(b14(j,1), b14(j,2), b14(j,3), b14(j,4))))
  }
  A_b15 <- function(w_1) {
    t(matrix(c(b15(w_1,1), b15(w_1,2), b15(w_1,3))))
  }
  A_b16 <- function(w_1) {
    t(matrix(c(b16(w_1,1), b16(w_1,2), b16(w_1,3))))
  }
  
  # Extract estimates and SEs
  plot_names <- c(
    "HR_mort_hiv_cal",
    "HR_mort_age",
    "HR_mort_cal",
    "HR_sero_cal",
    "HR_sero_age",
    "HR_init_age"
  )
  for (plot_name in plot_names) {
    
    # Set graph-specific variables
    if (plot_name=="HR_mort_hiv_cal") {
      
      title <- "HR of mortality, HIV-positive vs. HIV-negative individuals"
      par <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6", "beta_x7")
      A <- function(j, w_1) { t(matrix(c(
        b16(w_1,1), b16(w_1,2), b16(w_1,3),
        j*b16(w_1,1), j*b16(w_1,2), j*b16(w_1,3), 1
      ))) }
      par_F <- par_M <- par; A_M <- A_F <- A;
      
    } else if (plot_name=="HR_mort_age") {
      
      title <- "HR of mortality (age)"
      x_axis <- "age"
      par <- c("g_y1", "g_y2", "g_y3", "g_y4")
      A <- A_b13
      par_F <- par_M <- par; A_M <- A_F <- A;
      
    } else if (plot_name=="HR_mort_cal") {
      
      title <- "HR of mortality (calendar time)"
      x_axis <- "cal time"
      par <- c("t_y1", "t_y2", "t_y3", "t_y4")
      A <- A_b14
      par_F <- par_M <- par; A_M <- A_F <- A;
      
    } else if (plot_name=="HR_sero_cal") {
      
      title <- "HR of seroconversion (calendar time)"
      par <- c("t_x1")
      A <- function(j) { matrix(j) }
      par_F <- par_M <- par; A_M <- A_F <- A;
      x_axis <- "cal time"
      
    } else if (plot_name=="HR_sero_age") {
      
      title <- "HR of seroconversion (age)"
      x_axis <- "age"
      par_F <- c("g_x1", "g_x2", "g_x3", "g_x4")
      A_F <- A_b13
      par_M <- c("g_x1", "g_x2", "g_x3")
      A_M <- A_b15
      
    } else if (plot_name=="HR_init_age") {
      
      title <- "HR of HIV+ initial status (age)"
      x_axis <- "age"
      par <- c("g_s1", "g_s2", "g_s3", "g_s4")
      A <- A_b13
      par_F <- par_M <- par; A_M <- A_F <- A;
      
    }
    
    indices_M <- which(names(cfg$ests_M$opt$par) %in% par_M)
    beta_M <- matrix(cfg$ests_M$opt$par[indices_M])
    Sigma_M <- cfg$ests_M$hessian_inv[indices_M,indices_M]
    indices_F <- which(names(cfg$ests_F$opt$par) %in% par_F)
    beta_F <- matrix(cfg$ests_F$opt$par[indices_F])
    Sigma_F <- cfg$ests_F$hessian_inv[indices_F,indices_F]
    
    hr <- function(x, which, sex, log=F) {
      if (which=="cal time") {
        x <- scale_time(x, st=cfg$w_start)
      } else if (which=="age") {
        x <- scale_age(x)
      }
      if (sex=="M") {
        est <- c(A_M(x) %*% beta_M)
        se <- c(sqrt(A_M(x) %*% Sigma_M %*% t(A_M(x))))
      } else if (sex=="F") {
        est <- c(A_F(x) %*% beta_F)
        se <- c(sqrt(A_F(x) %*% Sigma_F %*% t(A_F(x))))
      }
      if (log) {
        return(est + c(0,-1.96,1.96)*se)
      } else {
        return(exp(est + c(0,-1.96,1.96)*se))
      }
    }
    
    hr2 <- function(j, w_1, sex, log=F) {
      j <- scale_time(j, st=cfg$w_start)
      w_1 <- scale_age(w_1)
      if (sex=="M") {
        est <- c(A_M(j, w_1) %*% beta_M)
        se <- c(sqrt(A_M(j, w_1) %*% Sigma_M %*% t(A_M(j, w_1))))
      } else if (sex=="F") {
        est <- c(A_F(j, w_1) %*% beta_F)
        se <- c(sqrt(A_F(j, w_1) %*% Sigma_F %*% t(A_F(j, w_1))))
      }
      if (log) {
        return(est + c(0,-1.96,1.96)*se)
      } else {
        return(exp(est + c(0,-1.96,1.96)*se))
      }
    }
    
    if (plot_name=="HR_mort_hiv_cal") {
      
      log <- F
      
      # First, make graph with x_axis=="cal time"
      grid <- seq(cfg$w_start,cfg$w_end,0.1)
      breaks <- seq(cfg$w_start, cfg$w_end, length.out=5)
      df_plot <- data.frame(
        x = rep(grid,6),
        y = c(
          sapply(grid, function(j) { hr2(j=j, w_1=20, sex="M", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=20, sex="F", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=35, sex="M", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=35, sex="F", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=50, sex="M", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=50, sex="F", log=log)[1] })
        ),
        ci_lo = c(
          sapply(grid, function(j) { hr2(j=j, w_1=20, sex="M", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=20, sex="F", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=35, sex="M", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=35, sex="F", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=50, sex="M", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=50, sex="F", log=log)[2] })
        ),
        ci_up = c(
          sapply(grid, function(j) { hr2(j=j, w_1=20, sex="M", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=20, sex="F", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=35, sex="M", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=35, sex="F", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=50, sex="M", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=50, sex="F", log=log)[3] })
        ),
        sex = rep(rep(c("Male", "Female"),3), each=length(grid)),
        age = rep(
          rep(c("Age: 20", "Age: 35", "Age: 50"), each=2),
          each=length(grid)
        )
      )
      plot <- ggplot(
        data = df_plot,
        aes(x=x, y=y)) +
        geom_line(color="forestgreen") +
        geom_ribbon(
          aes(ymin=ci_lo, ymax=ci_up),
          alpha = 0.2,
          fill = "forestgreen"
        ) +
        scale_x_continuous(
          breaks = breaks,
          minor_breaks = seq(cfg$w_start, cfg$w_end, 1)
        ) +
        scale_y_continuous(
          breaks = c(1:16),
          minor_breaks = NULL
        ) +
        coord_cartesian(ylim=c(1,16)) +
        facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(age)) +
        labs(
          x = "Year",
          y = "Hazard Ratio",
          title = title
        )
      ggsave(
        filename = paste0("../Figures + Tables/", c_date, " ", plot_name,
                          " (by year).pdf"),
        plot=plot, device="pdf", width=9, height=6
      )
      plot_paper <- plot + labs(title=NULL)
      ggsave(
        filename = paste0("../Figures + Tables/", c_date, " paper_", plot_name,
                          " (by year).pdf"),
        plot=plot_paper, device="pdf", width=9, height=6
      )
      
      # Second, make graph with x_axis=="age"
      grid <- seq(13,60,0.1)
      breaks <- seq(20,60, length.out=5)
      df_plot <- data.frame(
        x = rep(grid,6),
        y = c(
          sapply(grid, function(w_1) { hr2(j=2012, w_1=w_1, sex="M", log=log)[1] }),
          sapply(grid, function(w_1) { hr2(j=2012, w_1=w_1, sex="F", log=log)[1] }),
          sapply(grid, function(w_1) { hr2(j=2016, w_1=w_1, sex="M", log=log)[1] }),
          sapply(grid, function(w_1) { hr2(j=2016, w_1=w_1, sex="F", log=log)[1] }),
          sapply(grid, function(w_1) { hr2(j=2020, w_1=w_1, sex="M", log=log)[1] }),
          sapply(grid, function(w_1) { hr2(j=2020, w_1=w_1, sex="F", log=log)[1] })
        ),
        ci_lo = c(
          sapply(grid, function(w_1) { hr2(j=2012, w_1=w_1, sex="M", log=log)[2] }),
          sapply(grid, function(w_1) { hr2(j=2012, w_1=w_1, sex="F", log=log)[2] }),
          sapply(grid, function(w_1) { hr2(j=2016, w_1=w_1, sex="M", log=log)[2] }),
          sapply(grid, function(w_1) { hr2(j=2016, w_1=w_1, sex="F", log=log)[2] }),
          sapply(grid, function(w_1) { hr2(j=2020, w_1=w_1, sex="M", log=log)[2] }),
          sapply(grid, function(w_1) { hr2(j=2020, w_1=w_1, sex="F", log=log)[2] })
        ),
        ci_up = c(
          sapply(grid, function(w_1) { hr2(j=2012, w_1=w_1, sex="M", log=log)[3] }),
          sapply(grid, function(w_1) { hr2(j=2012, w_1=w_1, sex="F", log=log)[3] }),
          sapply(grid, function(w_1) { hr2(j=2016, w_1=w_1, sex="M", log=log)[3] }),
          sapply(grid, function(w_1) { hr2(j=2016, w_1=w_1, sex="F", log=log)[3] }),
          sapply(grid, function(w_1) { hr2(j=2020, w_1=w_1, sex="M", log=log)[3] }),
          sapply(grid, function(w_1) { hr2(j=2020, w_1=w_1, sex="F", log=log)[3] })
        ),
        sex = rep(rep(c("Male", "Female"),3), each=length(grid)),
        age = rep(
          rep(c("2012", "2016", "2020"), each=2),
          each=length(grid)
        )
      )
      plot <- ggplot(
        data = df_plot,
        aes(x=x, y=y)) +
        geom_line(color="forestgreen") +
        geom_ribbon(
          aes(ymin=ci_lo, ymax=ci_up),
          alpha = 0.2,
          fill = "forestgreen"
        ) +
        scale_x_continuous(
          breaks = breaks,
          minor_breaks = seq(cfg$w_start, cfg$w_end, 1)
        ) +
        scale_y_continuous(
          breaks = c(1:16),
          minor_breaks = NULL
        ) +
        coord_cartesian(ylim=c(1,16)) +
        facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(age)) +
        labs(
          x = "Age",
          y = "Hazard Ratio",
          title = title
        )
      ggsave(
        filename = paste0("../Figures + Tables/", c_date, " ", plot_name,
                          " (by age).pdf"),
        plot=plot, device="pdf", width=9, height=6
      )
      plot_paper <- plot + labs(title=NULL)
      ggsave(
        filename = paste0("../Figures + Tables/", c_date, " paper_", plot_name,
                          " (by age).pdf"),
        plot=plot_paper, device="pdf", width=9, height=6
      )
      
    } else {
      
      if (x_axis=="cal time") {
        grid <- seq(cfg$w_start,cfg$w_end,0.1)
        breaks <- seq(cfg$w_start, cfg$w_end, length.out=5)
      } else if (x_axis=="age") {
        grid <- seq(13,60,0.1)
        breaks <- seq(20,60, length.out=5)
      }
      
      df_plot <- data.frame(
        x = rep(grid,4),
        y = c(
          sapply(grid, function(x) { hr(x, x_axis, sex="M", log=T)[1] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="M", log=F)[1] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="F", log=T)[1] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="F", log=F)[1] })
        ),
        ci_lo = c(
          sapply(grid, function(x) { hr(x, x_axis, sex="M", log=T)[2] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="M", log=F)[2] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="F", log=T)[2] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="F", log=F)[2] })
        ),
        ci_up = c(
          sapply(grid, function(x) { hr(x, x_axis, sex="M", log=T)[3] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="M", log=F)[3] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="F", log=T)[3] }),
          sapply(grid, function(x) { hr(x, x_axis, sex="F", log=F)[3] })
        ),
        sex = rep(c("Male", "Female"), each=2*length(grid)),
        log = rep(c("log(HR)", "HR", "log(HR)", "HR"), each=length(grid))
      )
      
      plot <- ggplot(
        data = df_plot,
        aes(x=x, y=y)) +
        geom_line(color="forestgreen") +
        geom_ribbon(
          aes(ymin=ci_lo, ymax=ci_up),
          alpha = 0.2,
          fill = "forestgreen"
        ) +
        scale_x_continuous(breaks=breaks) +
        facet_grid(rows=dplyr::vars(log), cols=dplyr::vars(sex), scales="free") +
        labs(
          x = ifelse(x_axis=="cal time", "Year", "Age"),
          y = "Hazard Ratio",
          title = title
        )

      ggsave(
        filename = paste0("../Figures + Tables/", c_date, " ", plot_name,
                          ".pdf"),
        plot=plot, device="pdf", width=8, height=4
      )
      
    }
    
  }
  
}
