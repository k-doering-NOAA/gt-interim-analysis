# install and load packages ----
devtools::install_github("nmfs-fish-tools/SSMSE", ref = "interim_dev")
library(SSMSE)
library(r4ss)
library(ggplot2) # use install.packages("ggplot2") to install package if needed
library(tidyr) # use install.packages("tidyr") to install package if needed
library(dplyr)

# hard coded values -----
base_mod_path <- file.path("C:/Users/Nathan/Documents/GitHub/gt-interim-analysis","model_files", "GT_3.30_base_run")

#  create folders ----
EM_mod_path <- file.path("C:/Users/Nathan/Documents/GitHub/gt-interim-analysis","SS_runs", "GT_no_var_adj_EM")
dir.create(EM_mod_path,recursive = TRUE)
dir.create("figures")

# read in files datafile ---
# to read in ctl file, and to use in creating the sampling structure
dat <- SS_readdat(file.path(base_mod_path, 
                            "trigger_EonlyShrimp.dat"),
                  verbose = FALSE)
# to self tests, turn off data weighting through input variance and extra SD
ctl <- SS_readctl(file.path(base_mod_path, "Trigger_age_ctl_Eonly.ss"),
                  use_datlist = TRUE,
                  datlist = dat,
                  verbose = FALSE)
# adjust params in the EM model ----
ctl$Variance_adjustment_list$Value <- 1
ctl$Q_parms #not no extra SD so no need to adjusts
#write new model to a different folder.
r4ss::copy_SS_inputs(base_mod_path, EM_mod_path)
SS_writectl(ctl,
            outfile = file.path(EM_mod_path, "Trigger_age_ctl_Eonly.ss"), 
            overwrite = TRUE,
            verbose = FALSE)

# get sampling ---
sampling <- create_sample_struct(dat = dat, nyrs = 3)

# remove NAs for now
sampling$agecomp <- NULL

# modify sampling as needed
sampling_mod <- lapply(sampling, function(x) {
  x <- x[!is.na(x$Yr), ]
  if("SE" %in% colnames(x)) x$SE  <- 0.01
  if("Nsamp" %in% colnames(x))  x$Nsamp  <- 20
  x
}
)
sampling_mod$lencomp[sampling_mod$lencomp$FltSvy == 1, "Part"] <- 2

sample_struct_list <- list()
interim_struct_list <- list()
scenario_name <- vector()
MA_opts <- c(1)#,3,5)
Beta_opts <- c(0)#,1,3)
Index_opts <- c(5)#,6,7,8,9)
Ref_opts <- c(0)#,0,2035)
assess_freq <- 10
for(i in seq_along(MA_opts)){
  MA_years<-MA_opts[i]
  for(j in seq_along(Beta_opts)){
    Beta<-rep(Beta_opts[j],10)
    for(k in seq_along(Index_opts)){
      Index_weights<-rep(0,10)
      Index_weights[Index_opts[k]]<-1
      for(l in seq_along(Ref_opts)){
        Ref_years<-rep(Ref_opts[l],10)
        
        index_val <- (i-1)*length(Beta_opts)*length(Index_opts)*length(Ref_opts) +
                     (j-1)*length(Index_opts)*length(Ref_opts) +
                     (k-1)*length(Ref_opts) + l
        
        interim_struct_list[[index_val]] <- list(MA_years=MA_years, assess_freq=assess_freq,
                                                 Beta=Beta, Index_weights=Index_weights,
                                                 Ref_years=Ref_years)
        
        sample_struct_list[[index_val]] <- sampling_mod
        
        scenario_name[index_val] <- paste(MA_years,Beta,Index_opts[k],Ref_opts[l],sep="_")
      }
    }
  }
}

out_dir <- file.path("C:/Users/Nathan/Documents/GitHub/gt-interim-analysis","SS_runs")
dir.create(out_dir)
out_dir <- file.path("C:/Users/Nathan/Documents/GitHub/gt-interim-analysis","SS_runs", "triggerfish")
dir.create(out_dir)
# run the grey triggerfish base case ----
sampling_list <- run_SSMSE(scen_name_vec = scenario_name,
                           out_dir_scen_vec = out_dir,
                           iter_vec = 5,
                           OM_in_dir_vec = base_mod_path,
                           OM_name_vec = NULL,
                           EM_in_dir_vec = EM_mod_path,
                           EM_name_vec = NULL,
                           MS_vec = "Interim",
                           use_SS_boot_vec = TRUE,
                           scope = 2,
                           rec_dev_pattern = "AutoCorr_rand",
                           nyrs_vec = 5,
                           nyrs_assess_vec = 1,
                           impl_error_pattern = "none",
                           seed = 12345,
                           interim_struct_list = interim_struct_list,
                           sample_struct_list = sample_struct_list)

# plot some results ----
results <- SSMSE_summary_all(out_dir)

# plot SSB by year and model run. Note that forecast years are also plotted.
ggplot2::ggplot(data = results$ts, 
                ggplot2::aes(x = year, y = SpawnBio)) +
  ggplot2::geom_vline(xintercept = 2017, color = "gray") +
  ggplot2::geom_line(ggplot2::aes(linetype = as.character(iteration), color = model_run))+
  #ggplot2::scale_color_manual(values = c("#D65F00", "black")) +
  ggplot2::scale_linetype_manual(values = rep("solid", 10)) +
  ggplot2::guides(linetype = FALSE) +
  #ggplot2::facet_wrap(. ~ scenario) +
  ggplot2::theme_classic()


ggsave(file.path("figures", "SSB_base.png"))
