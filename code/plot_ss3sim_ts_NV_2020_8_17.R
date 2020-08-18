# load files, packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
dat <- read.csv(file.path("results", "ss3sim_ts_NV_2020_8_17.csv"))

# filter data, calc SSB ratio ----

dat_filter <-  dat %>% 
                separate(scenario, into = c("MA", "Beta", "Index", "Ref_Yr"), 
                         sep = "_", remove = FALSE) %>% 
                filter(MA != 5) %>%  # these scenarios didn't finish so remove. Would want to remove this eventually.
                select(year, iteration, model_run, scenario, SpawnBio, SPRratio)
unique(dat_filter$model_run) #verify OM_error_check not included (b/c would indicate an error)
dat_filter$model_type <- ifelse(grepl("_EM_", dat_filter$model_run), "EM", "OM")
dat_filter_OM <- dat_filter %>% 
                   filter(model_type == "OM") %>% 
                   select(year, iteration, scenario, SpawnBio)
dat_filter_EM <- dat_filter %>% 
                    filter(model_type == "EM") %>% 
                    select(year, iteration, scenario, SpawnBio, SPRratio) %>% 
                    rename(SSB_approach_SPR30 = SpawnBio)

dat_plot <- full_join(dat_filter_OM, dat_filter_EM) %>% 
                    mutate(SSBratio = SpawnBio/SSB_approach_SPR30) %>% 
                    separate(scenario, into = c("MA", "Beta", "Index", "Ref_Yr"), 
                            sep = "_", remove = FALSE)
dat_plot$Index_fac <- factor(dat_plot$Index, levels = c("5", "6", "7", "8", "9"), 
                                labels = c("MRFSS", "HB_E", "Larval", "SEAMAP", "Video"))
dat_plot$Beta_fac <- factor(dat_plot$Beta, levels = c("0", "1","3"),
                               labels = c("Beta = 0", "Beta = 1", "Beta = 3"))
dat_plot$MA_fac <- factor(dat_plot$MA, levels = c("1", "3"))
dat_plot$Ref_Yr_fac <- factor(dat_plot$Ref_Yr, levels = c("0", "2013", "2017", "2035"), 
                                 labels = c("Ref index - pred", "Ref index - yr 2013", 
                                            "Ref index - yr 2017", "Ref index - equil" ))
# calculate RMSE ----
RMSE <- dat_plot %>% 
              mutate(squared_res_SSB = (SpawnBio - SSB_approach_SPR30)^2) %>% 
              filter(year >= 2020) %>% 
              filter(year <= 2025) %>% 
              filter(!is.na(SpawnBio)) %>% 
              select(scenario, squared_res_SSB) %>% 
              group_by(scenario) %>% 
              summarize(mean_squared_SSB = mean(squared_res_SSB)) %>% 
              mutate(RMSE_SSB = sqrt(mean_squared_SSB))

RE <- dat_plot %>% 
        mutate(RE = (SpawnBio - SSB_approach_SPR30)/SpawnBio) %>%
        filter(year >= 2020) %>% 
        filter(year <= 2025) %>% 
        filter(!is.na(SpawnBio)) 
MRE <- RE %>% 
        group_by(scenario) %>% 
        summarize(MRE = median(RE)) %>% 
        separate(scenario, into = c("MA", "Beta", "Index", "Ref_Yr"), 
                 sep = "_", remove = FALSE)
  
# create RE plot ----
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
ggplot(RE, aes(x = Index, y = RE)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(aes(fill = MA), color = "grey50") +
  stat_summary(aes(color = MA), fun.data=data_summary, 
               position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(Beta), cols = vars(Ref_Yr))+
  xlab("Index of Abundance")+
  scale_color_manual(values = rep("black", length.out = length(unique(RE$MA)))) +
  labs(fill = "Moving Avg", 
       title = "RE of SSB compared to SSB@SPR30 in years 2020-2025") +
  scale_fill_brewer(palette = "Set2", direction = -1)+
  guides(color = FALSE)+
  theme()+
  theme_classic()
ggplot2::ggsave(file.path("figures", "RE_2020_2025.png"), height = 12, width = 15, units = "in")     

ggplot(MRE, aes(x = Index, y = MRE)) +
  geom_bar(aes(fill = MA), stat = "identity", position = "dodge" ) +
  facet_grid(rows = vars(Beta), cols = vars(Ref_Yr))+
  theme_classic()
ggplot2::ggsave(file.path("figures", "MedianRE_2020_2025.png"), height = 12, width = 15, units = "in")     
# create plots: ratio of SSB/SSBatSPR30 ----
# create factors for plotting
dat_vioplot <-  dat_plot %>% 
                  filter(!is.na(SpawnBio))
  
SSBratio_2025_plot <-   ggplot(subset(dat_vioplot, year == 2025), aes(x = Index_fac, y = SSBratio)) +
                          geom_hline(yintercept = 1, color = "black") +
                          geom_violin(aes(fill = MA_fac), color = "grey50") +
                          stat_summary(aes(color = MA_fac), fun.data=data_summary, position = position_dodge(width = 0.9)) +
                          facet_grid(rows = vars(Beta_fac), cols = vars(Ref_Yr_fac))+
                          xlab("Index of Abundance")+
                          scale_color_manual(values = rep("black", length.out = length(unique(dat_vioplot$MA)))) +
                          labs(fill = "Moving Avg", 
                               title = "SSB/SSB@SPR30 in year 2025", 
                               subtitle = "Points are medians w/ stdev error bars") +
                          scale_fill_brewer(palette = "Set2", direction = -1)+
                          guides(color = FALSE)+
                          theme()+
                          theme_classic()
ggplot2::ggsave(file.path("figures", "SSBratio_2025.png"), height = 12, width = 15, units = "in")                          

# All years

SSBratio_2020_2025 <- ggplot(subset(dat_vioplot, year >= 2020 & year <= 2025),
                             aes(x = Index_fac, y = SSBratio)) +
  geom_hline(yintercept = 1, color = "black") +
  geom_violin(aes(fill = MA_fac), color = "grey50") +
  stat_summary(aes(color = MA_fac), fun.data=data_summary, 
               position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(Beta_fac), cols = vars(Ref_Yr_fac))+
  xlab("Index of Abundance")+
  scale_color_manual(values = rep("black", length.out = length(unique(dat_vioplot$MA)))) +
  labs(fill = "Moving Avg", 
       title = "SSB/SSB@SPR30, yrs 2020-2025", 
       subtitle = "Points are medians w/ stdev error bars") +
  scale_fill_brewer(palette = "Set2", direction = -1)+
  guides(color = FALSE)+
  theme()+
  theme_classic()
ggplot2::ggsave(file.path("figures", "SSBratio_2020_2025.png"), height = 12, width = 15, units = "in")         

# create plots: SPR30 ----
# Plot SSB patterns. Note that they are the same across iterations, so only need to plot the first iteration
# why is the SPRratio got to 0 in yr 2018? figure this out....Seems to be the same across iterations and 
# scenarios....
plot_vals <-  ggplot(subset(dat_plot, iteration == 1 & scenario == "1_0_5_0"), aes(x = year, y = SPRratio)) +
                geom_line() +
                geom_point() +
                theme_classic()
ggplot2::ggsave(file.path("figures", "SPR_ratio.png"), height = 8, width = 12, units = "in")



