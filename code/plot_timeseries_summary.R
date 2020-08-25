# load files, packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
dat1 <- read.csv(file.path("results", "SSMSE_ts_part1.csv"))
dat2 <- read.csv(file.path("results", "SSMSE_ts_part2.csv"))

# preprocess data ----
# reassign scenario numbers based on if there already iterations with the same
# scenario numbers in dat1.
scens <- unique(c(dat1$scenario,dat2$scenario))
for (i in scens) {
  dat_1_scens_iters <- as.integer(unique(dat1[dat1$scenario == i, "iteration"]))
  dat_2_scens_iters <- as.integer(unique(dat2[dat2$scenario == i, "iteration"]))
  if(length(dat_1_scens_iters) == 0 | length(dat_2_scens_iters) == 0) next
  if(min(dat_2_scens_iters) < max(dat_1_scens_iters)) {
     dat2$iteration <- as.integer(dat2$iteration) + max(dat_1_scens_iters)
  }
}
dat <- rbind(dat1, dat2)

# filter data, calc SSB ratio ----

dat_filter <-  dat %>% 
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
                    mutate(SSBratio = SpawnBio/SSB_approach_SPR30)
dat_ctl <- dat_plot %>% 
  filter(scenario == "control")
dat_ctl$MA <- "1"
dat_ctl$Beta <- "control"
dat_ctl$Index <- "5"
dat_ctl$Ref_Yr <- "0"
dat_no_ctl <- dat_plot %>% 
                filter( scenario != "control") %>% 
                separate(scenario, into = c("MA", "Beta", "Index", "Ref_Yr"), 
                         sep = "_", remove = FALSE)
dat_ctl <- dat_ctl[,colnames(dat_no_ctl)]
dat_plot <- rbind(dat_no_ctl, dat_ctl)

dat_plot$Index_fac <- factor(dat_plot$Index, levels = c("5", "6", "7", "8", "9"), 
                                labels = c("MRFSS", "HB_E", "Larval", "SEAMAP", "Video"))
dat_plot$Beta_fac <- factor(dat_plot$Beta, levels = c("control", "0", "1","3"),
                               labels = c("control", "0", "1", "3"))
dat_plot$MA_fac <- factor(dat_plot$MA, levels = c("1", "3", "5"))
dat_plot$Ref_Yr_fac <- factor(dat_plot$Ref_Yr, levels = c("0", "2013", "2017", "2035"), 
                                 labels = c("Ref index - pred", "Ref index - yr 2013", 
                                            "Ref index - yr 2017", "Ref index - equil" ))
# create plots: ratio of SSB/SSBatSPR30 ----
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# create factors for plotting
# plot where facet_grid x is the survey index
dat_vioplot <- dat_plot %>% 
                 filter(!is.na(SpawnBio))
SSBratio_2020_2025_by_index <- ggplot(subset(dat_vioplot, year >= 2020 & year <= 2025),
                             aes(x = MA_fac, y = SSBratio)) +
  geom_hline(yintercept = 1, color = "black") +
  geom_violin(aes(fill = Beta_fac), color = "grey50") +
  stat_summary(aes(color = Beta_fac), fun.data=data_summary, 
               position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(Ref_Yr_fac), cols = vars(Index_fac))+
  xlab("Moving average")+
  scale_color_manual(values = rep("black", length.out = length(unique(dat_vioplot$Beta_fac)))) +
  labs(fill = "Beta", 
       title = "SSB/SSB@SPR30, yrs 2020-2025", 
       subtitle = "Points are medians w/ stdev error bars") +
  scale_fill_brewer(palette = "Set2", direction = -1)+
  guides(color = FALSE)+
  theme()+
  theme_classic()
ggplot2::ggsave(file.path("figures", "SSBratio_2020_2025_violins.png"), height = 12, width = 15, units = "in")

data_summary_CI <- function(x) {
  m <- median(x)
  ymin <- unname(quantile(x, probs = 0.025))
  ymax <- unname(quantile(x, probs = 0.975))
  return(c(y = m,ymin = ymin,ymax = ymax))
}

SSBratio_2020_2025_error_bars <- ggplot(subset(dat_vioplot, year >= 2020 & year <= 2025),
                                      aes(x = MA_fac, y = SSBratio)) +
  geom_hline(yintercept = 1, color = "black") +
  stat_summary(aes(color = Beta_fac), fun.data = data_summary_CI, 
               position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(Ref_Yr_fac), cols = vars(Index_fac))+
  xlab("Moving average")+
  #scale_color_manual(values = rep("black", length.out = length(unique(dat_vioplot$Beta_fac)))) +
  labs(color = "Beta", 
       title = "SSB/SSB@SPR30, yrs 2020-2025", 
       subtitle = "Points are medians w/ 95% CI") +
  scale_color_brewer(palette = "Set2", direction = -1)+
 #guides(color = TRUE)+
  theme()+
  theme_classic()
ggplot2::ggsave(file.path("figures", "SSBratio_2020_2025_CIs.png"), height = 12, width = 15, units = "in")

# calculate SSB/SSB_ctl ----
head(dat_plot)
dat_ctl <- dat_plot %>% 
                  filter(scenario == "control") %>% 
                  group_by(year) %>% 
                  summarize(SpawnBio_ctl_mean = mean(SpawnBio))


dat_no_ctl <- dat_plot %>% 
                filter(scenario != "control")
SSB_rel_ctl_plot_dat <- full_join(dat_no_ctl, dat_ctl, by = "year") %>% 
                          mutate(SSB_relative_ctl = SpawnBio/SpawnBio_ctl_mean) %>% 
                          filter(!is.na(SpawnBio)) %>% 
                          filter(year >= 2020) %>% 
                          filter(year <= 2025)
# plot SSB/SSB_ctl ----

SSBratio_ctl_2020_2025 <- ggplot(SSB_rel_ctl_plot_dat,
                                      aes(x = MA_fac, y = SSB_relative_ctl)) +
  geom_hline(yintercept = 1, color = "black") +
  geom_violin(aes(fill = Beta_fac), color = "grey50") +
  stat_summary(aes(color = Beta_fac), fun.data=data_summary, 
               position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(Ref_Yr_fac), cols = vars(Index_fac))+
  xlab("Moving average")+
  ylab("SSB/SSB_control")+
  scale_color_manual(values = rep("black", length.out = length(unique(dat_vioplot$Beta_fac)))) +
  labs(fill = "Beta", 
       title = "SSB/SSB_control, yrs 2020-2025", 
       subtitle = "Points are medians w/ std dev error bars") +
  scale_fill_brewer(palette = "Set2", direction = -1)+
  guides(color = FALSE)+
  theme()+
  theme_classic()
ggplot2::ggsave(file.path("figures", "SSBratio_ctl_2020_2025_violins.png"), height = 12, width = 15, units = "in")

SSBratio_ctl_2020_2025_error_bars <- ggplot(SSB_rel_ctl_plot_dat,
                                        aes(x = MA_fac, y = SSB_relative_ctl)) +
  geom_hline(yintercept = 1, color = "black") +
  stat_summary(aes(color = Beta_fac), fun.data = data_summary_CI, 
               position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(Ref_Yr_fac), cols = vars(Index_fac))+
  xlab("Moving average")+
  ylab("SSB/SSB_control")+
  #scale_color_manual(values = rep("black", length.out = length(unique(dat_vioplot$Beta_fac)))) +
  labs(color = "Beta", 
       title = "SSB/SSB_control, yrs 2020-2025", 
       subtitle = "Points are medians w/ 95% CI") +
  scale_color_brewer(palette = "Set2", direction = -1)+
  #guides(color = TRUE)+
  theme()+
  theme_classic()
ggplot2::ggsave(file.path("figures", "SSBratio_ctl_2020_2025_CIs.png"), height = 12, 
                width = 15, units = "in")
