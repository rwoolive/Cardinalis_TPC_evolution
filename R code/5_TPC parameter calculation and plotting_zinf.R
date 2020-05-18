
#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Calculate group mean parameters with credible 
#### intervals and plot TPCs
#### DATE LAST MODIFIED: 2020-18-11 by rcw


#load performr
devtools::install_github("silastittes/performr", local = FALSE, 
                         ref="zin", force=FALSE) 

#load other packages 
library(dplyr)
library(devtools)
library(performr)
library(tidyverse)
library(ggridges)
library(cowplot)
library(ggpubr)
library(agricolae)
library(GGally)

# Other settings
theme_set(theme_cowplot())
options(mc.cores = parallel::detectCores())
extract <- rstan::extract


############ 
# Read the model output/data back in. 
# For some of the columns, we need to make sure they 
# are characters, factors, or continuous numbers. We 
# also need to add Population and Year columns into 
# avDat.
############ 

# TPC prediction intervals
creds_groups_av <- read.csv("Analysis output/creds_groups_av_zinf.csv")[,-1]
creds_groups_av$species <- as.character(creds_groups_av$species)

# posterior draws
tidy_perf_groups_av <- read.csv("Analysis output/tidy_perf_groups_av_zinf.csv")[,-1]

# group means for each tpc parameter
mean_df <- read.csv("Analysis output/mean_df_groups_av_zinf.csv")[,-1]
mean_df$species <- as.factor(mean_df$species)

# group credible intervals for each tpc parameter
mean_df_ci <- read.csv("Analysis output/mean_df_ci_groups_av_zinf.csv")[,-1]
mean_df_ci$species <- as.factor(mean_df_ci$species)

# family-averaged data used in the model
avDat <- read.csv("Analysis output/avDat_zinf.csv")[,-1]
avDat$daytimeTemp <- as.numeric(avDat$daytimeTemp)
# Group 1 is N1 2010, group 2 is N1 2017, etc.
avDat$Pop <- rep("N1", dim(avDat)[1])
avDat$Pop[which(avDat$Group.ord %in% c(3,4))] <- "N2"
avDat$Pop[which(avDat$Group.ord %in% c(5,6))] <- "C1"
avDat$Pop[which(avDat$Group.ord %in% c(7,8))] <- "C2"
avDat$Pop[which(avDat$Group.ord %in% c(9,10))] <- "S1"
avDat$Pop[which(avDat$Group.ord %in% c(11,12))] <- "S2"
avDat$Pop <- as.factor(avDat$Pop)
avDat$Pop <- factor(avDat$Pop, levels(avDat$Pop)[c(3,4,1,2,5,6)])
avDat$Year <- rep(2010, dim(avDat)[1])
avDat$Year[which(avDat$Group.ord %in% c(2,4,6,8,10,12))] <- 2017
avDat$Year <- as.factor(avDat$Year)

# mean RGR and temperature from the family-averaged dataset
meansTot_avDat <- read.csv("Processed data/TPC data_cleaned_av.csv") %>% 
  dplyr::select(daytimeTemp, FamID, Group.ord, RGR) %>% 
  drop_na() %>% 
  summarize(RGR = mean(RGR, na.rm=TRUE),
            Temp = mean(daytimeTemp, na.rm=TRUE))








############ 
############ PLOT TPC'S FROM OUTPUT DATA  
############ FOR BOTH GROUP AND FAMILY RUNS
############ 




############ 
# In creds, back-transform temperature and RGR:
creds_groups_av$x <- creds_groups_av$x + meansTot_avDat$Temp
creds_groups_av$mu <- creds_groups_av$mu * meansTot_avDat$RGR
creds_groups_av$upper <- creds_groups_av$upper * meansTot_avDat$RGR
creds_groups_av$lower <- creds_groups_av$lower * meansTot_avDat$RGR


############ 
# In creds, assign population and year to "species" (groups) for plotting, then order "species"
# REMEMBER that group 1 is N1 2010, group 2 is N2 2010, etc
creds_groups_av$Pop <- rep("N1",dim(creds_groups_av)[1])
creds_groups_av$Pop[which(creds_groups_av$species %in% c("3", "4"))] <- "N2"
creds_groups_av$Pop[which(creds_groups_av$species %in% c("5", "6"))] <- "C1"
creds_groups_av$Pop[which(creds_groups_av$species %in% c("7", "8"))] <- "C2"
creds_groups_av$Pop[which(creds_groups_av$species %in% c("9", "10"))] <- "S1"
creds_groups_av$Pop[which(creds_groups_av$species %in% c("11", "12"))] <- "S2"
creds_groups_av$Pop <- as.factor(creds_groups_av$Pop)
creds_groups_av$Pop <- factor(creds_groups_av$Pop, levels(creds_groups_av$Pop)[c(3,4,1,2,5,6)])
creds_groups_av$Year <- rep("2010",dim(creds_groups_av)[1])
creds_groups_av$Year[which(creds_groups_av$species %in% c("2","4","6","8","10","12"))] <- "2017"
creds_groups_av$Year <- as.factor(creds_groups_av$Year)
creds_groups_av$species <- as.factor(creds_groups_av$species)
creds_groups_av$species <- factor(creds_groups_av$species, levels(creds_groups_av$species)[c(1,5:12,2:4)])


############ 
# In mean_df and mean_df_ci, add pop and year columns
# Group 1 is N1 2010, group 2 is N1 2017, etc.
mean_df$Pop <- factor(rep(c("N1", "N2", "C1", "C2", "S1", "S2"), each=2)) 
mean_df$year <- as.factor(rep(c("2010","2017"), 6)) 
write.csv(mean_df, "Analysis output/group_mean_params_av_zinf.csv")

mean_df_ci <- as.data.frame(mean_df_ci[order(as.numeric(mean_df_ci$species)),])
mean_df_ci$Pop <- factor(rep(c("N1", "N2", "C1", "C2", "S1", "S2"), each=2)) 
mean_df_ci$year <- as.factor(rep(c("2010","2017"), 6)) 



############ 
# Write supplementary table showing mean and CI for tpc parameters 
mean_ci_table <- mean_df %>% 
  right_join(mean_df_ci, by=c("species","species"))
mean_ci_table <- as.data.frame(mean_ci_table)
mean_ci_table <- mean_ci_table[,c("Pop.x", "year.x",
                                  "maximaBT", "maximaBT_lci", "maximaBT_uci",
                                  "B50", "B50_lci", "B50_uci",
                                  "B80", "B80_lci", "B80_uci",
                                  "breadthBT", "breadthBT_lci", "breadthBT_uci",
                                  "x_minBT", "x_minBT_lci", "x_minBT_uci",
                                  "x_maxBT", "x_maxBT_lci", "x_maxBT_uci",
                                  "max_RGR", "max_RGR_lci", "max_RGR_uci",
                                  "area", "area_lci", "area_uci")]
mean_ci_table <- mean_ci_table %>% mutate_at(vars(-c(Pop.x,year.x)), funs(round(., 2)))
mean_ci_table$ToptCI <- paste(mean_ci_table$maximaBT, " [", mean_ci_table$maximaBT_lci, ", ", mean_ci_table$maximaBT_uci,"]", sep="")
mean_ci_table$B50CI <- paste(mean_ci_table$B50, " [", mean_ci_table$B50_lci, ", ", mean_ci_table$B50_uci,"]", sep="")
mean_ci_table$B80CI <- paste(mean_ci_table$B80, " [", mean_ci_table$B80_lci, ", ", mean_ci_table$B80_uci,"]", sep="")
mean_ci_table$breadthCI <- paste(mean_ci_table$breadthBT, " [", mean_ci_table$breadthBT_lci, ", ", mean_ci_table$breadthBT_uci,"]", sep="")
mean_ci_table$x_minCI <- paste(mean_ci_table$x_minBT, " [", mean_ci_table$x_minBT_lci, ", ", mean_ci_table$x_minBT_uci,"]", sep="")
mean_ci_table$x_maxCI <- paste(mean_ci_table$x_maxBT, " [", mean_ci_table$x_maxBT_lci, ", ", mean_ci_table$x_maxBT_uci,"]", sep="")
mean_ci_table$max_RGRCI <- paste(mean_ci_table$max_RGR, " [", mean_ci_table$max_RGR_lci, ", ", mean_ci_table$max_RGR_uci,"]", sep="")
mean_ci_table$areaCI <- paste(mean_ci_table$area, " [", mean_ci_table$area_lci, ", ", mean_ci_table$area_uci,"]", sep="")
mean_ci_table <- mean_ci_table[,c("Pop.x", "year.x",
                                  "ToptCI", 
                                  "B50CI", 
                                  "B80CI", 
                                  "breadthCI",
                                  "x_minCI", 
                                  "x_maxCI", 
                                  "max_RGRCI", 
                                  "areaCI")]
mean_ci_table <- mean_ci_table[c(1,7,2,8,3,9,4,10,5,11,6,12),]
write.csv(mean_ci_table, "Analysis output/Table S4_zinf.csv")


########
# look at ranges in topt, thermal limits, breadth, and pmax
range(mean_df$maximaBT)
range(mean_df$x_minBT)
range(mean_df$x_maxBT)
range(mean_df$B50)
range(mean_df$max_RGR)




############ 
# In avDat, back-transform temperature and RGR, 
avDat$daytimeTemp <- avDat$daytimeTemp + meansTot_avDat$Temp
avDat$RGR <- avDat$RGR * meansTot_avDat$RGR



############
# In tidy dataset, add pop and year
# Group 1 is N1 2010, group 2 is N1 2017, etc.
tidy_perf_groups_av$Pop <- rep("N1", dim(tidy_perf_groups_av)[1])
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("3","4"))] <- "N2"
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("5","6"))] <- "C1"
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("7","8"))] <- "C2"
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("9","10"))] <- "S1"
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("11","12"))] <- "S2"
tidy_perf_groups_av$Pop <- as.factor(tidy_perf_groups_av$Pop)
tidy_perf_groups_av$Pop <- factor(tidy_perf_groups_av$Pop, levels(tidy_perf_groups_av$Pop)[c(3,4,1,2,5,6)])
tidy_perf_groups_av$speciesF <- as.factor(tidy_perf_groups_av$species)
tidy_perf_groups_av$Year <- rep(2010, dim(tidy_perf_groups_av)[1])
tidy_perf_groups_av$Year[which(tidy_perf_groups_av$species %in% c(2,4,6,8,10,12))] <- 2017
tidy_perf_groups_av$Year <- as.factor(tidy_perf_groups_av$Year)





############ 
# Plot combinations of marginal posteriors for parameters using pairs plots: 


# add group as a column in tidy dataframe
tidy_perf_groups_av$group <- paste(tidy_perf_groups_av$Pop, tidy_perf_groups_av$Year, sep=" ")
tidy_perf_groups_av$group <- as.factor(tidy_perf_groups_av$group)
tidy_perf_groups_av$group <- factor(tidy_perf_groups_av$group, levels(tidy_perf_groups_av$group)[c(5:8,1:4,9:12)])

# create new column for labelling facets in the plot
tidy_perf_groups_av$optimum <- tidy_perf_groups_av$maxima

# now plot
pairsplot <- tidy_perf_groups_av %>%
  dplyr::filter(draw %in% sample(1:ndraws_groups_av, 200)) %>%
  dplyr::select(shape1, shape2, x_min, optimum, x_max, stretch, nu, group) %>%
  GGally::ggpairs(1:7, mapping = aes(color = group), switch = "both",
                  diag = list(continuous=wrap("densityDiag", alpha=0.5)),
                  lower = list(continuous = wrap("points", alpha = 0.1, size = 1, shape = 1)), 
                  upper = "blank") +
  theme_bw() +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        strip.text.x = element_text(size=10),strip.text.y = element_text(size=10),
        legend.position="none",
        axis.text.x = element_text(angle = 90, hjust = 1))
pdf(file="Figures/Figure S3_zinf.pdf", height=8, width=8)
pairsplot
dev.off()




# Let's make a plot of tpc's for 2010 (blue) and 2017 (red) cohorts, 
# with 95% credible intervals, for each population, including  
# thermal optima (dotted vertical lines), 
# upper and lower thermal limits (notches on x axis), and
# breadth (dashed horizontal lines),

plwd <- 1.25 # parameter lwd
pla <- 1 # parameter line alpha
tpclwd <- 1.5 # tpc lwd
tpca <- 1 # tpc line alpha
bayesFit_groups_av <- ggplot(data = filter(creds_groups_av, level == 95)) +
  geom_point(data = avDat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35,
             aes(daytimeTemp, RGR, group=Pop, color=Year), inherit.aes = F, size=1.25) + # data for family averages at each temperature
  labs(x = expression(paste("Temperature (", degree, "C)")), 
       y = expression(paste("RGR (leaves leaves"^-1, " day"^-1,")"))) + 
  scale_fill_manual(values=c("blue", "red")) + # fill for year (2010=blue, 2017=red)
  geom_line(aes(x, mu, group=species, color=Year), inherit.aes = F, lwd = tpclwd, alpha=tpca) + # tpc for each group
  lemon::facet_rep_wrap(~Pop, nrow=3, ncol=2, repeat.tick.labels = FALSE) +  # panel for each population
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18),
        legend.title=element_blank(), legend.text = element_text(size = 11), legend.background = element_blank(), legend.position = c(20, 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(size=18, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_segment(data=mean_df, aes(x=maximaBT, y=max_RGR, xend=maximaBT, yend=0,group=Pop, color=year), 
               alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=3) + # vertical lines for thermal optima
  scale_color_manual(values=c("blue", "red"), labels = c("Ancestors", "Descendants"))  +
  geom_segment(data=mean_df, aes(x=maximaBT, y=max_RGR*0.5, xend=B50_high-0.5, yend=max_RGR*0.5, group=Pop, color=year), 
               alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=2,  arrow=arrow(length = unit(0.1, "cm")), lineend="round", linejoin="mitre") + # horizontal arrows for B50 (middle to right end)
  geom_segment(data=mean_df, aes(x=maximaBT, y=max_RGR*0.5, xend=B50_low+0.5, yend=max_RGR*0.5,  group=Pop, color=year), 
               alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=2, arrow=arrow(length = unit(0.1, "cm")), lineend="round", linejoin="mitre") + # horizontal arrows for B50 (middle to left end)
  geom_segment(data=mean_df, aes(x=x_minBT, y=0.075, xend=x_minBT, yend=0, group=Pop, color=year), 
               alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df, aes(x=x_maxBT, y=0.075, xend=x_maxBT, yend=0, group=Pop, color=year), 
               alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # upper thermal limit
  scale_x_continuous(limits=c(6,49), expand = c(0, 0)) +
  scale_y_continuous(limits=c(-0.01,1), expand = c(0, 0), breaks = seq(0,1,by=0.2), minor_breaks = waiver()) 
bayesFit_groups_av <- lemon::reposition_legend(bayesFit_groups_av, 'top left', panel = 'panel-1-1') # put the legend in the first panel







# export tpc plot
ggsave("Figures/Figure 3_zinf.png", 
       ggarrange(bayesFit_groups_av, ncol = 1, nrow = 1), 
       height=8, width=5, dpi=600)
