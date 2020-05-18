
#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Process TPC model output,
####          export model summary (Table S3) 
#### DATE LAST MODIFIED: 2020-05-18 by rcw
#### WARNING: MODEL OUTPUT FILES ARE LARGE AND WILL TAKE TIME TO READ IN.
#### SOME CALCULATIONS (E.G. PERFORMANCE MAXIMUM, P-VALUES, AND B50) 
#### ARE COMPUTATIONALLY INTENSIVE AND WILL TAKE TIME TO RUN.



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
# load model fits 
model_fits_groups_av <- rstan::read_stan_csv(c("Analysis output/model/stan_example_groups_av_zinf.samples_1.csv", 
                                               "Analysis output/model/stan_example_groups_av_zinf.samples_2.csv",
                                               "Analysis output/model/stan_example_groups_av_zinf.samples_3.csv", 
                                               "Analysis output/model/stan_example_groups_av_zinf.samples_4.csv"))


############ 
# read in overall average rgr and temperature
meansTot_avDat <- read.csv("Processed data/meansTot_avDat.csv")


############ 
# read in family-averaged data
avDat <- read.csv("Processed data/avDat.csv")



############ 
# The output is the same as any Stan model. 
# n_eff should be very large and Rhat should be close to 1
sumtab <- round(rstan::summary(model_fits_groups_av)$summary, digits=2)
write.csv(sumtab, file="Analysis output/Table S3_zinf.csv")




############ 
# We can access the posterior draws using rstan::extract(), which produces a list 
# containing draws for each parameter of the model.
draws_groups_av <- rstan::extract(model_fits_groups_av)
ndraws_groups_av <- length(draws_groups_av$lp__)





############ 
# Generate a tidy data frame with all the parameters, plus some 
# valuable derived parameters, like the optimum, area, and breadth for each group 
# We will use this data frame for other tasks below as well.

head(
  tidy_perf_groups_av <- 
    performr::perform_df(
      model_fits_groups_av, 
      species_order = c(1:12)
    ) 
)

# tidy_perf_groups_av is the dataframe that we can use for pairwise 
# comparisons and plotting.
# shape1: First of the two parameters that modifies curve asymmetry; when shape1
# is larger than shape2 , the curve will skew right
# shape2: Second parameter that modifies curve asymmetry; when shape2
# is larger than shape1 , the curve will skew left
# stretch: Dictates the maximum expected value of the response trait (maximum RGR)
# min_max: Performance breadth
# nu: Variance?
# x_min: Location along the environmental axis left of the optimum where the 
# response trait falls to 0 (low temperature threshold)
# x_max: Location along the environmental axis right of the optimum where the 
# response trait falls to 0 (high temperature threshold)




############ 
# Calculate performance maximum for each iteration of the model.
# To get the actual RGR max value, it should just be the value of the Kumaraswamy-ish 
# function at the optimum.  The values and scale will differ, but it should be very 
# correlated with the value of stretch. If you want to calculate it, the easiest 
# would be to use performr::performance_mu(). 
# This mutate() function seems to work to add in the variable for max_RGR in-place.
# Sadly, kinda slow, but do() adds a performance bar, so that's fun! 

tidy_perf_groups_av %<>%
  ungroup() %>% 
  mutate(idx = 1:n()) %>% 
  group_by(idx) %>% 
  do({
    mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))
  })







############ 
# Bayesian p-value: 
# A Bayesian p-value of 0.5 suggests the model generates data that look like the true data.

# the function:
bayes_p <- function(stan_df, raw_df, raw_group, raw_treatment, raw_response, ndraw = nrow(stan_df)){
  
  1:ndraw %>% 
    map_df(~{
      draw_i <- stan_df[.x, ]
      tidy_spp <- stan_df$species[.x]
      df_i <- raw_df[raw_df[[raw_group]] == tidy_spp,]
      
      
      mus <- performance_mu(xs = df_i[[raw_treatment]], 
                            shape1 = draw_i$shape1, 
                            shape2 = draw_i$shape2, 
                            stretch = draw_i$stretch, 
                            x_min = draw_i$x_min,
                            x_max = draw_i$x_max
      )
      
      pseudo <- posterior_predict(x = df_i[[raw_treatment]], draw_i)  
      
      data.frame(
        ssq_obs = sum((mus - df_i[[raw_response]])^2),
        ssq_pseudo = sum((mus - pseudo$trait)^2),
        species = tidy_spp
      )
      
    })
}

#example of function use
# stan_df -- the tidy stan data frame
# raw_df -- the data frame used as input to performr
# raw_group -- the column name that contains the groups the model was fit with *in quotes*
# raw_treatment -- the column name that contains the treatment variable the model was fit with *in quotes* 
# raw_response -- the column name that contains the response variable the model was fit with *in quotes* 

avDat$Group.ord <- as.character(avDat$Group.ord)
ssq_df <- bayes_p(
  stan_df = tidy_perf_groups_av, 
  raw_df = avDat, 
  raw_treatment = "daytimeTemp", 
  raw_response = "RGR", 
  raw_group = "Group.ord")
write.csv(ssq_df, "Analysis output/model/ssq_df_av_zinf.csv") 
ssq_df <- read.csv("Analysis output/model/ssq_df_av_zinf.csv") 

View(ssq_df)

#compute bayesian p value
ssq_df %>% 
  summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# overall bayesian p value = 0.8019292, that's good!

#compute bayesian p value for each group
ssq_group <- ssq_df %>% 
  dplyr::group_by(species) %>%
  summarise(b_pval=mean(ssq_obs > ssq_pseudo))
ssq_group <- as.data.frame(ssq_group)



#create a column for our groups in the overall ssq dataset
# Group 1 is N1 2010, group 2 is N1 2017, etc.
ssq_df$speciesF <- ssq_df$species
ssq_df$speciesF[which(ssq_df$species==1)] <- "N1 2010"
ssq_df$speciesF[which(ssq_df$species==2)] <- "N1 2017"
ssq_df$speciesF[which(ssq_df$species==3)] <- "N2 2010"
ssq_df$speciesF[which(ssq_df$species==4)] <- "N2 2017"
ssq_df$speciesF[which(ssq_df$species==5)] <- "C1 2010"
ssq_df$speciesF[which(ssq_df$species==6)] <- "C1 2017"
ssq_df$speciesF[which(ssq_df$species==7)] <- "C2 2010"
ssq_df$speciesF[which(ssq_df$species==8)] <- "C2 2017"
ssq_df$speciesF[which(ssq_df$species==9)] <- "S1 2010"
ssq_df$speciesF[which(ssq_df$species==10)] <- "S1 2017"
ssq_df$speciesF[which(ssq_df$species==11)] <- "S2 2010"
ssq_df$speciesF[which(ssq_df$species==12)] <- "S2 2017"
ssq_df$speciesF <- as.factor(ssq_df$speciesF)
ssq_df$speciesF <- factor(ssq_df$speciesF, levels(ssq_df$speciesF)[c(5:8,1:4,9:12)])

#create a column for groups and merge with their p values in the group p value dataset
ssq_group$ps <- round(ssq_group$b_pval,2)
ssq_group$group <- c("N1 2010",  "N1 2017","N2 2010","N2 2017",
                     "C1 2010", "C1 2017","C2 2010","C2 2017",
                     "S1 2010", "S1 2017","S2 2010","S2 2017")
ssq_group$group <- paste(ssq_group$group," (",ssq_group$ps,")", sep="")
# then write the file so we have it saved
write.csv(ssq_group, "Analysis output/p-values-group_zinf.csv")  

# in the overall ssq dataset, create a column of groups with their associated p value
# REMEMBER that group 1 is N1 2010, group 2 is N2 2010, etc
ssq_df$group <- as.character(ssq_df$species)
ssq_df <- ssq_df %>%
  mutate(group = fct_recode(group,
                            "N1 2010 (0.91)" = "1",
                            "N1 2017 (0.91)" = "2",
                            "N2 2010 (0.94)" = "3",
                            "N2 2017 (0.85)" = "4",
                            "C1 2010 (0.61)" = "5",
                            "C1 2017 (0.70)" = "6",
                            "C2 2010 (0.95)" = "7",
                            "C2 2017 (0.71)" = "8",
                            "S1 2010 (0.72)" = "9",
                            "S1 2017 (0.71)" = "10",
                            "S2 2010 (0.77)" = "11",
                            "S2 2017 (0.86)" = "12"))
# then reorder to N1 2010, N1 2017, etc.
ssq_df$group <- factor(ssq_df$group, levels(ssq_df$group)[c(1,5:12,2:4)])

# subsample to make a figure
subsamp <-  ssq_df%>%
  group_by(group)%>%
  sample_frac(size = 0.1, replace = F)

# figure of predicted vs. observed ssq
ssq_plot <- ggplot(subsamp, aes(y=log(ssq_pseudo), x=log(ssq_obs), color=group)) +
  geom_point(alpha=0.1) +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  xlab("Log observed ssq") +
  ylab("Log posterior predictive ssq") +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18),
        #legend.title=element_blank(), legend.text = element_text(size = 11), legend.background = element_blank(), legend.position = c(20, 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(size=18, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ssq_plot  


pdf("Figures/Figure S2_zinf.pdf", height=5, width=6.5)
ssq_plot
dev.off()
############ 







############ 
# back-transform posterior draws
tidy_perf_groups_av$maximaBT <- tidy_perf_groups_av$maxima + meansTot_avDat$Temp
tidy_perf_groups_av$x_minBT <- tidy_perf_groups_av$x_min + meansTot_avDat$Temp
tidy_perf_groups_av$x_maxBT <- tidy_perf_groups_av$x_max + meansTot_avDat$Temp
tidy_perf_groups_av$breadthBT <- tidy_perf_groups_av$x_maxBT-tidy_perf_groups_av$x_minBT
tidy_perf_groups_av$max_RGR <- tidy_perf_groups_av$max_RGR * meansTot_avDat$RGR



############ 
# CALCULATE B50, LOWER LIMITS, UPPER LIMITS
############ 
# THE FUNCTION #
# The higher prop_max, the narrower the output breadth will be as the interval is 
# moving nearer to the optimum. The function reports the x axis values the breadth 
# is calculated from: 'opt_breadth_low' & 'opt_breadth_high'. 'opt_breadth' is the 
# difference between low and high, i.e., the breadth. 

optimum_breadth <- 
  function(par_df, prop_max = 0.5, n_grid = 100){
    1:nrow(par_df) %>% map_df(function(i){
      xs = seq(par_df$x_min[i],
               par_df$x_max[i],
               length.out = n_grid)
      
      sim_mu <- performr::performance_mu(
        xs = xs, 
        shape1 = par_df$shape1[i],
        shape2 = par_df$shape2[i],
        stretch = par_df$stretch[i],
        x_min = par_df$x_min[i],
        x_max = par_df$x_max[i]
      )
      
      max_y <- which.max(sim_mu)
      prop_max_y <- sim_mu[max_y] * prop_max
      idx_max_low <- which.min((sim_mu[1:(max_y-1)] - prop_max_y)^2)
      idx_max_high <- which.min((sim_mu[max_y:length(sim_mu)] - prop_max_y)^2) + max_y
      tibble("opt_breadth_low" = xs[idx_max_low],
             "opt_breadth_high" = xs[idx_max_high], 
             "opt_breadth" = xs[idx_max_high] - xs[idx_max_low])
    })
  }

#########
# USAGE #
#########

#gotta do ungroup for now, unfortunately
#increase n_grid for more precise approximation
tidy_perf_groups_av %<>% 
  ungroup() %>% 
  bind_cols(
    optimum_breadth(tidy_perf_groups_av, prop_max = 0.5, n_grid = 100)
  )
# opt_breadth_low: lower limit for B50
# opt_breadth_high: upper limit for B50
# opt_breadth: difference between high and low


# Back-transorm
tidy_perf_groups_av$B50_low <- tidy_perf_groups_av$opt_breadth_low + meansTot_avDat$Temp
tidy_perf_groups_av$B50_high <- tidy_perf_groups_av$opt_breadth_high + meansTot_avDat$Temp
tidy_perf_groups_av$B50 <- tidy_perf_groups_av$B50_high - tidy_perf_groups_av$B50_low



############ 
# Now calculate B80
############ 
tidy_perf_groups_av %<>% 
  ungroup() %>% 
  bind_cols(
    optimum_breadth(tidy_perf_groups_av, prop_max = 0.8, n_grid = 100)
  )
# opt_breadth_low: lower limit for B80
# opt_breadth_high: upper limit for B80
# opt_breadth: difference between high and low


# Back-transorm
tidy_perf_groups_av$B80_low <- tidy_perf_groups_av$opt_breadth_low1 + meansTot_avDat$Temp
tidy_perf_groups_av$B80_high <- tidy_perf_groups_av$opt_breadth_high1 + meansTot_avDat$Temp
tidy_perf_groups_av$B80 <- tidy_perf_groups_av$B80_high - tidy_perf_groups_av$B80_low
############ 












############ 
# get mean values of each parameter
mean_df <- tidy_perf_groups_av %>% 
  group_by(species) %>% 
  summarise_if(is.numeric, .funs = c(mean))
mean_df$species <- as.integer(mean_df$species)
mean_df <- mean_df[order(mean_df$species),]
mean_df$species <- as.factor(mean_df$species)

mean_df_ci <- tidy_perf_groups_av %>% 
  group_by(species) %>% 
  summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)),
            maximaBT_uci = quantile(maximaBT, probs=c(0.975)),
            B50_lci = quantile(B50, probs=c(0.025)),
            B50_uci = quantile(B50, probs=c(0.975)),
            B80_lci = quantile(B80, probs=c(0.025)),
            B80_uci = quantile(B80, probs=c(0.975)),
            breadthBT_lci = quantile(breadthBT, probs=c(0.025)),
            breadthBT_uci = quantile(breadthBT, probs=c(0.975)),
            x_minBT_lci = quantile(x_minBT, probs=c(0.025)),
            x_minBT_uci = quantile(x_minBT, probs=c(0.975)),
            x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)),
            x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)),
            max_RGR_lci = quantile(max_RGR, probs=c(0.025)),
            max_RGR_uci = quantile(max_RGR, probs=c(0.975)),
            stretch_lci = quantile(stretch, probs=c(0.025)),
            stretch_uci = quantile(stretch, probs=c(0.975)),
            area_lci = quantile(area, probs=c(0.025)),
            area_uci = quantile(area, probs=c(0.975)))
mean_df_ci$species <- as.integer(mean_df_ci$species)
mean_df_ci <- mean_df_ci[order(mean_df_ci$species),]
mean_df_ci$species <- as.factor(mean_df_ci$species)







############ 
# The next step is to generate prediction intervals using the predict_interval() 
# function. The function generates posterior quantiles for each set of posterior 
# draws specified (x_draws), and averages over the quantiles. The argument p can 
# takes a vector of credible levels, which can be modified to consider other and/or 
# additional levels.

############ 
# set up sequence from smallest to largest x 

x_seq_groups_av = seq(min(draws_groups_av$x_min),
                      max(draws_groups_av$x_max),
                      length.out = 100)





############ 
# sub-sample draws randomly

poly_draws_groups_av <- sample(1:ndraws_groups_av, 100)



predict_interval <- function (x, spp, par_df, x_draws, p){
  if (missing(x)) {
    x <- seq(min(par_df$x_min), max(par_df$x_max), length.out = 100)
  }
  sub_df <- par_df %>% filter(draw %in% x_draws)
  
  p %>% map_df(~{
    posterior_quantile(x = x, par_df = sub_df, p = .x) %>%
      group_by(species, x) %>%
      summarise_all(.funs = mean) %>%
      mutate(level = .x) %>%
      arrange(x) %>%
      dplyr::select(-draw) %>%
      mutate(level = round(.x * 100, 0))
  }) %>%
    arrange(species, x)
}


head(
  creds_groups_av <- 
    predict_interval(
      x = x_seq_groups_av,
      spp = species,
      par_df = tidy_perf_groups_av,
      x_draws = poly_draws_groups_av,
      p = c(0.95, 0.5))
)



# Notice the output of predict_interval() also produces a “mu” column, which 
# contains the mean preduction of the curve for each input species.



############ 
# Write the model output/data to csv files
############ 

write.csv(creds_groups_av, "Analysis output/creds_groups_av_zinf.csv")
write.csv(mean_df, "Analysis output/mean_df_groups_av_zinf.csv")
write.csv(mean_df_ci, "Analysis output/mean_df_ci_groups_av_zinf.csv")
write.csv(avDat, "Analysis output/avDat_zinf.csv")
write.csv(tidy_perf_groups_av, "Analysis output/tidy_perf_groups_av_zinf.csv")

