
#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Create TPCs for each of the 12 M. cardinalis groups
#### and each of the 216 M. cardinalis families
#### using code from Silas Tittes et al. 2019

#### This is the main code we are working with for tpc evolution paper 
#### as of 14 February, 2010 In this code, we model a tpc for each group,
#### using data that are averaged for each family at each temperature
#### total number of groups: 12
#### total number of families: 216 (18 per group)
#### total number of temperatures: 8


# performr vignette by Silas Tittes: https://silastittes.github.io/performr/
# performr implements a probabilistic Bayesian hierarchical model (using Stan) 
# to predict tolerance/performance curves for a set of input taxa. The manuscript 
# that describes the method is Tittes et al. 2019 Am Nat (doi: 10.1086/701827)




#load performr
devtools::install_github("silastittes/performr", local = FALSE) 

#load other libraries 
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
############ LOADING IN AND SCALING DATA
############ 


############ 
# Summary data
dat <- read.csv("Processed data/TPC data_cleaned.csv") 
# how many plants do we have an RGR value for?
dat2 <- (dat[which(dat$outRGR==0),]) # 6033

# look at proportion surviving chamber runs
length(dat$RGR[which(dat$surv==0)])
survproptemp <- dat %>%
  dplyr::group_by(daytimeTemp) %>%
  dplyr::summarize(mean = mean(surv, na.rm = TRUE)) 
1-0.411 # proportion mortality at 10°C: 0.589
1-0.747 # proportion mortality at 10°C: 0.253






# leaf number going into chambers: range, mean, and se
# calculated for all individuals we have an RGR value for
range(na.omit(dat$Leaf.number.going.in[which(dat$outRGR==0)])) 
mean(na.omit(dat$Leaf.number.going.in[which(dat$outRGR==0)])) 
sd(na.omit(dat$Leaf.number.going.in[which(dat$outRGR==0)]))/sqrt(length(na.omit(dat$Leaf.number.going.in[which(dat$outRGR==0)])) )

# leaf number coming out of chambers: range, mean, and se
# calculated for all individuals we have an RGR value for
# first add brown and green leaves 
dat$Leaf.number.coming.out..green.2 <- dat$Leaf.number.coming.out..green.
dat$Leaf.number.coming.out..green.2[which(is.na(dat$Leaf.number.coming.out..green.2))] <- 0
dat$Leaf.number.coming.out..brown.2 <- dat$Leaf.number.coming.out..brown.
dat$Leaf.number.coming.out..brown.2[which(is.na(dat$Leaf.number.coming.out..brown.2))] <- 0
dat$leafout <- dat$Leaf.number.coming.out..green.2+dat$Leaf.number.coming.out..brown.2
range(na.omit(dat$leafout[which(dat$outRGR==0)])) 
mean(na.omit(dat$leafout[which(dat$outRGR==0)])) 
sd(na.omit(dat$leafout[which(dat$outRGR==0)]))/sqrt(length(na.omit(dat$leafout[which(dat$outRGR==0)])) )

# RGR: range, mean, and se
# calculated for all individuals we have an RGR value for
range(na.omit(dat$RGR[which(dat$outRGR==0)])) 
mean(na.omit(dat$RGR[which(dat$outRGR==0)])) 
sd(na.omit(dat$RGR[which(dat$outRGR==0)]))/sqrt(length(na.omit(dat$RGR[which(dat$outRGR==0)])) )




############ 
# Read in data and average by family in each temperature
dat <- read.csv("Processed data/TPC data_cleaned.csv") %>% 
  dplyr::group_by(FamID, daytimeTemp) %>%
  dplyr::summarize(RGR = mean(RGR, na.rm = TRUE)) 
dat$Fam.ord <- as.factor(as.integer(dat$FamID))
dat$Group.ord <- rep(c(1:12), each=(8*18)) # group 1 is N1 2010, group 2 is N2 2010, etc
dat <- filter(dat, !is.na(RGR))
write.csv(dat, "Processed data/TPC data_cleaned_av.csv")
hist(dat$RGR) # raw values range from 0 to 1 
table(dat$daytimeTemp) # raw temp values range from 10 to 15, with up to 216 RGR values (families) per temp

# summary info for RGR
range(na.omit(dat$RGR)) 
mean(na.omit(dat$RGR)) 
sd(na.omit(dat$RGR))/sqrt(length(na.omit(dat$RGR)) )


############ 
# Get raw overall mean of RGR and Temp from family-averaged dataset
# We will use this to scale RGR and center Temp data for the bayesian model
meansTot_avDat <- read.csv("Processed data/TPC data_cleaned_av.csv") %>% 
  dplyr::select(daytimeTemp, FamID, Group.ord, RGR) %>% 
  drop_na() %>% 
  summarize(RGR = mean(RGR, na.rm=TRUE),
            Temp = mean(daytimeTemp, na.rm=TRUE))


############ 
# Scale data (center temp around zero, scale RGR by mean)
avDat <- mutate(dat,
                daytimeTemp = daytimeTemp - meansTot_avDat$Temp,
                RGR = RGR / meansTot_avDat$RGR)
hist(avDat$RGR) # scaled values range from 0 to 4 
table(avDat$daytimeTemp) # centered temp values range from -17.52 to 17.49, with up to 216 RGR values (families) per temp



############ 
# Look at sample sizes of groups 
table(avDat$Group.ord) # across all temps
table(avDat[which(avDat$daytimeTemp==min(avDat$daytimeTemp)),c("Group.ord")]) # at 10/-5 °C
table(avDat[which(avDat$daytimeTemp==max(avDat$daytimeTemp)),c("Group.ord")]) # at 45/30 °C


############ 
# Look at mean and variation in temp and rgr across groups
avDat %>% 
  group_by(Group.ord) %>% 
  summarise(mean_temp = mean(daytimeTemp),
            var_tem = var(daytimeTemp),
            mean_gr = mean(RGR),
            var_gr = var(RGR),
            sampsize = length(RGR))



############ 
# Look at rgr values across temps for each family (color) within each group (panels)
ggplot(avDat, aes(x = daytimeTemp, y = RGR, group = FamID)) +
  geom_point(aes(x=daytimeTemp,y=RGR,fill=FamID),
             position = position_jitter(width=0.7), 
             size=1.5, shape=21, alpha=0.5) +
  guides(fill="none") +
  scale_y_continuous(name="RGR", limits=c(-0.1, 4)) + 
  scale_x_continuous(name="Temperature", limits=c(-20,20), breaks=seq(17.5,17.5,by=5)) + 
  theme_classic() +
  scale_fill_manual(values=rep(rainbow(18),12)) + 
  facet_wrap(~Group.ord, ncol=4) +
  theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  stat_summary(fun.y = mean, geom = "line", lwd=0.25, aes(color=Fam.ord)) +
  scale_color_manual(values=rep(rainbow(18),12))





############ 
############ GROUP RUN WITH DATA AVERAGED BY FAMILY/TEMP COMBO
############ 

perf_out <-
  stan_performance(
    df = avDat,
    response = RGR,
    treatment = daytimeTemp,
    group_ids = Group.ord,
    seed = 1234, max_treedepth = 12,
    file_id =  "Analysis output/model/stan_example_groups_av", 
    iter = 10000
  )




############ 
############ THE FOLLOWING CODE IS FOR PROCESSING OUTPUT  
############ FROM GROUP RUN WITH FAMILY-AVERAGED DATA
############ 


############ 
# load model fits 
model_fits_groups_av <- rstan::read_stan_csv(c("Analysis output/model/stan_example_groups_av.samples_1.csv", 
                                               "Analysis output/model/stan_example_groups_av.samples_2.csv",
                                               "Analysis output/model/stan_example_groups_av.samples_3.csv", 
                                               "Analysis output/model/stan_example_groups_av.samples_4.csv"))



############ 
# The output is the same as any Stan model. n_eff is very large and Rhat is close to 1
knitr::kable(rstan::summary(model_fits_groups_av)$summary, digits = 3)

sumtab <- round(rstan::summary(model_fits_groups_av)$summary, digits=2)
write.csv(sumtab, file="Analysis output/Table S3.csv")




############ 
# We can access the posterior draws using rstan::extract(), which produces a list 
# containing draws for each parameter of the model.
draws_groups_av <- rstan::extract(model_fits_groups_av)
ndraws_groups_av <- length(draws_groups_av$lp__)





############ 
# Visualize uncertainty in curves predictions
# First, let’s generate a tidy data frame with all the parameters, plus some 
# valuable derived parameters, like the optimum, area, and breadth for each species. 
# We will use this data frame for other tasks below as well.

head(
  tidy_perf_groups_av <- 
    performr::perform_df(
      model_fits_groups_av, 
      species_order = c(1:12)
    ) 
)




############ 
# Calculate performance maximum for each iteration of the model
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


###### tidy_perf_groups_av is the dataframe that we can use for pairwise 
# comparisons and plotting
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
write.csv(ssq_df, "Analysis output/model/ssq_df_av.csv") 
ssq_df <- read.csv("Analysis output/model/ssq_df_av.csv") 

View(ssq_df)

#compute bayesian p value
ssq_df %>% 
  summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# overall bayesian p value = 0.1862667

#compute bayesian p value for each group
ssq_group <- ssq_df %>% 
  dplyr::group_by(species) %>%
  summarise(b_pval=mean(ssq_obs > ssq_pseudo))
ssq_group <- as.data.frame(ssq_group)



#create a column for our groups in the overall ssq dataset
ssq_df$speciesF <- ssq_df$species
ssq_df$speciesF[which(ssq_df$species==1)] <- "N1 2010"
ssq_df$speciesF[which(ssq_df$species==7)] <- "N1 2017"
ssq_df$speciesF[which(ssq_df$species==2)] <- "N2 2010"
ssq_df$speciesF[which(ssq_df$species==8)] <- "N2 2017"
ssq_df$speciesF[which(ssq_df$species==3)] <- "C1 2010"
ssq_df$speciesF[which(ssq_df$species==9)] <- "C1 2017"
ssq_df$speciesF[which(ssq_df$species==4)] <- "C2 2010"
ssq_df$speciesF[which(ssq_df$species==10)] <- "C2 2017"
ssq_df$speciesF[which(ssq_df$species==5)] <- "S1 2010"
ssq_df$speciesF[which(ssq_df$species==11)] <- "S1 2017"
ssq_df$speciesF[which(ssq_df$species==6)] <- "S2 2010"
ssq_df$speciesF[which(ssq_df$species==12)] <- "S2 2017"
ssq_df$speciesF <- as.factor(ssq_df$speciesF)
ssq_df$speciesF <- factor(ssq_df$speciesF, levels(ssq_df$speciesF)[c(5:8,1:4,9:12)])

#create a column for groups and merge with their p values in the group p value dataset
ssq_group$ps <- round(ssq_group$b_pval,2)
ssq_group$group <- c("N1 2010", "N2 2010",
                     "C1 2010", "C2 2010",
                     "S1 2010", "S2 2010",
                     "N1 2017","N2 2017",
                     "C1 2017","C2 2017",
                     "S1 2017","S2 2017")
ssq_group$group <- paste(ssq_group$group," (",ssq_group$ps,")", sep="")
# then write the file so we have it saved
write.csv(ssq_group, "Analysis output/p-values-group.csv")  

# in the overall ssq dataset, create a column of groups with their associated p value
ssq_df$group <- as.character(ssq_df$species)
ssq_df <- ssq_df %>%
  mutate(group = fct_recode(group,
                            "N1 2010 (0.34)" = "1",
                            "N2 2010 (0.42)" = "2",
                             "C1 2010 (0.01)" = "3",
                            "C2 2010 (0.55)" = "4",
                            "S1 2010 (0.09)" = "5",
                            "S2 2010 (0.08)" = "6",
                            "N1 2017 (0.42)" = "7",
                            "N2 2017 (0.12)" = "8",
                            "C1 2017 (0.02)" = "9",
                            "C2 2017 (0.03)" = "10",
                            "S1 2017 (0.03)" = "11",
                            "S2 2017 (0.14)" = "12"))
# then reorder to N1 2010, N2 2010, C1 2020, etc.
ssq_df$group <- factor(ssq_df$group, levels(ssq_df$group)[c(1,5,6,7,8,9,10,11,12,2,3,4)])

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
  
  
pdf("Figures/Figure S2.pdf", height=5, width=6.5)
ssq_plot
dev.off()







############ 
# back-transform posterior draws
tidy_perf_groups_av$maximaBT <- tidy_perf_groups_av$maxima + meansTot_avDat$Temp
tidy_perf_groups_av$x_minBT <- tidy_perf_groups_av$x_min + meansTot_avDat$Temp
tidy_perf_groups_av$x_maxBT <- tidy_perf_groups_av$x_max + meansTot_avDat$Temp
tidy_perf_groups_av$breadthBT <- tidy_perf_groups_av$x_maxBT-tidy_perf_groups_av$x_minBT
tidy_perf_groups_av$max_RGR <- tidy_perf_groups_av$max_RGR * meansTot_avDat$RGR



############ 
# CALCULATE B50, LOWER LIMITS, UPPER LIMITS
################
# THE FUNCTION #
# The higher prop_max, the narrower the output breadth will be as the interval is 
# moving nearer to the optimum. The function reports the x axis values the breadth 
# is calculated from: 'opt_breadth_low' & 'opt_breadth_high'. 'opt_breadth' is the 
# difference between low and high, i.e., the breadth. 
################

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
# Write the model output/data to csv. 
############ 

write.csv(creds_groups_av, "Analysis output/creds_groups_av.csv")
write.csv(mean_df, "Analysis output/mean_df_groups_av.csv")
write.csv(mean_df_ci, "Analysis output/mean_df_ci_groups_av.csv")
write.csv(avDat, "Analysis output/avDat.csv")
write.csv(tidy_perf_groups_av, "Analysis output/tidy_perf_groups_av.csv")







############ 
# Read the model output/data back in. 
# Start from here if amending plots, but be sure to load in packages at the top of this R script
############ 

creds_groups_av <- read.csv("Analysis output/creds_groups_av.csv")
creds_groups_av <- creds_groups_av[,-1]
creds_groups_av$species <- as.character(creds_groups_av$species)

tidy_perf_groups_av <- read.csv("Analysis output/tidy_perf_groups_av.csv")
tidy_perf_groups_av <- tidy_perf_groups_av[,-1]

mean_df <- read.csv("Analysis output/mean_df_groups_av.csv")
mean_df <- mean_df[,-1]
mean_df$species <- as.factor(mean_df$species)

mean_df_ci <- read.csv("Analysis output/mean_df_ci_groups_av.csv")
mean_df_ci <- mean_df_ci[,-1]
mean_df_ci$species <- as.factor(mean_df_ci$species)

avDat <- read.csv("Analysis output/avDat.csv")
avDat <- avDat[,-1]
avDat$daytimeTemp <- as.numeric(avDat$daytimeTemp)
avDat$Pop <- rep("N1", dim(avDat)[1])
avDat$Pop[which(avDat$Group.ord %in% c(2,8))] <- "N2"
avDat$Pop[which(avDat$Group.ord %in% c(3,9))] <- "C1"
avDat$Pop[which(avDat$Group.ord %in% c(4,10))] <- "C2"
avDat$Pop[which(avDat$Group.ord %in% c(5,11))] <- "S1"
avDat$Pop[which(avDat$Group.ord %in% c(6,12))] <- "S2"
avDat$Pop <- as.factor(avDat$Pop)
avDat$Pop <- factor(avDat$Pop, levels(avDat$Pop)[c(3,4,1,2,5,6)])
avDat$Year <- rep(2010, dim(avDat)[1])
avDat$Year[which(avDat$Group.ord %in% c(7:12))] <- 2017
avDat$Year <- as.factor(avDat$Year)

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
# In creds, assign population and year to "species" for plotting, then order "species"
creds_groups_av$Pop <- rep("N1",dim(creds_groups_av)[1])
creds_groups_av$Pop[which(creds_groups_av$species %in% c("2", "8"))] <- "N2"
creds_groups_av$Pop[which(creds_groups_av$species %in% c("3", "9"))] <- "C1"
creds_groups_av$Pop[which(creds_groups_av$species %in% c("4", "10"))] <- "C2"
creds_groups_av$Pop[which(creds_groups_av$species %in% c("5", "11"))] <- "S1"
creds_groups_av$Pop[which(creds_groups_av$species %in% c("6", "12"))] <- "S2"
creds_groups_av$Pop <- as.factor(creds_groups_av$Pop)
creds_groups_av$Pop <- factor(creds_groups_av$Pop, levels(creds_groups_av$Pop)[c(3,4,1,2,5,6)])
creds_groups_av$Year <- rep("2010",dim(creds_groups_av)[1])
creds_groups_av$Year[which(creds_groups_av$species %in% c("7","8","9","10","11","12"))] <- "2017"
creds_groups_av$Year <- as.factor(creds_groups_av$Year)
creds_groups_av$species <- as.factor(creds_groups_av$species)
creds_groups_av$species <- factor(creds_groups_av$species, levels(creds_groups_av$species)[c(1,5:12,2:4)])


############ 
# In mean_df and mean_df_ci, add pop and year columns
mean_df$Pop <- factor(rep(c("N1", "N2", "C1", "C2", "S1", "S2"), 2)) 
mean_df$year <- as.factor(rep(c("2010","2017"), each=6)) 
write.csv(mean_df, "Analysis output/group_mean_params_av.csv")

mean_df_ci <- as.data.frame(mean_df_ci[order(as.numeric(mean_df_ci$species)),])
mean_df_ci$Pop <- factor(rep(c("N1", "N2", "C1", "C2", "S1", "S2"), 2)) 
mean_df_ci$year <- as.factor(rep(c("2010","2017"), each=6)) 



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
write.csv(mean_ci_table, "Analysis output/Table S4.csv")


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
tidy_perf_groups_av$Pop <- rep("N1", dim(tidy_perf_groups_av)[1])
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("2","8"))] <- "N2"
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("3","9"))] <- "C1"
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("4","10"))] <- "C2"
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("5","11"))] <- "S1"
tidy_perf_groups_av$Pop[which(tidy_perf_groups_av$species %in% c("6","12"))] <- "S2"
tidy_perf_groups_av$Pop <- as.factor(tidy_perf_groups_av$Pop)
tidy_perf_groups_av$Pop <- factor(tidy_perf_groups_av$Pop, levels(tidy_perf_groups_av$Pop)[c(3,4,1,2,5,6)])
tidy_perf_groups_av$speciesF <- as.factor(tidy_perf_groups_av$species)
tidy_perf_groups_av$Year <- rep(2010, dim(tidy_perf_groups_av)[1])
tidy_perf_groups_av$Year[which(tidy_perf_groups_av$species %in% c(7:12))] <- 2017
tidy_perf_groups_av$Year <- as.factor(tidy_perf_groups_av$Year)





############ 
# Plot combinations of marginal posteriors for parameters using pairs plots: 


# add group as a column in tidy dataframe
tidy_perf_groups_av$group <- paste(tidy_perf_groups_av$Pop, tidy_perf_groups_av$Year, sep=" ")
tidy_perf_groups_av$group <- as.factor(tidy_perf_groups_av$group)
tidy_perf_groups_av$group <- factor(tidy_perf_groups_av$group, levels(tidy_perf_groups_av$group)[c(5:8,1:4,9:12)])

# create new columns for labelling facets in the plot
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
pdf(file="Figures/Figure S3.pdf", height=8, width=8)
pairsplot
dev.off()




# Let's make a plot shows the tpc's for 2010 (blue) and 2017 (red) cohorts, 
# with 95% confidence intervals, for each population, including  
# thermal optima (dotted vertical lines), 
# upper and lower thermal limits (notches on x axis),
# breadth (dashed horizontal lines),
# performance maxima (dotted horizontal lines)

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
  # geom_segment(data=mean_df, aes(x=maximaBT, y=max_RGR, xend=6, yend=max_RGR, group=Pop, color=year), 
  #              alpha=1, lwd=0.75, inherit.aes=FALSE, lty=3) +
  geom_segment(data=mean_df, aes(x=x_minBT, y=0.075, xend=x_minBT, yend=0, group=Pop, color=year), 
               alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df, aes(x=x_maxBT, y=0.075, xend=x_maxBT, yend=0, group=Pop, color=year), 
               alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # upper thermal limit
  scale_x_continuous(limits=c(6,49), expand = c(0, 0)) +
  scale_y_continuous(limits=c(-0.01,1), expand = c(0, 0), breaks = seq(0,1,by=0.2), minor_breaks = waiver()) 
bayesFit_groups_av <- lemon::reposition_legend(bayesFit_groups_av, 'top left', panel = 'panel-1-1') # allows me to put the legend in the first panel







# export plots
pdf("Figures/Figure 3.pdf", height=8, width=5)
figure <- ggarrange(bayesFit_groups_av, ncol = 1, nrow = 1)
figure
dev.off()






###### 
###### PAIRWISE COMPARISONS
###### 

###### params (character vector of the performance parameters) 
params <- c("stretch", "maxima", "maximaBT", "breadth", "breadthBT", 
            "x_min", "x_minBT", "x_max", "x_maxBT", "area", 
            "B50_low", "B50_high", "B50", "B80_low", "B80_high","B80",
            "max_RGR")
###### level (the probability of difference threshold. the table will print 
###### a "*" if the probability of one group having a smaller parameter value 
###### than another is less than 0.05 or greater than 0.95, and two "**" if it 
###### above the level value specified  -- default is 0.99. The values in the 
###### table are the posterior average of the column group minus the row group.)
hi <- 0.99
lo <- 1 - hi

###### pairwise comparison function

pairwise_probs_groups_av <- function(stan_df = stan_df, group, params, level = 0.99){
  group <- enquo(group)
  params %>% walk(function(param){
    
    var_df <- stan_df %>%
      dplyr::select(!!group, param, draw) %>%
      spread(!!group, param) %>%
      dplyr::select(-draw)
    
    spp_names <- names(var_df)
    
    comp_mat <- 
      seq_along(var_df) %>%
      map(~{
        seq_along(var_df) %>%
          map_chr(function(x){
            c1 <- var_df[[x]] 
            c2 <- var_df[[.x]] 
            comp_p <- mean(c1 > c2)
            mean_diff <- round(mean(c1 - c2), 3)
            
            case_when(
              (comp_p > hi | comp_p < lo) & mean_diff != 0 ~ str_glue("{mean_diff}**"),
              (comp_p > 0.975 | comp_p < 0.025) & mean_diff != 0 ~ str_glue("{mean_diff}*"), ### RCW altered this from ST's code that originally used the 0.95/0.05 thresholds
              TRUE ~ str_glue("{mean_diff}")
            )
            
          }) %>%
          rbind()
      }) %>%
      do.call(rbind, .)
    
    dimnames(comp_mat) <- list(spp_names, spp_names)
    
    write.csv(
      x = comp_mat, 
      quote = F, 
      row.names = T,
      file = str_glue("Analysis output/pairwise_tables_groups_av/{param}_pairwise_newthresh.csv") ### RCW altered this from ST's code that originally used the 0.95/0.05 thresholds (original file did not have the "_newthresh")
    )
  })
}

pairwise_probs_groups_av(stan_df = tidy_perf_groups_av, 
                      group = species, 
                      params = params)








###### 
###### PAIRWISE COMPARISON TABLE FOR TOPT AND BREADTH
###### 


# calculate how many iterations estimated greater topt for N1 2017 vs. N1 2010
N1comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==7)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==1)]
length(which(N1comp==TRUE))/length(N1comp)
## then, the average difference with 95% CI
N1compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==7)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==1)]
mN1 <- mean(N1compNum)
qN1 <- quantile(N1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for N2 2017 vs. N2 2010
N2comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==8)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==2)]
length(which(N2comp==TRUE))/length(N2comp)
## then, the average difference with 95% CI
N2compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==8)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==2)]
mN2 <- mean(N2compNum)
qN2 <- quantile(N2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for C1 2017 vs. C1 2010
C1comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==9)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==3)]
length(which(C1comp==TRUE))/length(C1comp)
## then, the average difference with 95% CI
C1compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==9)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==3)]
mC1 <- mean(C1compNum)
qC1 <- quantile(C1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for C2 2017 vs. C2 2010
C2comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==10)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==4)]
length(which(C2comp==TRUE))/length(C2comp)
## then, the average difference with 95% CI
C2compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==10)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==4)]
mC2 <- mean(C2compNum)
qC2 <- quantile(C2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for S1 2017 vs. S1 2010
S1comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==11)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==5)]
length(which(S1comp==TRUE))/length(S1comp)
## then, the average difference with 95% CI
S1compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==11)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==5)]
mS1 <- mean(S1compNum)
qS1 <- quantile(S1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for S2 2017 vs. S2 2010
S2comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==12)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==6)]
length(which(S2comp==TRUE))/length(S2comp)
## then, the average difference with 95% CI
S2compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==12)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==6)]
mS2 <- mean(S2compNum)
qS2 <- quantile(S2compNum, probs = c(0.025,0.975))

# combine into a dataframe
means <- round(c(mN1, mN2, mC1, mC2, mS1, mS2),3)
lowers <- as.vector(round(c(qN1[1], qN2[1], qC1[1], qC2[1], qS1[1], qS2[1]),3))
uppers <- as.vector(round(c(qN1[2], qN2[2], qC1[2], qC2[2], qS1[2], qS2[2]),3))
pwctopt <- paste(means, " [", lowers, ", ", uppers, "]", sep="")



# calculate how many iterations estimated greater b50 for N1 2017 vs. N1 2010
N1comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==7)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==1)]
length(which(N1comp==TRUE))/length(N1comp)
## then, the average difference with 95% CI
N1compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==7)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==1)]
mN1 <- mean(N1compNum)
qN1 <- quantile(N1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for N2 2017 vs. N2 2010
N2comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==8)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==2)]
length(which(N2comp==TRUE))/length(N2comp)
## then, the average difference with 95% CI
N2compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==8)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==2)]
mN2 <- mean(N2compNum)
qN2 <- quantile(N2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for C1 2017 vs. C1 2010
C1comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==9)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==3)]
length(which(C1comp==TRUE))/length(C1comp)
## then, the average difference with 95% CI
C1compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==9)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==3)]
mC1 <- mean(C1compNum)
qC1 <- quantile(C1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for C2 2017 vs. C2 2010
C2comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==10)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==4)]
length(which(C2comp==TRUE))/length(C2comp)
## then, the average difference with 95% CI
C2compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==10)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==4)]
mC2 <- mean(C2compNum)
qC2 <- quantile(C2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for S1 2017 vs. S1 2010
S1comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==11)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==5)]
length(which(S1comp==TRUE))/length(S1comp)
## then, the average difference with 95% CI
S1compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==11)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==5)]
mS1 <- mean(S1compNum)
qS1 <- quantile(S1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for S2 2017 vs. S2 2010
S2comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==12)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==6)]
length(which(S2comp==TRUE))/length(S2comp)
## then, the average difference with 95% CI
S2compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==12)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==6)]
mS2 <- mean(S2compNum)
qS2 <- quantile(S2compNum, probs = c(0.025,0.975))

# combine into a dataframe
means <- round(c(mN1, mN2, mC1, mC2, mS1, mS2),3)
lowers <- as.vector(round(c(qN1[1], qN2[1], qC1[1], qC2[1], qS1[1], qS2[1]),3))
uppers <- as.vector(round(c(qN1[2], qN2[2], qC1[2], qC2[2], qS1[2], qS2[2]),3))
pwcb50 <- paste(means, " [", lowers, ", ", uppers, "]", sep="")


# now combine pwc's into a dataframe and export to csv
pwc <- data.frame(pwctopt=pwctopt,
                  pwcb50=pwcb50)
write.csv(pwc, "Analysis output/Table 1 tpc.csv")



