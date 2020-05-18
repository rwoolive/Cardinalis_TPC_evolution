
#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Model TPCs across 8 temperature regimes for each  
####          of the 12 M. cardinalis groups (18 families/group) 
####          using the heirarchical bayesian performance curve  
####          model of Tittes et al. 2019 (doi: 10.1086/701827)
####          that accounts for zero-inflation
#### DATE LAST MODIFIED: 2020-05-18 by rcw

#### WARNING: THE MODEL (LINE 189) IS COMPUTATIONALLY INTENSIVE AND WILL TAKE TIME TO RUN
# Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3. http://mc-stan.org
# Package performr version 0.2


#load performr
devtools::install_github("silastittes/performr", local = FALSE, 
                         ref="zin", force=FALSE) 

#load other packages 
library(dplyr)
library(devtools)
library(performr)
library(tidyverse)
library(ggridges)
library(ggpubr)


# Other settings
theme_set(theme_cowplot())
options(mc.cores = parallel::detectCores())
extract <- rstan::extract




############ 
############ LOADING IN, AVERAGING, AND SCALING DATA
############ 


############ 
# Summary data
dat <- read.csv("Processed data/TPC data_cleaned.csv") 
# how many plants do we have an RGR value for?
dat <- (dat[which(dat$outRGR==0),]) # 6033

# look at proportion surviving chamber runs
length(dat2$RGR[which(dat$surv==0)]) # 699 died
survproptemp <- dat %>%
  dplyr::group_by(daytimeTemp) %>%
  dplyr::summarize(mean = mean(surv, na.rm = TRUE)) 
1-survproptemp$mean[which(survproptemp$daytimeTemp==10)] # proportion mortality at 10°C: 0.593
1-survproptemp$mean[which(survproptemp$daytimeTemp==45)] # proportion mortality at 45°C: 0.257






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
# Export a reference dataset for the families within each group
# Group 1 is N1 2010, group 2 is N1 2017, etc.
famref <- data.frame(FamID=sort(unique(dat$FamID)),
                     Group=rep(NA,length(unique(dat$FamID))),
                     Group.ord=rep(NA,length(unique(dat$FamID))))
for(i in 1:length(unique(dat$FamID))){
  famref$Group[i] <- as.character(unique(dat$Group[which(dat$FamID==famref$FamID[i])]))
  famref$Group.ord[i] <- as.character(unique(dat$Group.ord[which(dat$FamID==famref$FamID[i])]))
}
write.csv(famref, "Processed data/famref.csv")

############ 
# Read in data and average by family in each temperature, then add group/year info
dat <- read.csv("Processed data/TPC data_cleaned.csv") %>% 
  dplyr::group_by(FamID, daytimeTemp) %>%
  dplyr::summarize(RGR = mean(RGR, na.rm = TRUE)) 
dat$Group <- rep(NA, dim(dat)[1])
dat$Group.ord <- rep(NA, dim(dat)[1])
famref <- read.csv("Processed data/famref.csv")
# Group 1 is N1 2010, group 2 is N1 2017, etc.
for(i in 1:dim(dat)[1]){
  dat$Group[i] <- as.character(famref$Group[which(famref$FamID==dat$FamID[i])])
  dat$Group.ord[i] <- famref$Group.ord[which(famref$FamID==dat$FamID[i])]
}
dat <- filter(dat, !is.na(RGR))
write.csv(dat, "Processed data/TPC data_cleaned_av.csv")
hist(dat$RGR) # raw values range from 0 to 1 
table(dat$daytimeTemp) # raw temp values range from 10 to 15, with up to 216 RGR values (families) per temp
table(dat$Group) # raw temp values range from 10 to 15, with up to 216 RGR values (families) per temp

# summary info for RGR
range(na.omit(dat$RGR)) 
mean(na.omit(dat$RGR)) 
sd(na.omit(dat$RGR))/sqrt(length(na.omit(dat$RGR)) )


############
# Get raw overall mean of RGR and Temp from family-averaged dataset
# We will use this to scale RGR and center Temp data for the bayesian model
meansTot_avDat <- read.csv("Processed data/TPC data_cleaned_av.csv") %>% 
  dplyr::select(daytimeTemp, RGR) %>% 
  drop_na() %>% 
  summarize(RGR = mean(RGR, na.rm=TRUE),
            Temp = mean(daytimeTemp, na.rm=TRUE))
write.csv(meansTot_avDat, "Processed data/meansTot_avDat.csv")


############ 
# Scale data (center temp around zero, scale RGR by mean)
avDat <- mutate(dat,
                daytimeTemp = daytimeTemp - meansTot_avDat$Temp,
                RGR = RGR / meansTot_avDat$RGR)
hist(avDat$RGR) # scaled values range from 0 to 4 
table(avDat$daytimeTemp) # centered temp values range from -17.52 to 17.49, with up to 216 RGR values (families) per temp
write.csv(avDat, "Processed data/avDat.csv")


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
  stat_summary(fun = mean, geom = "line", lwd=0.25, aes(color=FamID)) +
  scale_color_manual(values=rep(rainbow(18),12))





############ 
############ GROUP RUN WITH DATA AVERAGED BY FAMILY/TEMP COMBO
############ WARNING: this is a data hungry analysis! 
############ 

perf_out <-
  stan_performance(
    df = avDat, # dataset: rgr averaged by family x temp combo
    response = RGR, # RGR that has been scaled by the grand mean
    treatment = daytimeTemp, # temperature that has been centered around zero
    group_ids = Group.ord, # groups: Group 1 is N1 2010, group 2 is N1 2017, etc.
    seed = 1234, max_treedepth = 12,
    file_id =  "Analysis output/model/stan_example_groups_av_zinf", 
    iter = 10000
  )
