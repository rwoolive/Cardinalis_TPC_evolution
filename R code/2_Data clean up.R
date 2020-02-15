
#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Clean up TPC data for M. cardinalis study populations
#### & make preliminary figures showing RGR & survival across temperature
#### AUTHOR: Rachel Wooliver



library(ggplot2)
library(dplyr)
library(tidyverse)

###
### DATASET
###
data <- read.csv("Raw data/Master TPC 2019.csv", na.strings=c("NA","","--"))
str(data)
###
### Measurements
###
# 1. Leaf number going in
# 2. Leaf number coming out (green)
# 3. Leaf number coming out (brown)



### For leaf number going in AND leaf number coming out green:
# 0 = live seedling with only cotyledons (no true leaves)
# NA = no seedling there
### For leaf number coming out brown:
# 0 = dead seedling with only cotyledons (no true leaves)
# NA = no brown leaves



###
### Make response variables numeric
###
data$Leaf.number.going.in <- as.numeric(data$Leaf.number.going.in); hist(data$Leaf.number.going.in)
data$Leaf.number.coming.out..green. <- as.numeric(as.character(data$Leaf.number.coming.out..green.)); hist(data$Leaf.number.coming.out..green.)
data$Leaf.number.coming.out..brown. <- as.numeric(as.character(data$Leaf.number.coming.out..brown.)); hist(data$Leaf.number.coming.out..brown.)


###
### Remove trays 55-60 (all died during run) and 103-108 (not used in experiment)
###
data <- data[-which(data$Tray %in% c(55:60, 103:108)),] 
# original sample size was 7776, now 6912






###
### Plants that we will not have an RGR value for:
###
# 1. NA leaf number going in
data$NAleafIn <- rep(0, dim(data)[1])
data$NAleafIn[which(is.na(data$Leaf.number.going.in)==TRUE)] <- 1
sum(data$NAleafIn) # 545 did not germinate or died prior to chamber runs
# 2. 0 leaf number going in
data$zeroLeafIn <- rep(0, dim(data)[1])
data$zeroLeafIn[which(data$Leaf.number.going.in==0)] <- 1
sum(data$zeroLeafIn) # 321 only had cotyledons prior to chamber runs
# 3. Accidentally sampled before we counted leaf number out (for gene expression or picture)
data$sampledB4leafOut <- rep(0, dim(data)[1])
data$sampledB4leafOut[which(data$Notes=="Samples for GE 4/10 - Didn't get leaf # out" | data$Notes=="sampled early" | data$Notes=="Sampled early" | data$Notes=="Sampled Early")] <- 1
sum(data$sampledB4leafOut) # 13 were accidentally damaged before final leaf number was measured
# So these are the plants that we will not include in our RGR calculation:
data$outRGR <- rep(0, dim(data)[1])
data$outRGR[which(data$NAleafIn==1 | data$zeroLeafIn==1 | data$sampledB4leafOut==1)] <- 1
dim(data)[1]-sum(as.numeric(data$outRGR)) # 6033 plants we have RGR data for


###
## Make sure leaf number coming out >= leaf number going in, then calculate RGR
###
totalLeavesOut <- rowSums(data[,c("Leaf.number.coming.out..green.", "Leaf.number.coming.out..brown.")], na.rm=TRUE)
totalLeavesOut[which(data$outRGR==1)] <- NA
data$zeros <- rep(0, dim(data)[1])
totalLeavesIn <- rowSums(data[,c("Leaf.number.going.in","zeros")], na.rm=TRUE)
inspect <- which(totalLeavesIn > totalLeavesOut)
data[inspect,c("Tray","Cell","Leaf.number.going.in", "Leaf.number.coming.out..green.", "Leaf.number.coming.out..brown.")]
# Assume that leaves turned brown if not recorded coming out
# (verified using pictures of trays coming out of the chambers)
data$Leaf.number.coming.out..brown.[inspect]<-c(2,4,2,4,8,6,2,2,2,6,6,22,8,8)
# Also need to look at cases where there were zero leaves going in, NA coming out
# These were probably overlooked dead seedlings
# (verified using pictures of trays coming out of the chambers)
zeroIn <- rep(0, dim(data)[1])
zeroIn[which(data$Leaf.number.going.in==0)] <- 1
naOutG <- rep(0, dim(data)[1])
naOutG[which(is.na(data$Leaf.number.coming.out..green.)==TRUE)] <- 1
naOutB <- rep(0, dim(data)[1])
naOutB[which(is.na(data$Leaf.number.coming.out..brown.)==TRUE)] <- 1
naOutBG <- rep(0, dim(data)[1])
naOutBG[which((naOutG + naOutB)==2)] <- 1
data[which((naOutBG+zeroIn)==2),]
data$Leaf.number.coming.out..brown.[which((naOutBG+zeroIn)==2)] <- 0


###
## Calculate RGR
###
data$RGR <- rep(NA, dim(data)[1])
keep <- seq(1, dim(data)[1])[-which(data$outRGR==1)]
data$RGR[keep] <- (rowSums(data[keep,c("Leaf.number.coming.out..green.","Leaf.number.coming.out..brown.")], na.rm=TRUE) - data$Leaf.number.going.in[keep])/(data$Leaf.number.going.in[keep]*7)
# Look at histogram of RGR to make sure all values are positive 
hist(data$RGR); range(na.omit(data$RGR))




###
### Interesting things in addition to RGR and Fv/Fm we might analyze:
###

## 1. Survival
# survival = 0 if the plant died
# and survival = NA if Leaf.number.going.in, Leaf.number.coming.out..green., & Leaf.number.coming.out..brown. are all NA
data$surv <- rep(1,dim(data)[1])
died <- rep(0,dim(data)[1])
died[which(naGoingIn==0 & data$Leaf.number.going.in==data$Leaf.number.coming.out..brown.)] <- 1
died[which(data$Leaf.number.coming.out..green.>0)] <- 0
data$surv[which(died==1)] <- 0
data$surv[which(naGoingIn==1 & data$newGerm==0)]  <- NA








### 
###  Reorganize levels and assign colors to regions
### 
data$Group <- factor(data$Group, levels(data$Group)[c(5:8, 1:4, 9:12)])
data$Pop <- factor(data$Pop, levels(data$Pop)[c(3,4,1,2,5,6)])
data$Reg <- factor(data$Reg, levels(data$Reg)[c(2,1,3)])
data$regCol <- rep("#9E0142", dim(data)[1]) # red
data$regCol[which(data$Reg=="C")] <- "#FEF0A5" # yellow
data$regCol[which(data$Reg=="N")] <- "#5E4FA2" # purple






### 
### Our response variables
### 

### 1. RGR (numeric continuous)
hist(data$RGR)

### 2. surv (numeric 0=died / 1=survived)
data %>% group_by(Temp) %>% summarize(mean = mean(surv, na.rm = TRUE))



# Create new variable for temperature
data$daytimeTemp <- floor(as.numeric(as.character(data$Temp)))
data$daytimeTemp[which(data$Temp.ord==0)] <- 10





### PLOT OVERALL RGR DATA WITH TREATMENT MEANS
ggplot(data, aes(x = daytimeTemp, y = RGR, group=Group)) +
  geom_point(aes(x=daytimeTemp,y=RGR,fill=Group),
             position = position_jitter(width=0.7), 
             size=1.5, shape=21, alpha=0.3) + 
  scale_y_continuous(name="RGR", limits=c(-0.05, 1.1)) + 
  scale_x_continuous(name="Temperature", limits=c(7.5,47.5), breaks=seq(10,45,by=5)) + 
  theme_classic() +
  scale_fill_manual(values=c(rep("#5E4FA2", 4), rep("#FEF0A5", 4), rep("#9E0142", 4))) + 
  facet_wrap(~Group, ncol=4) +
  theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  stat_summary(fun.y = mean, colour = "black", geom = "line", lwd=1.5)



### PLOT OVERALL RGR DATA WITH TREATMENT MEANS BY FAMILY
ggplot(data, aes(x = daytimeTemp, y = RGR, group = FamID)) +
  geom_point(aes(x=daytimeTemp,y=RGR,fill=FamID),
                        position = position_jitter(width=0.7), 
                        size=1.5, shape=21) +
  scale_y_continuous(name="RGR by family", limits=c(-0.05, 1.1)) + 
  scale_x_continuous(name="Temperature", limits=c(7.5,47.5), breaks=seq(10,45,by=5)) + 
  theme_classic() +
  scale_fill_manual(values=rep(rainbow(18),12)) + 
  facet_wrap(~Group, ncol=4) +
  theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  stat_summary(fun.y = mean, geom = "line", lwd=0.25, aes(color=FamID)) +
  scale_color_manual(values=rep(rainbow(18),12))



### PLOT OVERALL SURVIVAL DATA WITH TREATMENT MEANS
ggplot(data, aes(x = daytimeTemp, y = surv, group = Group)) +
  geom_point(aes(colour=Group), alpha=0.55) +
  scale_y_continuous(name="Survival", limits=c(-0.05, 1.05)) + 
  scale_x_continuous(name="Temperature", limits=c(7.5,47.5), breaks=seq(10,45,by=5)) + 
  theme_classic() +
  scale_colour_manual(values=c(rep("#5E4FA2", 4), rep("#FEF0A5", 4), rep("#9E0142", 4))) + 
  facet_wrap(~Group, ncol=4) +
  theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  geom_hline(yintercept=c(0,1), linetype="dashed", color = "black") +
  geom_hline(yintercept=seq(.25,0.75,by=0.25), color = "gray", size=0.15) +
  stat_summary(fun.y = mean, colour = "black", geom = "line", lwd=1.5)





### PLOT OVERALL SURVIVAL DATA WITH FAMILY MEANS
ggplot(data, aes(x = daytimeTemp, y = surv, colour = FamID)) +
  geom_point(aes(colour=FamID), size=0.35) +
  scale_y_continuous(name="Survival by family", limits=c(-0.05, 1.05)) + 
  scale_x_continuous(name="Temperature", limits=c(7.5,47.5), breaks=seq(10,45,by=5)) + 
  theme_classic() +
  scale_colour_manual(values=rep(rainbow(18),12)) + 
  facet_wrap(~Group, ncol=4) +
  theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  geom_hline(yintercept=c(0,1), linetype="dashed", color = "black") +
  geom_hline(yintercept=seq(.25,0.75,by=0.25), color = "gray", size=0.15) +
  stat_summary(fun.y = mean, geom = "line", lwd=0.25)





# Export data for analysis:
write.csv(data, "Processed data/TPC data_cleaned.csv")


