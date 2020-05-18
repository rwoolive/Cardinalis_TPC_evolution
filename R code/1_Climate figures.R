
#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Plot climate data for M. cardinalis study populations 
####          Range map with anomalies in seasonality and max July temp (Fig. 2)  
####          Variation in anomalies across study years (Fig. S5)
####          Table showing coordinates & climate info by population (Table S1)
####          Climate data were derived from Climate WNA v. 5.51 
#### AUTHOR: Seema Sheth
#### DATE LAST MODIFIED: 2020-02-18 by rcw



###########################
# load packages
###########################
library(plyr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(ggmap)
library(maps)
library(mapdata)
library(mapproj)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)





###########################
# read in .csv file of climate WNA data for each population and create new 
# grouping variable for year
###########################
dat=read_csv("Raw data/Population_localities_1901-2017MSYT.csv") %>%
  arrange(Year) %>%
  mutate(ID1=factor(ID1),ID2=factor(ID2))



###########################
# read in cardinalis occurrences from Baja 
# (from Angert et al. 2018 Am Nat, doi: 10.1086/695984)
###########################
localities_Baja=read_csv("Raw data/baja.csv") %>%
  filter(grepl(pattern="Baj", Point)) %>% # remove AZ and Cedros Island localities
  dplyr::select(Lat,Lon) %>%
  dplyr::rename(Latitude=Lat,Longitude=Lon)




###########################
# read in cardinalis occurrences from US 
# (from Angert et al. 2018 Am Nat, doi: 10.1086/695984)
###########################
localities=read_csv("Raw data/all.records.aug.31.csv") %>%
  dplyr::filter(DATASET=="herb"&PRESABS==1) %>% # remove absence data
  dplyr::select(Latitude,Longitude) %>%
  dplyr::bind_rows(localities_Baja) # merge localities from US and Baja




###########################
# Calculate seasonality:
# difference between monthly min and max of all months
###########################
dat$maxT <- NA
dat$minT <- NA
dat$seasonality <- NA
for(i in 1:dim(dat)[1]){
  dat$maxT[i] <- max(dat$Tmax01[i], dat$Tmax02[i], dat$Tmax03[i], dat$Tmax04[i], dat$Tmax05[i], dat$Tmax06[i], dat$Tmax07[i], dat$Tmax08[i], dat$Tmax09[i], dat$Tmax10[i], dat$Tmax11[i], dat$Tmax12[i])
  dat$minT[i] <- min(dat$Tmin01[i], dat$Tmin02[i], dat$Tmin03[i], dat$Tmin04[i], dat$Tmin05[i], dat$Tmin06[i], dat$Tmin07[i], dat$Tmin08[i], dat$Tmin09[i], dat$Tmin10[i], dat$Tmin11[i], dat$Tmin12[i])
  }
dat$seasonality <- dat$maxT - dat$minT







###########################
# subset climate data from historical (1951-2000) vs. study years (2010-2017)
###########################
clim_hist= dat %>%
  filter(Year>1950 & Year<2001) %>% 
  arrange(-Latitude)
clim_study=dat %>%
  filter(Year>2009) %>% 
  arrange(-Latitude)


###########################
# make file of lat/lon only
###########################
study_pops=clim_study %>%
  filter(Year==2010)


###########################
# order populations by latitude from north to south
###########################
clim_hist$ID1 = factor(clim_hist$ID1,levels = c("S2","S1","C2","C1","N2","N1"))
clim_study$ID1 = factor(clim_study$ID1,levels = c("S2","S1","C2","C1","N2","N1"))
study_pops$ID1= factor(study_pops$ID1,levels = c("S2","S1","C2","C1","N2","N1"))






###########################
# Calculate historical seasonality for each population and recent anomalies.
# Historical seasonality is greatest in Central pops, lowest in Northern pops
# average anomaly is negative in Central/Southern pops, positive in Northern pops
###########################
avS <- clim_hist %>%
  group_by(ID1) %>%
  dplyr::summarize(meanS = mean(seasonality, na.rm=TRUE))
avS <- avS[c(6:1),]
avS$dev2010S <- clim_study$seasonality[which(clim_study$Year==2010)] - avS$meanS
avS$dev2011S <- clim_study$seasonality[which(clim_study$Year==2011)] - avS$meanS
avS$dev2012S <- clim_study$seasonality[which(clim_study$Year==2012)] - avS$meanS
avS$dev2013S <- clim_study$seasonality[which(clim_study$Year==2013)] - avS$meanS
avS$dev2014S <- clim_study$seasonality[which(clim_study$Year==2014)] - avS$meanS
avS$dev2015S <- clim_study$seasonality[which(clim_study$Year==2015)] - avS$meanS
avS$dev2016S <- clim_study$seasonality[which(clim_study$Year==2016)] - avS$meanS
avS$dev2017S <- clim_study$seasonality[which(clim_study$Year==2017)] - avS$meanS
avS <- avS %>% mutate(avDevS = rowMeans(dplyr::select(., starts_with("dev"))))
avS

write.csv(avS, "Processed data/average seasonality by pop.csv")




###########################
# historical MAT for each population and recent anomalies
# historical MAT increases as you go south
# average anomaly in MAT is greatest in the center, lowest in the north
###########################
avMAT <- clim_hist %>%
  group_by(ID1) %>%
  dplyr::summarize(meanMAT = mean(MAT, na.rm=TRUE))
avMAT <- avMAT[c(6:1),]
avMAT$dev2010MAT <- clim_study$MAT[which(clim_study$Year==2010)] - avMAT$meanMAT
avMAT$dev2011MAT <- clim_study$MAT[which(clim_study$Year==2011)] - avMAT$meanMAT
avMAT$dev2012MAT <- clim_study$MAT[which(clim_study$Year==2012)] - avMAT$meanMAT
avMAT$dev2013MAT <- clim_study$MAT[which(clim_study$Year==2013)] - avMAT$meanMAT
avMAT$dev2014MAT <- clim_study$MAT[which(clim_study$Year==2014)] - avMAT$meanMAT
avMAT$dev2015MAT <- clim_study$MAT[which(clim_study$Year==2015)] - avMAT$meanMAT
avMAT$dev2016MAT <- clim_study$MAT[which(clim_study$Year==2016)] - avMAT$meanMAT
avMAT$dev2017MAT <- clim_study$MAT[which(clim_study$Year==2017)] - avMAT$meanMAT
avMAT <- avMAT %>% mutate(avDevMAT= rowMeans(dplyr::select(., starts_with("dev")))) 
avMAT

write.csv(avMAT, "Processed data/average MAT by pop.csv")




###########################
# historical july max temp for each population and recent anomalies
# historical july max temp increases as you go south
# average anomaly in july max temp is greatest in the center & north, decreases a little in the south
###########################
avJMT <- clim_hist %>%
  group_by(ID1) %>%
  dplyr::summarize(meanJMT = mean(Tmax07, na.rm=TRUE))
avJMT <- avJMT[c(6:1),]
avJMT$dev2010 <- clim_study$Tmax07[which(clim_study$Year==2010)] - avJMT$meanJMT
avJMT$dev2011 <- clim_study$Tmax07[which(clim_study$Year==2011)] - avJMT$meanJMT
avJMT$dev2012 <- clim_study$Tmax07[which(clim_study$Year==2012)] - avJMT$meanJMT
avJMT$dev2013 <- clim_study$Tmax07[which(clim_study$Year==2013)] - avJMT$meanJMT
avJMT$dev2014 <- clim_study$Tmax07[which(clim_study$Year==2014)] - avJMT$meanJMT
avJMT$dev2015 <- clim_study$Tmax07[which(clim_study$Year==2015)] - avJMT$meanJMT
avJMT$dev2016 <- clim_study$Tmax07[which(clim_study$Year==2016)] - avJMT$meanJMT
avJMT$dev2017 <- clim_study$Tmax07[which(clim_study$Year==2017)] - avJMT$meanJMT
avJMT <- avJMT %>% mutate(avDevJMT= rowMeans(dplyr::select(., starts_with("dev")))) 
avJMT

write.csv(avJMT, "Processed data/average JMT by pop.csv")





###########################
# calculate historical precipitation 
# historical MAP is greatest in N, then C, and lowest in S
precip_hist=clim_hist %>%
  group_by(ID1) %>%
  summarize(mean_precip=mean(MAP))

# calculate recent precipitation (by rcw on 20191001)
# recent MAP is greatest in N, then C, and lowest in S
precip_study=clim_study %>%
  group_by(ID1) %>%
  summarize(mean_precip=mean(MAP))

# calculate precipitation anomaly (by rcw on 20191001)
# recent MAP is greatest in N, then C, and lowest in S
avPrecip= data.frame(ID1 = precip_hist$ID1,
                        hist = precip_hist$mean_precip,
                        study = precip_study$mean_precip,
                        anom = precip_study$mean_precip - precip_hist$mean_precip)
avPrecip <- avPrecip[c(6:1),]
write.csv(avPrecip, "Processed data/average MAP by pop and anom.csv")
###########################




###########################
# Table S1: Coordinates and climate data for each population
###########################
popdat <- data.frame(pop=localities$ID1,
                     lat=localities$,
                     long=localities$,
                     elev=Localities$,
                     mat.anom=avMAT$,
                     jmt.anom=avJMT$,
                     precip.anom=avPrecip$,
                     seasonality.anom=avS$)
write.csv(popdat, "Processed data/Table S1")






###########################
### FIGURE 2 BY SNS ON 20190930: updated by RCW 20200127
###########################
clim_study$ID3 <- as.numeric(clim_study$ID2)
clim_study$Year2 <- rep(0,dim(clim_study)[1])
clim_study$Year2[which(clim_study$Year==2011)] <- 1
clim_study$Year2[which(clim_study$Year==2012)] <- 2
clim_study$Year2[which(clim_study$Year==2013)] <- 3
clim_study$Year2[which(clim_study$Year==2014)] <- 4
clim_study$Year2[which(clim_study$Year==2015)] <- 5
clim_study$Year2[which(clim_study$Year==2016)] <- 6
clim_study$Year2[which(clim_study$Year==2017)] <- 7

#********************************************
# Fig. 2a: Map of cardinalis focal populations and locality data from  
# herbarium specimens with clean and simple background
# Cite Angert et al. 2018 AmNat for source of locality data
#********************************************

# get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# create map of states
states <- map_data("state")

# make map object of CA and OR
CA_OR <- states %>%
  filter(region %in% c("california", "oregon"))

# make bounding box
ymin = min(localities$Latitude) - 1.5
ymax = max(localities$Latitude) + 1.5
xmin = min(localities$Longitude) - 1 # latitude labels look nicer bounding by demography sites vs. all locality data
xmax = max(localities$Longitude) + 1 # latitude labels look nicer bounding by demography sites vs. all locality data
e = extent(xmin, xmax,ymin, ymax)

# make map; modified from M_cardinalis_range_map.R; unlike map above, this map includes locality data for M. cardinalis
card_map=ggplot(data=world,fill="lightgrey",col="black",size=0.3) + 
  geom_sf() +
  geom_polygon(aes(x = long, y = lat, group = group), data=states,fill="transparent",col="black",size=0.2) +
  theme_minimal() +
  coord_sf(xlim = c(xmin,xmax), ylim = c(ymin+0.5,ymax), expand = FALSE) +
  geom_point(aes(x=Longitude,y=Latitude),data=localities,alpha=0.75,shape=21,size=3,fill="white") +
  geom_point(data=study_pops,aes(x=Longitude,y=Latitude,fill=ID1),size=8,alpha=0.75,shape=21) +
  geom_text(data=study_pops, aes(x=Longitude,y=Latitude,label=ID1),hjust=-0.5, vjust=c(0.5,0.5,-0.75,1,-0.5,1),size=8,fontface="bold") +
  scale_fill_manual(values=c("#9E0142","#9E0142","#FEF0A5","#FEF0A5","#5E4FA2","#5E4FA2")) + # use this for colored lines
  labs(x="Longitude",y="Latitude",fill="ID1",data=study_pops) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=20),
        legend.position = "none",panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ggtitle('A') + theme(plot.title=element_text(hjust=0,size=24)) 
  

# Note: sometimes I get the following error: "Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : polygon edge not found"
# Which is discussed here: https://github.com/tidyverse/ggplot2/issues/2252
# Typically it works every few times I try, and others have a similar experience

#********************************************
# Fig. 2b: July maximum temperature with anomalies
#********************************************
jmt_plot2=ggplot(clim_hist, aes(group=ID1,x=ID1, y=Tmax07, fill=ID1,colour=ID1)) + 
  geom_boxplot(size=1.3,width=0.55,alpha=0.8,outlier.shape = NA) + # boxplots of historical values
  guides(fill=FALSE,color=FALSE) + # no legends
  coord_flip() + # x axis becomes y axis, y axis becomes x axis
  scale_fill_manual(values=c("#9E0142","#9E0142","#FEF0A5","#FEF0A5","#5E4FA2","#5E4FA2")) + # fill of boxplots
  scale_color_manual(values=rep("black",6)) + # black outline for boxplots
  geom_point(aes(fill=ID1), alpha=0.75, position=position_jitterdodge(jitter.width=1), size=2) + # points for historical values
  geom_point(aes(y=Tmax07,x=ID3+0.4),clim_study,shape=23,size=5,fill="grey",colour="black") + # points for recent values
  ylab(expression(paste("Maximum July temperature (°C)"))) + ylim(24, 36.2) + # y label and limits
  xlab("Population") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),axis.text=element_text(size=18),
                     axis.title=element_text(size=20),panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ggtitle('B') + theme(plot.title=element_text(hjust=0,size=24)) 


#********************************************
# Fig. 2c: Seasonality with anomalies
#********************************************
seasonality_plot=ggplot(clim_hist, aes(group=ID1,x=ID1, y=seasonality, fill=ID1,colour=ID1)) + 
  geom_boxplot(size=1.3,width=0.55,alpha=0.8,outlier.shape = NA) + # boxplots of historical values
  guides(fill=FALSE,color=FALSE) + # no legends
  coord_flip() + # x axis becomes y axis, y axis becomes x axis
  scale_fill_manual(values=c("#9E0142","#9E0142","#FEF0A5","#FEF0A5","#5E4FA2","#5E4FA2")) + # fill of boxplots
  scale_color_manual(values=rep("black",6)) + # black outline for boxplots
  geom_point(aes(fill=ID1), alpha=0.75, position=position_jitterdodge(jitter.width=1), size=2) + # points for historical values
  geom_point(aes(y=seasonality,x=ID3+0.4),clim_study,shape=23,size=5,fill="grey",colour="black") + # points for recent values
  ylab(expression(paste("Seasonality (°C)"))) + ylim(24, 36.8) + # y label and limits
  xlab("Population") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),axis.text=element_text(size=18),
                     axis.title=element_text(size=20),panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ggtitle('C') + theme(plot.title=element_text(hjust=0,size=24)) 




# print to pdf: ALT FIGURE 2 FOR TPC MANUSCRIPT 
# shows map, historical Jul Tmax with anomalies, historical seasonality with anomalies
ggsave("Figures/Figure 2.png", 
       grid.arrange(card_map,jmt_plot2,seasonality_plot,ncol=3, nrow=1, widths=c(1,1,1),heights=1), 
       width=15, height=7.5, dpi=600)








####################
## Historical mean+se temperatures and recent trends
####################

##### Look at variation in July maximum in the 7 recent years (2010-2017)

# ANCOVA
mod <- lm(Tmax07 ~ Year*ID1, data=clim_study)
mod1 <- lm(Tmax07 ~ Year+ID1, data=clim_study)
mod2 <- lm(Tmax07 ~ Year, data=clim_study)
mod3 <- lm(Tmax07 ~ Year+ID1+I(Year^2), data=clim_study)
MuMIn::AICc(mod, mod1, mod2, mod3)
summary(mod1) # b for year is 0.25774, r2adj is 0.8188
anova(mod1)


# Now create predicted lines for each population
N1pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("N1",500))
N1pred$fit <- predict(mod1, newdata = N1pred, re.form = NA, type="response")
N2pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("N2",500))
N2pred$fit <- predict(mod1, newdata = N2pred, re.form = NA, type="response")
C1pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("C1",500))
C1pred$fit <- predict(mod1, newdata = C1pred, re.form = NA, type="response")
C2pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("C2",500))
C2pred$fit <- predict(mod1, newdata = C2pred, re.form = NA, type="response")
S1pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("S1",500))
S1pred$fit <- predict(mod1, newdata = S1pred, re.form = NA, type="response")
S2pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("S2",500))
S2pred$fit <- predict(mod1, newdata = S2pred, re.form = NA, type="response")

avJMT$xs <- rep(2008,6)+c(-0.2,0.1,-0.1,0.2,0,0.3)
sds <- as.data.frame(clim_hist %>%
  group_by(ID1) %>%
  summarise(Tmax07.sd = sd(Tmax07)))
avJMT$sd <- rev(sds$Tmax07.sd)
avJMT$meanJMTuppersd <- avJMT$meanJMT+avJMT$sd
avJMT$meanJMTlowersd <- avJMT$meanJMT-avJMT$sd
avJMT$se <- rev(sds$Tmax07.sd)/sqrt(length(unique(clim_hist$Year)))
avJMT$meanJMTupperse <- avJMT$meanJMT+avJMT$se
avJMT$meanJMTlowerse <- avJMT$meanJMT-avJMT$se

JMTlmSE <- ggplot() +
  geom_point(data = clim_study, aes(y = Tmax07, x = Year, fill = ID1), 
             alpha=0.5, shape=21, colour="black", size=4) + 
  geom_line(data = N1pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5)  +
  geom_line(data = N2pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5)  +
  geom_line(data = C1pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5)  +
  geom_line(data = C2pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5)  +
  geom_line(data = S1pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5 )  +
  geom_line(data = S2pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5 )  +
  geom_line(data = N1pred, aes(x = Year, y = fit), size = 1.5, color="#5E4FA2", alpha=0.5)  +
  geom_line(data = N2pred, aes(x = Year, y = fit), size = 1.5, color="#5E4FA2", alpha=0.5)  +
  geom_line(data = C1pred, aes(x = Year, y = fit), size = 1.5, color="#FEF0A5", alpha=0.5)  +
  geom_line(data = C2pred, aes(x = Year, y = fit), size = 1.5, color="#FEF0A5", alpha=0.5)  +
  geom_line(data = S1pred, aes(x = Year, y = fit), size = 1.5, color="#9E0142", alpha=0.5 )  +
  geom_line(data = S2pred, aes(x = Year, y = fit), size = 1.5, color="#9E0142", alpha=0.5 )  +
  ylim(25, 35) + #xlim(2008, 2018) +
  scale_x_continuous(limits=c(2007.5, 2017), breaks=seq(2008,2016,by=2),  labels=c("","2010","2012","2014","2016")) +
  ylab(expression(paste("Maximum July temperature (°C)"))) + xlab("Year") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),axis.text=element_text(size=18),
                     axis.title=element_text(size=20),panel.border = element_rect(colour = "black", fill=NA, size=2),
                     plot.title = element_text(size = 26)) +
  ggtitle('A') + guides(fill=F) +
  scale_fill_manual(values=(c("#9E0142" , "#9E0142", "#FEF0A5", "#FEF0A5", "#5E4FA2", "#5E4FA2"))) +
  geom_segment(data=avJMT, mapping=aes(x=xs, y=meanJMTupperse, xend=xs, yend=meanJMTlowerse), 
               color="black", size=1.2, alpha=0.2) +
  geom_segment(data=avJMT, mapping=aes(x=xs, y=meanJMTupperse, xend=xs, yend=meanJMTlowerse), 
               color=rev(c("#9E0142" , "#9E0142", "#FEF0A5", "#FEF0A5", "#5E4FA2", "#5E4FA2")),
               size=1, alpha=0.2) +
  geom_point(data=avJMT, aes(x=xs, y=meanJMT), shape=21, size=6, alpha=0.3,
             fill=rev(c("#9E0142" , "#9E0142", "#FEF0A5", "#FEF0A5", "#5E4FA2", "#5E4FA2"))) +
  geom_vline(xintercept=2009) +
  annotate("text", label = c("N1", "N2", "C1", "C2", "S1", "S2"), y = avJMT$meanJMT, x = avJMT$xs, size = 3, colour = "black")

JMTlmSE





##### Look at variation in seasonality in the 7 recent years (2010-2017)
clim_study$Year2 <- clim_study$Year-2010
# ANCOVA
mod <- lm(seasonality ~ Year2*ID1, data=clim_study)
mod1 <- lm(seasonality ~ Year2+ID1, data=clim_study)
mod2 <- lm(seasonality ~ Year2, data=clim_study)
mod3 <- lm(seasonality ~ ID1+I(Year2^2)+Year2, data=clim_study)
MuMIn::AICc(mod, mod1, mod2, mod3)
summary(mod3) # bs for year are 5.826e+02/-1.446e-01, r2adj is 0.5143 
anova(mod3)

mod3 <- lm(seasonality ~ ID1+I(Year^2)+Year, data=clim_study)

# Now create predicted lines for each population
N1pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("N1",500))
N1pred$fit <- predict(mod3, newdata = N1pred, re.form = NA, type="response")
N2pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("N2",500))
N2pred$fit <- predict(mod3, newdata = N2pred, re.form = NA, type="response")
C1pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("C1",500))
C1pred$fit <- predict(mod3, newdata = C1pred, re.form = NA, type="response")
C2pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("C2",500))
C2pred$fit <- predict(mod3, newdata = C2pred, re.form = NA, type="response")
S1pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("S1",500))
S1pred$fit <- predict(mod3, newdata = S1pred, re.form = NA, type="response")
S2pred <- data.frame(Year = seq(from = min(clim_study$Year),to = max(clim_study$Year),length.out = 500),ID1 = rep("S2",500))
S2pred$fit <- predict(mod3, newdata = S2pred, re.form = NA, type="response")

avS$xs <- rep(2008,6)+c(-0.2,0.1,-0.1,0.2,0,0.3)
sds <- as.data.frame(clim_hist %>%
                       group_by(ID1) %>%
                       summarise(seasonality.sd = sd(seasonality)))
avS$sd <- rev(sds$seasonality.sd)
avS$meanSuppersd <- avS$meanS+avS$sd
avS$meanSlowersd <- avS$meanS-avS$sd
avS$se <- rev(sds$seasonality.sd)/sqrt(length(unique(clim_hist$Year)))
avS$meanSupperse <- avS$meanS+avS$se
avS$meanSlowerse <- avS$meanS-avS$se

SlmSE <- ggplot() +
  geom_point(data = clim_study, aes(y = seasonality, x = Year, fill = ID1), 
             alpha=0.5, shape=21, colour="black", size=4) +
  geom_line(data = N1pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5)  +
  geom_line(data = N2pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5)  +
  geom_line(data = C1pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5)  +
  geom_line(data = C2pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5)  +
  geom_line(data = S1pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5 )  +
  geom_line(data = S2pred, aes(x = Year, y = fit), size = 2, color="black", alpha=0.5 )  +
  geom_line(data = N1pred, aes(x = Year, y = fit), size = 1.5, color="#5E4FA2", alpha=0.5)  +
  geom_line(data = N2pred, aes(x = Year, y = fit), size = 1.5, color="#5E4FA2", alpha=0.5)  +
  geom_line(data = C1pred, aes(x = Year, y = fit), size = 1.5, color="#FEF0A5", alpha=0.5)  +
  geom_line(data = C2pred, aes(x = Year, y = fit), size = 1.5, color="#FEF0A5", alpha=0.5)  +
  geom_line(data = S1pred, aes(x = Year, y = fit), size = 1.5, color="#9E0142", alpha=0.5 )  +
  geom_line(data = S2pred, aes(x = Year, y = fit), size = 1.5, color="#9E0142", alpha=0.5 )  +
  ylim(26, 36) + #xlim(2008, 2018) +
  scale_x_continuous(limits=c(2007.5, 2017), breaks=seq(2008,2016,by=2),  labels=c("","2010","2012","2014","2016")) +
  ylab(expression(paste("Seasonality (°C)"))) + xlab("Year") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),axis.text=element_text(size=18),
                     axis.title=element_text(size=20),panel.border = element_rect(colour = "black", fill=NA, size=2),
                     plot.title = element_text(size = 26)) +
  ggtitle('B') + guides(fill=F) +
  scale_fill_manual(values=(c("#9E0142" , "#9E0142", "#FEF0A5", "#FEF0A5", "#5E4FA2", "#5E4FA2"))) +
  geom_segment(data=avS, mapping=aes(x=xs, y=meanSupperse, xend=xs, yend=meanSlowerse), 
               color="black", size=1.2, alpha=0.2) +
  geom_segment(data=avS, mapping=aes(x=xs, y=meanSupperse, xend=xs, yend=meanSlowerse), 
               color=rev(c("#9E0142" , "#9E0142", "#FEF0A5", "#FEF0A5", "#5E4FA2", "#5E4FA2")),
               size=1, alpha=0.2) +
  geom_point(data=avS, aes(x=xs, y=meanS), shape=21, size=6, alpha=0.3,
             fill=rev(c("#9E0142" , "#9E0142", "#FEF0A5", "#FEF0A5", "#5E4FA2", "#5E4FA2"))) +
  geom_vline(xintercept=2009) +
  annotate("text", label = c("N1", "N2", "C1", "C2", "S1", "S2"), y = avS$meanS, x = avS$xs, size = 3, colour = "black") 

SlmSE





# print to pdf: Figure S5
# shows recent trends in july max temp and intra-annual temp var
png("Figures/Figure S5.png", width=700,height=350)    
grid.arrange(JMTlmSE,SlmSE, ncol=2, nrow=1)
dev.off()








####################
# DO GROWTH CHAMBER TEMPS MATCH TEMPS IN THE FIELD?
####################
dat$Region <- as.character(dat$ID1)
dat$Region[which(dat$Region %in% c("N1", "N2"))] <- "North"
dat$Region[which(dat$Region %in% c("C1", "C2"))] <- "Central"
dat$Region[which(dat$Region %in% c("S1", "S2"))] <- "South"
dat$Region <- as.factor(dat$Region)
dat$Region <- factor(dat$Region, levels(dat$Region)[c(2,1,3)])
dat$ID1 <- factor(dat$ID1, levels(dat$ID1)[c(3,4,1,2,5,6)])



# Illustrate how monthly min and max temperatures of all months across 
# all populations compares to growth chamber temperatures (day and night)
# and facet by population
palpha <- 0.2
psize <- 0.5
Tcomp3 <- ggplot() +
  geom_point(data = dat, aes(y = Tmin01, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin02, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin03, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin04, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin05, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin06, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin07, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin08, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin09, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin10, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin11, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmin12, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax01, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax02, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax03, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax04, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax05, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax06, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax07, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax08, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax09, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax10, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax11, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_point(data = dat, aes(y = Tmax12, x = Year, fill = Region), alpha=palpha, shape=21, size=psize) + 
  geom_hline(yintercept=c(c(10,15,20,25,30,35,40,45)+0.2), linetype="dashed", color = "black", size=1, alpha=0.75) +
  geom_hline(yintercept=c(c(-5,0,5,10,15,20,25,30)-0.2), linetype="dashed", color = "gray", size=1, alpha=0.75) +
  ylim(-8, 46) + xlim(1900, 2017) +
  ylab(expression(paste("Min./Max. monthly temperatures (°C)"))) + xlab("Year") +
  theme_classic() +   
  theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=2), 
        panel.background = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", size = 0)) +
  scale_fill_manual(values=c("#5E4FA2", "#FEF0A5", "#9E0142")) +
  guides(fill = guide_legend( title.theme = element_text(size = 10))) +
  facet_wrap(~ID1, nrow=3)
Tcomp3




# print to pdf: FIG S3 FOR TPC MANUSCRIPT
quartz(height=6,width=4,type="pdf",dpi=600,file="Figures/Figure S1.pdf")
Tcomp3
dev.off()



