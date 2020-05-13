
#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Determine if TPC parameters are 
####          correlated with climatic data 
#### AUTHOR: Rachel Wooliver
#### DATE LAST MODIFIED: 2020-05-12 by rcw


###########################
# load packages and specify settings
###########################
library(ggplot2)
library(tidyverse)
theme_set(theme_bw())


###########################
# read in the data
###########################
popBayes_params <- read.csv("Analysis output/group_mean_params_av_zinf.csv") # these are the parameters for the group-level Bayes model of family-averaged data
avJMT <- read.csv("Processed data/average JMT by pop.csv") # population-averaged july max temp with recent anomalies
avS <- read.csv("Processed data/average seasonality by pop.csv") # population-averaged seasonality with recent anomalies



###########################
# add JMT and seasonality means and mean anomalies to tpc parameter dataset
###########################
pop_params <- popBayes_params
pop_params$meanJMT <- rep(avJMT$meanJMT, 2)
pop_params$devJMT <- rep(avJMT$avDevJMT, 2)
pop_params$meanS <- rep(avS$meanS, 2)
pop_params$devS <- rep(avS$avDevS, 2)


###########################
# add region and redo region, year, and population columns as factors
###########################
pop_params$reg <- rep(rep(c("N", "C", "S"), each=2),2)
pop_params$reg <- as.factor(pop_params$reg)
pop_params$reg <- factor(pop_params$reg, levels(pop_params$reg)[c(2,1,3)])
pop_params$year <- as.factor(pop_params$year)
pop_params$Pop <- as.factor(pop_params$Pop)
pop_params$Pop <- factor(pop_params$Pop, levels(pop_params$Pop)[rev(c(3,4,1,2,5,6))])


###########################
# shapes for years
# 2010=circle, 2017=triangle
###########################
shapes <- c(21,24)


###########################
# create datasets by year
###########################
dat2010 <- pop_params[which(pop_params$year==2010),]
dat2017 <- pop_params[which(pop_params$year==2017),]







##################################################################
##### Fig 3 of TPC manuscript
##################################################################
palpha <- 0.65 # alpha of points
psize <- 8 # size of points
tsize <- 4 # size of population text
############################
# is topt related to historical July max temps? 
# test of divergence across space
############################

# model of topt across historical July max temps
# no effects of JMT, year, or their interaction
mod <- lm(maximaBT~meanJMT*year, data=pop_params)
anova(mod)
shapiro.test(residuals(mod)) # normally distributed residuals

# remove year from model
mod <- lm(maximaBT~meanJMT, data=pop_params)
anova(mod) # df=1,10; ss=0.47142 ; msq=0.47142  ;  f=1.9989 ; p=0.1878
round(summary(mod)$adj.r.squared,digits=3) # rsqadj=0.083
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3) # b=0.095
shapiro.test(residuals(mod)) # normally distributed residuals

# remove temp from model
mod <- lm(maximaBT~year, data=pop_params)
anova(mod) # df=1,10; ss=0.10335; msq=0.10335;  f=0.379; p=0.5519
round(summary(mod)$adj.r.squared,digits=3) # rsqadj=-0.06
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3) # b=0.186
shapiro.test(residuals(mod)) # normally distributed residuals



# plot
plot1Bayes <- ggplot(dat2010, aes(x=meanJMT, y=maximaBT)) + 
  geom_smooth(data=pop_params, method="lm", se=T, color="gray", fill="gray", size=1, alpha=0.4, linetype=2) + 
  geom_point(data=dat2017, aes(fill=reg), shape=shapes[2], color="black", size=psize, alpha=palpha) +
  geom_point(aes(fill=reg, shape=year), color="black", size=psize, alpha=palpha) + 
  guides(size = FALSE, fill=FALSE, shape=FALSE) +
  scale_fill_manual(values=c("#5E4FA2", "#FEF0A5", "#9E0142")) +
  scale_shape_manual(values=c(shapes)[1]) +
  labs(y=expression(paste("Thermal optimum (°C)")),
       x=expression(paste("Maximum July temperature (°C)"))) +
  annotate("text", label = c("North", "Central", "South"), y = 34.25, x = c(28.1,30.1,32.1), size = 5, colour = "black", hjust = 0) +
  geom_point(aes(y=34.25, x=27.75), shape=22, size=8, fill="#5E4FA2") +
  geom_point(aes(y=34.25, x=29.75), shape=22, size=8, fill="#FEF0A5") +
  geom_point(aes(y=34.25, x=31.75), shape=22, size=8, fill="#9E0142") +
  annotate("text", label = c("Ancestors", "Descendants"), y = 34, x = c(28.1,30.6), size = 5, colour = "black", hjust = 0) +
  geom_point(aes(y=34, x=27.75), shape=shapes[1], size=8) +
  geom_point(aes(y=34, x=30.25), shape=shapes[2], size=8) +
  xlim(c(26.75,34.05)) + ylim(c(32,34.3)) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16),
        legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        strip.background = element_rect(colour=NA, fill=NA),
        strip.text.x = element_text(size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate("text", label = c("N1", "N2", "C1", "C2", "S1", "S2"), y = dat2010$maximaBT, x = dat2010$meanJMT, size = tsize, colour = "black")
  
  
  plot1Bayes
  
  
  


############################
# is B50 related to seasonality? 
# test of divergence across space
############################

# model of B50 across seasonality
# B50 increases with seasonality
# no effect of year or interaction between seasonality and year
mod <- lm(B50~meanS*year, data=pop_params)
anova(mod) 

# remove cohort
mod <- lm(B50~meanS, data=pop_params)
anova(mod) # df=1,10; ss=1.6625 ; msq=1.6625  ;  f=5.1997  ; p=0.04576
round(summary(mod)$adj.r.squared,digits=3) # rsqadj=0.276
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3) # b=0.198
shapiro.test(residuals(mod)) # normally distributed residuals

# remove seasonality
mod <- lm(B50~year, data=pop_params)
anova(mod) # df=1,10; ss=0.3851  ; msq=0.3851    ;  f=0.8607   ; p=0.3754
round(summary(mod)$adj.r.squared,digits=3) # rsqadj= -0.013
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3) # b= -0.358
shapiro.test(residuals(mod)) # normally distributed residuals


# plot
arrs <- data.frame(x=dat2010$meanS,
                   xend=dat2017$meanS,
                   y=dat2010$B50 + c(0,0,0,0,-0.13,0),
                   yend=dat2017$B50 + c(0,0,0,0,0.21,0))

plot2bBayes <- ggplot(dat2010, aes(x=meanS, y=B50)) + 
  geom_smooth(data=pop_params, method="lm", se=T, color="gray", fill="gray", size=1, alpha=0.4, linetype=1) +
  geom_point(data=dat2017, aes(fill=reg), shape=shapes[2], color="black", size=psize, alpha=palpha) +
  geom_segment(data = arrs[5,], aes(x=x, xend=xend, y=y, yend=yend), size = 1, color="black",
               arrow = arrow(length = unit(0.25, "cm"))) +
  xlim(c(28,34)) +
  ylim(c(21,24.2)) +
  geom_point(aes(fill=reg, shape=year), color="black", size=psize, alpha=palpha) + 
  guides(size = FALSE, fill=FALSE, shape=FALSE) +
  scale_fill_manual(values=c("#5E4FA2", "#FEF0A5", "#9E0142")) +
  scale_shape_manual(values=c(shapes)[1]) +
  labs(y=expression(paste("Breadth (°C)")),
       x=expression(paste("Seasonality (°C)"))) +
  theme(axis.text=element_text(size=14),
        axis.title.y=element_text(size=16),
        axis.title.x = element_text(size=16),
        #axis.text.y=element_blank(),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        strip.background = element_rect(colour=NA, fill=NA),
        strip.text.x = element_text(size=20, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate("text", label = c("N1", "N2", "C1", "C2", "S1", "S2"), y = dat2010$B50, x = dat2010$meanS, size = tsize, colour = "black")


plot2bBayes





# plot climate vs. tpc parameter regressions 
# export tpc plot
ggsave("Figures/Figure 4_zinf.png", 
       ggpubr::ggarrange(plot1Bayes, plot2bBayes, ncol=2, nrow=1, labels = c("A","B"), font.label = list(size=18)), 
       height=4.5, width=10, dpi=600)

figure <- ggpubr::ggarrange(plot1Bayes, plot2bBayes, ncol=2, nrow=1, labels = c("A","B"), font.label = list(size=18))
pdf("Figures/Figure 4_zinf.pdf", height=4.5, width=10)
figure
dev.off()


