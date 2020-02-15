#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Determine if TPC parameters have evolved, and whether they are 
#### correlated with climatic data 
#### AUTHOR: Rachel Wooliver 

library(ggplot2)
library(tidyverse)
theme_set(theme_bw())


###########################
# read in the data
###########################
popBayes_params <- read.csv("Analysis output/group_mean_params_av.csv") # these are the parameters for the Bayes model including family-averaged data
avJMT <- read.csv("Processed data/average JMT by pop.csv") # population-averaged july max temp
avS <- read.csv("Processed data/average seasonality by pop.csv") # population-averaged seasonality



###########################
# add JMT and seasonality means and anomalies to tpc parameter dataset
###########################
pop_params <- popBayes_params
pop_params$meanJMT <- rep(avJMT$meanJMT, 2)
pop_params$devJMT <- rep(avJMT$avDevJMT, 2)
pop_params$meanS <- rep(avS$meanS, 2)
pop_params$devS <- rep(avS$avDevS, 2)


###########################
# add region and make factors
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
palpha <- 0.65
psize <- 8
tsize <- 4
############################
# is topt related to historical July max temps? 
# test of divergence across space
############################

# model of topt across historical July max temps
# no effects
mod <- lm(maximaBT~meanJMT*year, data=pop_params)
anova(mod)
round(summary(mod)$adj.r.squared,digits=3)
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3)
shapiro.test(residuals(mod)) # normally distributed residuals

# remove year from model
mod <- lm(maximaBT~meanJMT, data=pop_params)
anova(mod) # df=1,10; ss=0.20412; msq=0.20412;  f=0.9029; p=0.3644
round(summary(mod)$adj.r.squared,digits=3) # rsqadj=-0.009
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3) # b=0.063
shapiro.test(residuals(mod)) # normally distributed residuals

# remove temp from model
mod <- lm(maximaBT~year, data=pop_params)
anova(mod) # df=1,10; ss=0.06488; msq=0.06488;  f=0.2703; p=0.6144
round(summary(mod)$adj.r.squared,digits=3) # rsqadj=-0.071
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3) # b=0.147
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
  annotate("text", label = c("North", "Central", "South"), y = 34.15, x = c(28.1,30.1,32.1), size = 5, colour = "black", hjust = 0) +
  geom_point(aes(y=34.15, x=27.75), shape=22, size=8, fill="#5E4FA2") +
  geom_point(aes(y=34.15, x=29.75), shape=22, size=8, fill="#FEF0A5") +
  geom_point(aes(y=34.15, x=31.75), shape=22, size=8, fill="#9E0142") +
  annotate("text", label = c("Ancestors", "Descendants"), y = 33.85, x = c(28.1,30.6), size = 5, colour = "black", hjust = 0) +
  geom_point(aes(y=33.85, x=27.75), shape=shapes[1], size=8) +
  geom_point(aes(y=33.85, x=30.25), shape=shapes[2], size=8) +
  xlim(c(26.75,34.05)) + ylim(c(31.67,34.2)) +
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
# B50 tends to increase with seasonality
mod <- lm(B50~meanS*year, data=pop_params)
anova(mod) 
# remove cohort
mod <- lm(B50~meanS, data=pop_params)
anova(mod) # df=1,10; ss=2.7252 ; msq=2.72523  ;  f=7.9093  ; p=0.0184
round(summary(mod)$adj.r.squared,digits=3) # rsqadj=0.386
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3) # b=0.253
shapiro.test(residuals(mod)) # normally distributed residuals
# remove seasonality
mod <- lm(B50~year, data=pop_params)
anova(mod) # df=1,10; ss=0.1564  ; msq=0.15643    ;  f=0.2601   ; p=0.6211
round(summary(mod)$adj.r.squared,digits=3) # rsqadj= -0.072
round(as.data.frame(summary(mod)$coefficients)$Estimate[2],digits=3) # b= -0.228
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
library(ggpubr)
theme_set(theme_pubr())
figure <- ggarrange(plot1Bayes, plot2bBayes, ncol=2, nrow=1, labels = c("A","B"), font.label = list(size=18))
pdf("Figures/Figure 4.pdf", height=4.5, width=10)
figure
dev.off()


