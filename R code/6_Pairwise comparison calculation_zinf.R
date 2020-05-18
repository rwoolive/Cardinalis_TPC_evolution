
#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Calculate pairwise comparisons between groups (Table 1)
#### DATE LAST MODIFIED: 2020-05-18 by rcw


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





# posterior draws
tidy_perf_groups_av <- read.csv("Analysis output/tidy_perf_groups_av_zinf.csv")[,-1]






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
###### than another is less than 0.025 or greater than 0.975, and two "**" if it 
###### above the level value specified  -- default is 0.99. The values in the 
###### table are the posterior average of the column group minus the row group.)
hi <- 0.99
lo <- 1 - hi

###### pairwise comparison function

pairwise_probs_groups_av <- function(stan_df, group, params, level = 0.99){
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
      file = str_glue("Analysis output/pairwise_tables_groups_av_zinf/{param}_pairwise_newthresh.csv") ### RCW altered this from ST's code that originally used the 0.95/0.05 thresholds (original file did not have the "_newthresh")
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
N1comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==2)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==1)]
length(which(N1comp==TRUE))/length(N1comp)
## then, the average difference with 95% CI
N1compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==2)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==1)]
mN1 <- mean(N1compNum)
qN1 <- quantile(N1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for N2 2017 vs. N2 2010
N2comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==4)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==3)]
length(which(N2comp==TRUE))/length(N2comp)
## then, the average difference with 95% CI
N2compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==4)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==3)]
mN2 <- mean(N2compNum)
qN2 <- quantile(N2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for C1 2017 vs. C1 2010
C1comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==6)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==5)]
length(which(C1comp==TRUE))/length(C1comp)
## then, the average difference with 95% CI
C1compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==6)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==5)]
mC1 <- mean(C1compNum)
qC1 <- quantile(C1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for C2 2017 vs. C2 2010
C2comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==8)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==7)]
length(which(C2comp==TRUE))/length(C2comp)
## then, the average difference with 95% CI
C2compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==8)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==7)]
mC2 <- mean(C2compNum)
qC2 <- quantile(C2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for S1 2017 vs. S1 2010
S1comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==10)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==9)]
length(which(S1comp==TRUE))/length(S1comp)
## then, the average difference with 95% CI
S1compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==10)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==9)]
mS1 <- mean(S1compNum)
qS1 <- quantile(S1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater topt for S2 2017 vs. S2 2010
S2comp <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==12)] > tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==11)]
length(which(S2comp==TRUE))/length(S2comp)
## then, the average difference with 95% CI
S2compNum <- tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==12)] - tidy_perf_groups_av$maximaBT[which(tidy_perf_groups_av$species==11)]
mS2 <- mean(S2compNum)
qS2 <- quantile(S2compNum, probs = c(0.025,0.975))

# combine into a dataframe
means <- round(c(mN1, mN2, mC1, mC2, mS1, mS2),3)
lowers <- as.vector(round(c(qN1[1], qN2[1], qC1[1], qC2[1], qS1[1], qS2[1]),3))
uppers <- as.vector(round(c(qN1[2], qN2[2], qC1[2], qC2[2], qS1[2], qS2[2]),3))
pwctopt <- paste(means, " [", lowers, ", ", uppers, "]", sep="")



# calculate how many iterations estimated greater b50 for N1 2017 vs. N1 2010
N1comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==2)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==1)]
length(which(N1comp==TRUE))/length(N1comp)
## then, the average difference with 95% CI
N1compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==2)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==1)]
mN1 <- mean(N1compNum)
qN1 <- quantile(N1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for N2 2017 vs. N2 2010
N2comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==4)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==3)]
length(which(N2comp==TRUE))/length(N2comp)
## then, the average difference with 95% CI
N2compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==4)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==3)]
mN2 <- mean(N2compNum)
qN2 <- quantile(N2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for C1 2017 vs. C1 2010
C1comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==6)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==5)]
length(which(C1comp==TRUE))/length(C1comp)
## then, the average difference with 95% CI
C1compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==6)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==5)]
mC1 <- mean(C1compNum)
qC1 <- quantile(C1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for C2 2017 vs. C2 2010
C2comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==8)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==7)]
length(which(C2comp==TRUE))/length(C2comp)
## then, the average difference with 95% CI
C2compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==8)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==7)]
mC2 <- mean(C2compNum)
qC2 <- quantile(C2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for S1 2017 vs. S1 2010
S1comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==10)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==9)]
length(which(S1comp==TRUE))/length(S1comp)
## then, the average difference with 95% CI
S1compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==10)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==9)]
mS1 <- mean(S1compNum)
qS1 <- quantile(S1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater b50 for S2 2017 vs. S2 2010
S2comp <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==12)] > tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==11)]
length(which(S2comp==TRUE))/length(S2comp)
## then, the average difference with 95% CI
S2compNum <- tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==12)] - tidy_perf_groups_av$B50[which(tidy_perf_groups_av$species==11)]
mS2 <- mean(S2compNum)
qS2 <- quantile(S2compNum, probs = c(0.025,0.975))

# combine into a dataframe
means <- round(c(mN1, mN2, mC1, mC2, mS1, mS2),3)
lowers <- as.vector(round(c(qN1[1], qN2[1], qC1[1], qC2[1], qS1[1], qS2[1]),3))
uppers <- as.vector(round(c(qN1[2], qN2[2], qC1[2], qC2[2], qS1[2], qS2[2]),3))
pwcb50 <- paste(means, " [", lowers, ", ", uppers, "]", sep="")



# calculate how many iterations estimated greater b80 for N1 2017 vs. N1 2010
N1comp <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==2)] > tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==1)]
length(which(N1comp==TRUE))/length(N1comp)
## then, the average difference with 95% CI
N1compNum <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==2)] - tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==1)]
mN1 <- mean(N1compNum)
qN1 <- quantile(N1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater B80 for N2 2017 vs. N2 2010
N2comp <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==4)] > tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==3)]
length(which(N2comp==TRUE))/length(N2comp)
## then, the average difference with 95% CI
N2compNum <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==4)] - tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==3)]
mN2 <- mean(N2compNum)
qN2 <- quantile(N2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater B80 for C1 2017 vs. C1 2010
C1comp <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==6)] > tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==5)]
length(which(C1comp==TRUE))/length(C1comp)
## then, the average difference with 95% CI
C1compNum <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==6)] - tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==5)]
mC1 <- mean(C1compNum)
qC1 <- quantile(C1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater B80 for C2 2017 vs. C2 2010
C2comp <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==8)] > tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==7)]
length(which(C2comp==TRUE))/length(C2comp)
## then, the average difference with 95% CI
C2compNum <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==8)] - tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==7)]
mC2 <- mean(C2compNum)
qC2 <- quantile(C2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater B80 for S1 2017 vs. S1 2010
S1comp <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==10)] > tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==9)]
length(which(S1comp==TRUE))/length(S1comp)
## then, the average difference with 95% CI
S1compNum <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==10)] - tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==9)]
mS1 <- mean(S1compNum)
qS1 <- quantile(S1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater B80 for S2 2017 vs. S2 2010
S2comp <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==12)] > tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==11)]
length(which(S2comp==TRUE))/length(S2comp)
## then, the average difference with 95% CI
S2compNum <- tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==12)] - tidy_perf_groups_av$B80[which(tidy_perf_groups_av$species==11)]
mS2 <- mean(S2compNum)
qS2 <- quantile(S2compNum, probs = c(0.025,0.975))

# combine into a dataframe
means <- round(c(mN1, mN2, mC1, mC2, mS1, mS2),3)
lowers <- as.vector(round(c(qN1[1], qN2[1], qC1[1], qC2[1], qS1[1], qS2[1]),3))
uppers <- as.vector(round(c(qN1[2], qN2[2], qC1[2], qC2[2], qS1[2], qS2[2]),3))
pwcb80 <- paste(means, " [", lowers, ", ", uppers, "]", sep="")



# calculate how many iterations estimated greater max_RGR for N1 2017 vs. N1 2010
N1comp <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==2)] > tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==1)]
length(which(N1comp==TRUE))/length(N1comp)
## then, the average difference with 95% CI
N1compNum <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==2)] - tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==1)]
mN1 <- mean(N1compNum)
qN1 <- quantile(N1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater max_RGR for N2 2017 vs. N2 2010
N2comp <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==4)] > tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==3)]
length(which(N2comp==TRUE))/length(N2comp)
## then, the average difference with 95% CI
N2compNum <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==4)] - tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==3)]
mN2 <- mean(N2compNum)
qN2 <- quantile(N2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater max_RGR for C1 2017 vs. C1 2010
C1comp <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==6)] > tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==5)]
length(which(C1comp==TRUE))/length(C1comp)
## then, the average difference with 95% CI
C1compNum <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==6)] - tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==5)]
mC1 <- mean(C1compNum)
qC1 <- quantile(C1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater max_RGR for C2 2017 vs. C2 2010
C2comp <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==8)] > tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==7)]
length(which(C2comp==TRUE))/length(C2comp)
## then, the average difference with 95% CI
C2compNum <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==8)] - tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==7)]
mC2 <- mean(C2compNum)
qC2 <- quantile(C2compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater max_RGR for S1 2017 vs. S1 2010
S1comp <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==10)] > tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==9)]
length(which(S1comp==TRUE))/length(S1comp)
## then, the average difference with 95% CI
S1compNum <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==10)] - tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==9)]
mS1 <- mean(S1compNum)
qS1 <- quantile(S1compNum, probs = c(0.025,0.975))

# calculate how many iterations estimated greater max_RGR for S2 2017 vs. S2 2010
S2comp <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==12)] > tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==11)]
length(which(S2comp==TRUE))/length(S2comp)
## then, the average difference with 95% CI
S2compNum <- tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==12)] - tidy_perf_groups_av$max_RGR[which(tidy_perf_groups_av$species==11)]
mS2 <- mean(S2compNum)
qS2 <- quantile(S2compNum, probs = c(0.025,0.975))

# combine into a dataframe
means <- round(c(mN1, mN2, mC1, mC2, mS1, mS2),3)
lowers <- as.vector(round(c(qN1[1], qN2[1], qC1[1], qC2[1], qS1[1], qS2[1]),3))
uppers <- as.vector(round(c(qN1[2], qN2[2], qC1[2], qC2[2], qS1[2], qS2[2]),3))
pwcmax_RGR <- paste(means, " [", lowers, ", ", uppers, "]", sep="")



# now combine pwc's into a dataframe and export to csv
pwc <- data.frame(pwctopt=pwctopt,
                  pwcb50=pwcb50,
                  pwcb80=pwcb80,
                  pwcmax_RGR=pwcmax_RGR)
write.csv(pwc, "Analysis output/Table 1 tpc_zinf.csv")

