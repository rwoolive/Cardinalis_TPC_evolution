#### PROJECT: Mimulus cardinalis TPC project
#### PURPOSE: Plot pairwise comparisons of TPC parameters 
#### between 12 M. cardinalis groups.
#### AUTHOR: Rachel Wooliver


library(reshape2)
library(ggplot2)


# read in data, where values represent posterior averages of the column group minus 
# the row group, * denotes probability of one group having a smaller parameter
# value than another is less than 0.05 or greater than 0.95, and ** denotes
# probability of one group having a smaller parameter value than anoter is less
# than 0.01 or greater than 0.99  

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for comparisons nonessential for hypothesis testing and 1's for comparisons 
# essential for hypothesis testing
boxes <- matrix(0, nrow=12, ncol=12)
boxes[1,2] <- 1
boxes[3,4] <- 1
boxes[5,6] <- 1
boxes[7,8] <- 1
boxes[9,10] <- 1
boxes[11,12] <- 1
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri_boxes <- get_upper_tri(boxes)
melted_cormat_boxes <- melt(upper_tri_boxes, na.rm = TRUE)
melted_cormat_boxes$value[which(melted_cormat_boxes$value==0)] <- NA
melted_cormat_boxes$value[which(melted_cormat_boxes$value==1)] <- "[  ]"




##### 
##### thermal optima (same results with raw and back-transformed values):
##### 

dat <- read.csv("Analysis output/pairwise_tables_groups_av/maximaBT_pairwise_newthresh.csv", 
                row.names=1)[c(1,7,2,8,3,9,4,10,5,11,6,12),c(1,7,2,8,3,9,4,10,5,11,6,12)] 
colnames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
rownames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
dat <- dat[12:1,12:1]

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for non-significant comparisons and 1's for significant comparisons
aster <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(aster)[1]){
  vals <- dat[,i]
  for(j in 1:dim(aster)[1]){
    val <- as.numeric(as.character(vals[j]))
    if(is.na(val)==FALSE){aster[j,i] <- 0} else{aster[j,i] <- 1}
  }
}



# create a matrix the same size as the pairwise comparison table, but containing 
# numeric values (no asterisks)
nums <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(nums)[1]){
  nums[,i] <- as.numeric(stringr::str_remove_all(as.character(dat[,i]), "[*]")) * (-1)
}
rownames(nums) <- rownames(dat)
colnames(nums) <- colnames(dat)



# plot heatmap with asterisks for significant comparisons
# datasets: aster (0=nonsignif/1=signif matrix), and nums (values of dat without asterisks)
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(nums)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Now print asterisks onto heatmap. First get the melted upper triangle of aster
upper_tri_aster <- get_upper_tri(aster)
melted_cormat_aster <- melt(upper_tri_aster, na.rm = TRUE)
melted_cormat_aster$value[which(melted_cormat_aster$value==0)] <- NA
melted_cormat_aster$value[which(melted_cormat_aster$value==1)] <- "*"
melted_cormat$aster <- melted_cormat_aster$value

# put boxes into melted_cormat
melted_cormat$boxes <- melted_cormat_boxes$value


# Set limits of values
limits <- c(min(melted_cormat$value),max(melted_cormat$value))

# Heatmap, where red means y group is greater than x group
ggheatmapTopt <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(limits), 
                       space = "Lab", 
                       name="Pairwise difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label=aster), color = "black", shape=8) +
  geom_text(aes(Var2, Var1, label=boxes), color = "black", size=5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))







##### 
##### x_min:
##### 

dat <- read.csv("Analysis output/pairwise_tables_groups_av/x_minBT_pairwise_newthresh.csv", 
                row.names=1)[c(1,7,2,8,3,9,4,10,5,11,6,12),c(1,7,2,8,3,9,4,10,5,11,6,12)] 
colnames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
rownames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
dat <- dat[12:1,12:1]

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for non-significant comparisons and 1's for significant comparisons
aster <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(aster)[1]){
  vals <- dat[,i]
  for(j in 1:dim(aster)[1]){
    val <- as.numeric(as.character(vals[j]))
    if(is.na(val)==FALSE){aster[j,i] <- 0} else{aster[j,i] <- 1}
  }
}

# create a matrix the same size as the pairwise comparison table, but containing 
# numeric values (no asterisks)
nums <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(nums)[1]){
  nums[,i] <- as.numeric(stringr::str_remove_all(as.character(dat[,i]), "[*]")) * (-1)
}
rownames(nums) <- rownames(dat)
colnames(nums) <- colnames(dat)



# plot heatmap with asterisks for significant comparisons
# datasets: aster (0=nonsignif/1=signif matrix), and nums (values of dat without asterisks)
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(nums)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Now print asterisks onto heatmap. First get the melted upper triangle of aster
upper_tri_aster <- get_upper_tri(aster)
melted_cormat_aster <- melt(upper_tri_aster, na.rm = TRUE)
melted_cormat_aster$value[which(melted_cormat_aster$value==0)] <- NA
melted_cormat_aster$value[which(melted_cormat_aster$value==1)] <- "*"
melted_cormat$aster <- melted_cormat_aster$value

# put boxes into melted_cormat
melted_cormat$boxes <- melted_cormat_boxes$value

# Set limits of values
limits <- c(min(melted_cormat$value),max(melted_cormat$value))

# Heatmap, where red means y group is greater than x group
ggheatmapXmin <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(limits), 
                       space = "Lab", 
                       name="Pairwise difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label=aster), color = "black", shape=8) +
  geom_text(aes(Var2, Var1, label=boxes), color = "black", size=5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))







##### 
##### x_max:
##### 

dat <- read.csv("Analysis output/pairwise_tables_groups_av/x_maxBT_pairwise_newthresh.csv", 
                row.names=1)[c(1,7,2,8,3,9,4,10,5,11,6,12),c(1,7,2,8,3,9,4,10,5,11,6,12)] 
colnames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
rownames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
dat <- dat[12:1,12:1]

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for non-significant comparisons and 1's for significant comparisons
aster <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(aster)[1]){
  vals <- dat[,i]
  for(j in 1:dim(aster)[1]){
    val <- as.numeric(as.character(vals[j]))
    if(is.na(val)==FALSE){aster[j,i] <- 0} else{aster[j,i] <- 1}
  }
}

# create a matrix the same size as the pairwise comparison table, but containing 
# numeric values (no asterisks)
nums <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(nums)[1]){
  nums[,i] <- as.numeric(stringr::str_remove_all(as.character(dat[,i]), "[*]")) * (-1)
}
rownames(nums) <- rownames(dat)
colnames(nums) <- colnames(dat)



# plot heatmap with asterisks for significant comparisons
# datasets: aster (0=nonsignif/1=signif matrix), and nums (values of dat without asterisks)
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(nums)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Now print asterisks onto heatmap. First get the melted upper triangle of aster
upper_tri_aster <- get_upper_tri(aster)
melted_cormat_aster <- melt(upper_tri_aster, na.rm = TRUE)
melted_cormat_aster$value[which(melted_cormat_aster$value==0)] <- NA
melted_cormat_aster$value[which(melted_cormat_aster$value==1)] <- "*"
melted_cormat$aster <- melted_cormat_aster$value

# put boxes into melted_cormat
melted_cormat$boxes <- melted_cormat_boxes$value

# Set limits of values
limits <- c(min(melted_cormat$value),max(melted_cormat$value))

# Heatmap, where red means y group is greater than x group
ggheatmapXmax <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(limits), 
                       space = "Lab", 
                       name="Pairwise difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label=aster), color = "black", shape=8) +
  geom_text(aes(Var2, Var1, label=boxes), color = "black", size=5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))














##### 
##### breadth (B50):
##### 

dat <- read.csv("Analysis output/pairwise_tables_groups_av/B50_pairwise_newthresh.csv", 
                row.names=1)[c(1,7,2,8,3,9,4,10,5,11,6,12),c(1,7,2,8,3,9,4,10,5,11,6,12)] 
colnames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
rownames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
dat <- dat[12:1,12:1]

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for non-significant comparisons and 1's for significant comparisons
aster <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(aster)[1]){
  vals <- dat[,i]
  for(j in 1:dim(aster)[1]){
    val <- as.numeric(as.character(vals[j]))
    if(is.na(val)==FALSE){aster[j,i] <- 0} else{aster[j,i] <- 1}
  }
}

# create a matrix the same size as the pairwise comparison table, but containing 
# numeric values (no asterisks)
nums <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(nums)[1]){
  nums[,i] <- as.numeric(stringr::str_remove_all(as.character(dat[,i]), "[*]")) * (-1)
}
rownames(nums) <- rownames(dat)
colnames(nums) <- colnames(dat)



# plot heatmap with asterisks for significant comparisons
# datasets: aster (0=nonsignif/1=signif matrix), and nums (values of dat without asterisks)
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(nums)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Now print asterisks onto heatmap. First get the melted upper triangle of aster
upper_tri_aster <- get_upper_tri(aster)
melted_cormat_aster <- melt(upper_tri_aster, na.rm = TRUE)
melted_cormat_aster$value[which(melted_cormat_aster$value==0)] <- NA
melted_cormat_aster$value[which(melted_cormat_aster$value==1)] <- "*"
melted_cormat$aster <- melted_cormat_aster$value

# put boxes into melted_cormat
melted_cormat$boxes <- melted_cormat_boxes$value

# Set limits of values
limits <- c(min(melted_cormat$value),max(melted_cormat$value))

# Heatmap, where red means y group is greater than x group
ggheatmapB50 <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(limits), 
                       space = "Lab", 
                       name="Pairwise difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label=aster), color = "black", shape=8) +
  geom_text(aes(Var2, Var1, label=boxes), color = "black", size=5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))












##### 
##### breadth (B80):
##### 

dat <- read.csv("Analysis output/pairwise_tables_groups_av/B80_pairwise_newthresh.csv", 
                row.names=1)[c(1,7,2,8,3,9,4,10,5,11,6,12),c(1,7,2,8,3,9,4,10,5,11,6,12)] 
colnames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
rownames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
dat <- dat[12:1,12:1]

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for non-significant comparisons and 1's for significant comparisons
aster <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(aster)[1]){
  vals <- dat[,i]
  for(j in 1:dim(aster)[1]){
    val <- as.numeric(as.character(vals[j]))
    if(is.na(val)==FALSE){aster[j,i] <- 0} else{aster[j,i] <- 1}
  }
}

# create a matrix the same size as the pairwise comparison table, but containing 
# numeric values (no asterisks)
nums <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(nums)[1]){
  nums[,i] <- as.numeric(stringr::str_remove_all(as.character(dat[,i]), "[*]")) * (-1)
}
rownames(nums) <- rownames(dat)
colnames(nums) <- colnames(dat)



# plot heatmap with asterisks for significant comparisons
# datasets: aster (0=nonsignif/1=signif matrix), and nums (values of dat without asterisks)
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(nums)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Now print asterisks onto heatmap. First get the melted upper triangle of aster
upper_tri_aster <- get_upper_tri(aster)
melted_cormat_aster <- melt(upper_tri_aster, na.rm = TRUE)
melted_cormat_aster$value[which(melted_cormat_aster$value==0)] <- NA
melted_cormat_aster$value[which(melted_cormat_aster$value==1)] <- "*"
melted_cormat$aster <- melted_cormat_aster$value

# put boxes into melted_cormat
melted_cormat$boxes <- melted_cormat_boxes$value

# Set limits of values
limits <- c(min(melted_cormat$value),max(melted_cormat$value))

# Heatmap, where red means y group is greater than x group
ggheatmapB80 <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(limits), 
                       space = "Lab", 
                       name="Pairwise difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label=aster), color = "black", shape=8) +
  geom_text(aes(Var2, Var1, label=boxes), color = "black", size=5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))







##### 
##### performance maximum:
##### 

dat <- read.csv("Analysis output/pairwise_tables_groups_av/max_RGR_pairwise_newthresh.csv", 
                row.names=1)[c(1,7,2,8,3,9,4,10,5,11,6,12),c(1,7,2,8,3,9,4,10,5,11,6,12)] 
colnames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
rownames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
dat <- dat[12:1,12:1]

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for non-significant comparisons and 1's for significant comparisons
aster <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(aster)[1]){
  vals <- dat[,i]
  for(j in 1:dim(aster)[1]){
    val <- as.numeric(as.character(vals[j]))
    if(is.na(val)==FALSE){aster[j,i] <- 0} else{aster[j,i] <- 1}
  }
}

# create a matrix the same size as the pairwise comparison table, but containing 
# numeric values (no asterisks)
nums <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(nums)[1]){
  nums[,i] <- as.numeric(stringr::str_remove_all(as.character(dat[,i]), "[*]")) * (-1)
}
rownames(nums) <- rownames(dat)
colnames(nums) <- colnames(dat)



# plot heatmap with asterisks for significant comparisons
# datasets: aster (0=nonsignif/1=signif matrix), and nums (values of dat without asterisks)
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(nums)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Now print asterisks onto heatmap. First get the melted upper triangle of aster
upper_tri_aster <- get_upper_tri(aster)
melted_cormat_aster <- melt(upper_tri_aster, na.rm = TRUE)
melted_cormat_aster$value[which(melted_cormat_aster$value==0)] <- NA
melted_cormat_aster$value[which(melted_cormat_aster$value==1)] <- "*"
melted_cormat$aster <- melted_cormat_aster$value

# put boxes into melted_cormat
melted_cormat$boxes <- melted_cormat_boxes$value

# Set limits of values
limits <- c(min(melted_cormat$value),max(melted_cormat$value))

# Heatmap, where red means y group is greater than x group
ggheatmapmax_RGR <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(limits), 
                       space = "Lab", 
                       name="Pairwise difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label=aster), color = "black", shape=8) +
  geom_text(aes(Var2, Var1, label=boxes), color = "black", size=5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))







##### 
##### area:
##### 

dat <- read.csv("Analysis output/pairwise_tables_groups_av/area_pairwise_newthresh.csv", 
                row.names=1)[c(1,7,2,8,3,9,4,10,5,11,6,12),c(1,7,2,8,3,9,4,10,5,11,6,12)] 
colnames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
rownames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
dat <- dat[12:1,12:1]

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for non-significant comparisons and 1's for significant comparisons
aster <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(aster)[1]){
  vals <- dat[,i]
  for(j in 1:dim(aster)[1]){
    val <- as.numeric(as.character(vals[j]))
    if(is.na(val)==FALSE){aster[j,i] <- 0} else{aster[j,i] <- 1}
  }
}

# create a matrix the same size as the pairwise comparison table, but containing 
# numeric values (no asterisks)
nums <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(nums)[1]){
  nums[,i] <- as.numeric(stringr::str_remove_all(as.character(dat[,i]), "[*]")) * (-1)
}
rownames(nums) <- rownames(dat)
colnames(nums) <- colnames(dat)



# plot heatmap with asterisks for significant comparisons
# datasets: aster (0=nonsignif/1=signif matrix), and nums (values of dat without asterisks)
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(nums)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Now print asterisks onto heatmap. First get the melted upper triangle of aster
upper_tri_aster <- get_upper_tri(aster)
melted_cormat_aster <- melt(upper_tri_aster, na.rm = TRUE)
melted_cormat_aster$value[which(melted_cormat_aster$value==0)] <- NA
melted_cormat_aster$value[which(melted_cormat_aster$value==1)] <- "*"
melted_cormat$aster <- melted_cormat_aster$value

# put boxes into melted_cormat
melted_cormat$boxes <- melted_cormat_boxes$value

# Set limits of values
limits <- c(min(melted_cormat$value),max(melted_cormat$value))

# Heatmap, where red means y group is greater than x group
ggheatmapArea <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(limits), 
                       space = "Lab", 
                       name="Pairwise difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label=aster), color = "black", shape=8) +
  geom_text(aes(Var2, Var1, label=boxes), color = "black", size=5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))






##### 
##### critical breadth:
##### 

dat <- read.csv("Analysis output/pairwise_tables_groups_av/breadthBT_pairwise_newthresh.csv", 
                row.names=1)[c(1,7,2,8,3,9,4,10,5,11,6,12),c(1,7,2,8,3,9,4,10,5,11,6,12)] 
colnames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
rownames(dat) <- c("N1 2010", "N1 2017", "N2 2010", "N2 2017",
                   "C1 2010", "C1 2017", "C2 2010", "C2 2017",
                   "S1 2010", "S1 2017", "S2 2010", "S2 2017")
dat <- dat[12:1,12:1]

# create a matrix the same size as the pairwise comparison table, but containing 
# 0's for non-significant comparisons and 1's for significant comparisons
aster <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(aster)[1]){
  vals <- dat[,i]
  for(j in 1:dim(aster)[1]){
    val <- as.numeric(as.character(vals[j]))
    if(is.na(val)==FALSE){aster[j,i] <- 0} else{aster[j,i] <- 1}
  }
}

# create a matrix the same size as the pairwise comparison table, but containing 
# numeric values (no asterisks)
nums <- matrix(NA, nrow=12, ncol=12)
for(i in 1:dim(nums)[1]){
  nums[,i] <- as.numeric(stringr::str_remove_all(as.character(dat[,i]), "[*]")) * (-1)
}
rownames(nums) <- rownames(dat)
colnames(nums) <- colnames(dat)



# plot heatmap with asterisks for significant comparisons
# datasets: aster (0=nonsignif/1=signif matrix), and nums (values of dat without asterisks)
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(nums)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Now print asterisks onto heatmap. First get the melted upper triangle of aster
upper_tri_aster <- get_upper_tri(aster)
melted_cormat_aster <- melt(upper_tri_aster, na.rm = TRUE)
melted_cormat_aster$value[which(melted_cormat_aster$value==0)] <- NA
melted_cormat_aster$value[which(melted_cormat_aster$value==1)] <- "*"
melted_cormat$aster <- melted_cormat_aster$value

# put boxes into melted_cormat
melted_cormat$boxes <- melted_cormat_boxes$value

# Set limits of values
limits <- c(min(melted_cormat$value),max(melted_cormat$value))

# Heatmap, where red means y group is greater than x group
ggheatmapBreadth <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(limits), 
                       space = "Lab", 
                       name="Pairwise difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label=aster), color = "black", shape=8) +
  geom_text(aes(Var2, Var1, label=boxes), color = "black", size=5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))











################ Plot everything together
# Heatmaps indicate pairwise difference between groups,
# where red means x group is greater than y group 
# and, blue means x group is less than y group
# comparisons with brackets indicate those that compare
# years within each population.




library(ggpubr)
theme_set(theme_minimal())
pdf("Figures/Figure S4.pdf", height=15, width=15)
figure <- ggarrange(ggheatmapTopt, ggheatmapB50, ggheatmapB80, 
                    ggheatmapBreadth, ggheatmapXmin, ggheatmapXmax,
                    ggheatmapmax_RGR, ggheatmapArea, 
                    labels = c("A) Thermal optimum", "B) B50", "C) B80",
                               "D) Critical breadth","E) Lower thermal limit", "F) Upper thermal limit", 
                               "G) Performance maximum", "H) Area"),
                    label.x = c(0,0.3,0.3,0,0,0,0,0.3), label.y=1, 
                    ncol = 3, nrow = 3)
figure
dev.off()
 







