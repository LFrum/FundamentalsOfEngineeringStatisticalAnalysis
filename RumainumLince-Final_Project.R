# Lince Rumainum
# Fundamentals of Engineering Statistical Analysis
# Final Project

# Libraries for the Final Project
library(tidyverse)
library(tidyr)
library(dplyr)
library(outliers)
library(EnvStats)
library(ggplot2)
library(mice)
library(caret)
library(earth)
library(Metrics) 
# ggbiplot is attached later on 
# to avoid conflict of plyr and dplyr package
#library(ggbiplot)

# upload data
setwd("C:/Users/Lince/Documents/Fall 2019/DSA5013-FundEngrStatistics/FinalProject/")
dat <- read.csv("pmsm_temperature_data.csv")

# view the number of row and column
nrow(dat) # 998070
ncol(dat) # 13

# make a copy of original data containing unmodified data
all_dat <- dat

# look at missing data
myMeanNAsFun <- function(x) mean(is.na(x))
apply(dat, 2, myMeanNAsFun)
# missing data using mice package
md.pattern(dat) # no missing data

# view info about the data set
summary(dat)
str(dat)

# create 4x3 for display for plots
par(mfrow=c(4,3)) 
# create histogram for all attributes of original data
for (i in 1:(ncol(dat)-1)){
  curCol <-  toString(colnames(dat[i]))
  hist(dat[,i], xlab = curCol, main =paste(c("Histogram of ", curCol)))
}
par(mfrow=c(1,1)) # reset display

# create 4x3 for display for plots
par(mfrow=c(4,3)) 
# create box plot for all attributes of original data
for (i in 1:(ncol(dat)-1)){
  curCol <-  toString(colnames(dat[i]))
  boxplot(dat[i],  main =paste(c("Box Plot of ", curCol)))
}
par(mfrow=c(1,1)) # reset display

# since the every profile ID is an independent measurement session, 
# we are going to restrict it to data for profile_id = 4
# create data frame for profile_id = 4
dfProfile4 <- dat # copy data
dfProfile4 <- dfProfile4[which(dfProfile4$profile_id == 4),] 

# view the number of row and column
nrow(dfProfile4) # 33423
ncol(dfProfile4) # 13

# remove proflie_id since we know all the profile_id = 4
dfProfile4 <- dfProfile4[,-c(13)]

# view the histogram and boxplot of the data of profile_id = 4
# create 4x3 for display for plots
par(mfrow=c(4,3)) 
# create histogram for all attributes of profile ID 4
for (i in 1:ncol(dfProfile4)){
  curCol <-  toString(colnames(dfProfile4[i]))
  hist(dfProfile4[,i], xlab = curCol, main =paste(c("Histogram of ", curCol)))
}
par(mfrow=c(1,1)) # reset display

# create 4x3 for display for plots
par(mfrow=c(4,3)) 
# create box plot for all attributes of profile ID 4
for (i in 1:ncol(dfProfile4)){
  curCol <-  toString(colnames(dfProfile4[i]))
  boxplot(dfProfile4[i],  main =paste(c("Box Plot of ", curCol)))
}
par(mfrow=c(1,1)) # reset display

# since the data set was already normalized by the author,
# normalizing by taking the log transformation would actually create
# a lot of NAs NANs values so I proceed in finding the outliers instead


# using rotnerTest to detect outliers in the data
outliersPID4 <- c(0,0,0,0,0,0,0,0,0,0,0,0)
for(i in 1:ncol(dfProfile4)){
  rt1 <- rosnerTest(dfProfile4[,i], k = round(nrow(dfProfile4)*0.25,0), warn = F)
  outliersPID4[i] <- length(rt1$all.stats$Outlier [rt1$all.stats$Outlier == TRUE])
}

# show how many outliers for each columns
outliersPID4

#colnames(dfProfile4)
# Grubbs test for the single outlier in coolant
grubbs.test(dfProfile4[,2])
# plot the pm vs coolant to see the outlier
ggplot(dfProfile4, aes(x = coolant, y = pm)) + geom_point()
# update data frame without the outlier value
log_pid4 <- log_pid4 [-which(log_pid4$coolant < -1.1565),]

# get the column with the highest number of outliers
col_out <- which.max(outliersPID4)

# create an updated dfProfile4 for the column with the large outliers
log_pid4 <- dfProfile4
# perform a log transformation on this data
log_pid4$i_q<- log(log_pid4[,col_out]+1)
# see the update box plot of i_q
curCol <-  toString(colnames(log_pid4[col_out]))
boxplot(log_pid4[col_out],  main =paste(c("Box Plot of ", curCol)))
# summary of the new i_q values
summary(log_pid4[col_out])

# look at the updated data
ggplot(log_pid4, aes(x = coolant, y = pm)) + geom_point()

# check for outliers on i_q variable after updated
rt_iq <- rosnerTest(log_pid4[,col_out], k = round(nrow(log_pid4)*0.25,0), warn = F)
length(rt_iq$all.stats$Outlier [rt_iq$all.stats$Outlier == TRUE]) # 0 outliers

# run PCA on the data
dfProfile4_pca <- prcomp(log_pid4[,], scale = T)
summary(dfProfile4_pca) # summary of PCA

# 7 principal components needed to reach a ~99.59% cummulative proportion 
# 4 principal components needed to reach a ~93.63% cummulative proportion 
# 3 principal components needed to reach a ~86.37% cummulative proportion 
# 2 principal components needed to reach a ~75.40% cummulative proportion 
# 1 principal component needed to reach a ~54.95% cummulative proportion 

# create PCA plot
library(ggbiplot)
# biplot for PCA 1 and PCA 2
ggbiplot(dfProfile4_pca, obs.scale = 1, var.scale = 1, varname.size = 6, varname.adjust = 1.1, varname.abbrev = FALSE,
         labels.size = 10, alpha = 0.1, choices = (c(1,2)), circle = TRUE, ellipse = TRUE) +
  xlim(-5,5) +
  ylim(-5,5)
# to resolve error for dplyr
detach(package:ggbiplot)
detach(package:plyr)


# separate data into training and test data
df_Train1 <-  sample_frac(log_pid4, 0.70) # 70% of data for training data 
rowsNum <- as.numeric(rownames(log_pid4)) # turn rownames() to num
df_Test  <-  log_pid4[-rowNums,] # 30% of data for testing data

# turn pm of Test to NAs
dfTest_wNAs <- df_Test
dfTest_wNAs$pm <- NA

########################################################################################################################

# LINEAR MODEL

########################################################################################################################
# create linear model using all attributes as predictors
lm1 = lm(pm ~ ., data = df_Train1) 

# analysis the behavior of the linear model 1 
summary(lm1)
par(mfrow=c(2,2)) 
plot(lm1)
par(mfrow=c(1,1)) 
order(varImp(lm1))
# variable importance of the linear model
varImp(lm1)

plot(varimp_lm, main="Variable Importance with Linear Model")
dfTest_wNAs$pm <- NA
# predict the shares and replace anything < 0 with 0
result_shares = predict(lm1, newdata = dfTest_wNAs)  

# RMSE and MSE using the first linear model
rmse(df_Test$pm, result_shares) # 0.2398643
mean((df_Test$pm - result_shares) ^ 2) # 0.05753489


# from the result of the variance importane in first linear model,
# I reduce the dimensionality to the top 4 importance attributes
# as pm predictors to create a second linear model
lm2 = lm(pm ~ stator_winding + stator_tooth + u_q + motor_speed, data = df_Train1) 

# analysis the behavior of the linear model 2
summary(lm2)
par(mfrow=c(2,2)) 
plot(lm2)
par(mfrow=c(1,1)) 

# put NAs for pm
dfTest_wNAs$pm <- NA
# predict the pm
result_pm2 = predict(lm2, newdata = dfTest_wNAs) 

# RMSE and MSE using the second linear model
rmse(df_Test$pm, result_pm2) # 0.2279918
mean((df_Test$pm - result_pm2) ^ 2) # 0.05198024



########################################################################################################################

# NON- LINEAR MODEL

########################################################################################################################

# MARS model
marsFit <- earth(pm ~ ., df_Train1, degree = 3, nfold = 3)
# look at summary, the behavior of the model, and its variable importance
summary(marsFit)
plot(marsFit)
varImp(marsFit)

# predict result for MARS model
dfTest_wNAs$pm <- NA
result_mars <- predict(marsFit, dfTest_wNAs)

# RMSE and MSE for the MARS fit
rmse(df_Test$pm, result_mars) # 0.0448497
mean((df_Test$pm - result_mars) ^ 2) # 0.002011496

# using the same attributes as linear model 2
df_Train2 <- df_Train1[,c(4,5,9,11,12)]
marsFit2 <- earth(pm ~ ., df_Train2, degree = 3, nfold = 3)
# look at summary, the behavior of the model, and its variable importance
summary(marsFit2)
plot(marsFit2)
plot(varImp(marsFit2))

# predict result for second MARS model
dfTest_wNAs$pm <- NA
result_mars2 <- predict(marsFit2, dfTest_wNAs)

# RMSE and MSE for the second MARS fit
rmse(df_Test$pm, result_mars2) # 0.05933646
mean((df_Test$pm - result_mars2) ^ 2) # 0.003520815

# using the same attributes as linear model 2 actually 
# resulted in a higher RMSE and MSE values