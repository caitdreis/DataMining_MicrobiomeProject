#Microbiome Project

#----------------- Set-up #####

#----------------- Packages
library(plyr)
library(tidyverse) # Load the core tidyverse packages: ggplot2, tibble, 
# tidyr, readr, purrr, and dplyr
library(vegan) #calculation of diversity metrics
library(psych) #descriptive statistics
library(class) #package for KNN model
library(randomForest) #package for Random Forest model
library (gbm) #package for boosting model
library(caret)
library(pscl) #For logistic regression R^2
library(ROCR) #For ROC curves


#----------------- Working Directory
setwd("~/Documents/GitHub/DataMining_MicrobiomeProject")

#----------------- Read in Data
microbiome <- read.csv("MicrobiomeWithMetadata.csv", encoding = 'utf-8', stringsAsFactors = FALSE)
#View(microbiome)

#----------------- Grading Criteria #####
#Data exploration
#Data cleaning (missing values, outliers, etc.)
#Rationale for the selected statistical modeling methods
#Correct implementation and use of statistical modeling methods
#Appropriate model selection approach (train/test, cross-validation, etc.)
#Thoroughly documented code (R comments)
#Appropriate use of functions to reduce code complexity and redundancy
#Writing quality for final report, evaluated in terms of conformance to process outline, level of 
#detail, and correctness.

#----------------- Data Exploration & Cleaning ####
#This dataset was pre-curated from the original Science article

nrow(microbiome) #675 observations

sum(is.na(microbiome)) #No missing data requiring imputation

colnames(microbiome[1:5])
#"Diet" "Source" "Donor" "CollectionMet" "Sex"  
#Beyond this, all column names are OTU groups

table(microbiome$Donor) #Distribution of donor groups among the observations (see below for key)
#0   1   2   3   4   5   6 
#313  64   6 151  92  46   3 

#Donor	
# 0 - HMouseLFPP -- mice
# 1 - CONVR
# 2 - Human
# 3 - Fresh -- mice
# 4 - Frozen -- mice
# 5 - HMouseWestern -- mice
# 6 - CONVD

#We will be using only non-human (mouse) donors for our analysis...
microbiome <- subset(microbiome, Donor == 0 | Donor == 3 | Donor == 4 | Donor == 5)

nrow(microbiome) #602 observations left

#Further explore the distribution of variables in this set...
table(microbiome$Diet)
# 0   1   4 
# 350 243 9 

#Diet	
#0	LFPP: low fat, high-plant polysaccharide diet
#1	Western: high-fat, high-sugar Western diet
#4	Suckling

#Remove all rows with suckling diets (want to compare only high and low fat diets)...
microbiome <- subset(microbiome, Diet == 0 | Diet == 1)

nrow(microbiome) #593 -- FINAL NUMBER OF ROWS

#We will be using only Diet, Source, Sex, and OTU groups going forward
microbiome <- microbiome[,c(1,2,5:100)]

#Check for class types of each non-OTU column
class((microbiome[,1])) #integer
class((microbiome[,2])) #integer
class((microbiome[,3])) #integer

#All of these should be factors
microbiome[,1] <- as.factor(microbiome[,1])
microbiome[,2] <- as.factor(microbiome[,2])
microbiome[,3] <- as.factor(microbiome[,3])

#Double check for correct conversion
class((microbiome[,1])) #factor - good

#Explore the distribution of sex among the observations...
table(microbiome$Sex)
#  0   1 
# 543  50  -- Uneven distribution (see below for key)

#Sex
#0 - Male
#1 - Female

plot(Diet~Sex, data=microbiome)
#Females appear to be more unevenly split between diets (favoring the low fat diet)

#   Diet	
#0	LFPP: low fat, high-plant polysaccharide diet
#1	Western: high-fat, high-sugar Western diet

#Repeat for source (again, refer to key below)
table(microbiome$Source)
#  0   1   2   3   4   5   6   7   8   9  10  11  12 
# 13  13  15  10 444  10  12  12   9  13   9   7  26

plot(Diet~Source, data=microbiome)
#11 is interesting, as it shows that samples from the stomach are classified as high fat only
#This is likely a coincidence

#Sources
# 0 - Cecum1
# 1 - Cecum2
# 2 - Colon1
# 3 - Colon2
# 4 - Feces
# 5 - SI1
# 6 - SI13
# 7 - SI15
# 8 - SI2
# 9 - SI5
# 10 - SI9
# 11 - Stomach
# 12 - Cecum


#----------------- Feature Engineering of Diversity Metrics ####

#Please refer to our paper for an explanation of these new variables

#Shannon Diversity Index
microbiome$ShannonIndex <- NULL
microbiome$ShannonIndex <- diversity(microbiome[,4:98])
summary(microbiome$ShannonIndex)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000001 0.348900 0.693300 1.214000 1.168000 4.433000

#Reyyi Index
#specnumber fucntion finds the number of species in the sample
specnumber(microbiome[,4:98]) #95 species present in each sample
h <- diversity(microbiome[,4:98])
J <- h/log(specnumber(microbiome[,4:98])) #assessing Pielou's eveness

#now take a random subset of six sites to assess Renyi diversity index
set.seed(2)
k <- sample(nrow(microbiome[,4:98]), 6) 
R <- renyi(microbiome[,4:98][k,])
plot(R) #can visualize these six sites, label sites chosen

#also make new column for Renyi index 
microbiome$Renyi <- NULL
microbiome$Renyi <- renyi(microbiome[,4:98], scales = 32)
summary(microbiome$Renyi)
# Min.  1st Qu.   Median   Mean  3rd Qu.    Max. 
# 0.0000  0.1216  0.5273  0.9932  0.9235  4.0790  

#Is there a significant difference between these indexes?
comp.index <- aov(microbiome$ShannonIndex ~ microbiome$Renyi, data=microbiome)
plot(comp.index)
summary(comp.index)
#                 Df Sum Sq Mean Sq F value Pr(>F)    
#microbiome$Renyi   1 1183.2    1183   37225 <2e-16 ***
#Residuals        591   18.8       0  
#the signifcant p-value tells us that there is a significant difference between the mean
#Shannon Index and Renyi index.

#----------------- Further Descriptive and ANOVA Exploration with Featured Variables 
#----------------- Shannon Index
describeBy(microbiome$ShannonIndex, microbiome$Diet)
#Descriptive statistics by group 
#group: 0
#   vars   n mean   sd median trimmed  mad min  max range skew kurtosis   se
# X1    1 350 1.26 1.47   0.69    1.03 0.62   0 4.43  4.43 1.44     0.53 0.08

#   group: 1
# vars   n mean   sd median trimmed  mad min  max range skew kurtosis   se
# X1    1 243 1.15 1.36   0.71    0.89 0.75   0 4.42  4.42 1.65     1.44 0.09

#significance testing with ANOVA using diet
fit <- aov(microbiome$ShannonIndex ~ microbiome$Diet, data=microbiome)
plot(fit)
summary(fit)
#               Df Sum Sq Mean Sq F value Pr(>F)
# microbiome$Diet   1    1.7   1.720   0.847  0.358
# Residuals       591 1200.2   2.031  
#No significance of diet on the Shannon Index

#posthoc TukeyHSD (Honestly Significant Differences)
TukeyHSD(fit) #where fit comes from aov()
#Tukey multiple comparisons of means
#95% family-wise confidence level
#Fit: aov(formula = microbiome$ShannonIndex ~ microbiome$Diet, data = microbiome)
#$`microbiome$Diet`
#       diff        lwr       upr     p adj
#1-0 -0.1095192 -0.3432244 0.1241861 0.3577576
#Tukey's posthoc testing tells us that there is truly no significance in diet's
#effect on the mean diversity index measure.

#significance testing with ANOVA using sex
fit.sex <- aov(microbiome$ShannonIndex ~ microbiome$Sex, data=microbiome)
plot(fit.sex)
summary(fit.sex)
#               Df Sum Sq Mean Sq F value Pr(>F)
#microbiome$Sex   1      0  0.0013   0.001   0.98
#Residuals      591   1202  2.0338 

#posthoc TukeyHSD 
TukeyHSD(fit.sex) 
#         diff        lwr       upr     p adj
#1-0 0.005360838 -0.4085735 0.4192952 0.9797162
#Tukey's posthoc testing tells us that there is no significance in sex's
#effect on the mean diversity index measure.

#----------------- Renyi Index
describeBy(microbiome$Renyi, microbiome$Diet)
#Descriptive statistics by group 
#group: 0
# vars   n mean   sd median trimmed  mad min  max range skew kurtosis   se
# X1    1 350 1.04 1.37   0.56     0.8 0.56   0 4.07  4.07 1.53     0.66 0.07

#group: 1
# vars   n mean   sd median trimmed  mad min  max range skew kurtosis   se
# X1    1 243 0.93 1.26   0.53    0.67 0.62   0 4.08  4.08 1.76     1.64 0.08

#significance testing with ANOVA using diet
fit.ren <- aov(microbiome$Renyi ~ microbiome$Diet, data=microbiome)
plot(fit.ren)
summary(fit.ren)
#               Df Sum Sq Mean Sq F value Pr(>F)
#microbiome$Diet   1    1.8   1.765   1.007  0.316
# Residuals       591 1036.4   1.754    
#No significant on the effect of diet on the Renyi diversity index

#posthoc TukeyHSD (Honestly Significant Differences)
TukeyHSD(fit.ren) #where fit comes from aov()
#Tukey multiple comparisons of means
#95% family-wise confidence level
#Fit: aov(formula = microbiome$Renyi ~ microbiome$Diet, data = microbiome)
#$`microbiome$Diet`
#       diff        lwr       upr     p adj
#1-0 -0.1109416 -0.3281099 0.1062267 0.3161206
#Tukey's posthoc testing tells us that there is truly no significance in diet's
#effect on the mean Renyi diversity index measure.

#significance testing with ANOVA using sex
fit.sex.ren <- aov(microbiome$Renyi ~ microbiome$Sex, data=microbiome)
plot(fit.sex.ren)
summary(fit.sex.ren)
#               Df Sum Sq Mean Sq F value Pr(>F)
#microbiome$Sex   1    0.1  0.0685   0.039  0.844
# Residuals      591 1038.1  1.7565 

#posthoc TukeyHSD 
TukeyHSD(fit.sex.ren) 
#         diff        lwr       upr     p adj
#1-0 0.038679 -0.3460048 0.4233628 0.8435246


#----------------- Training and Testing Sets ######

#Cross validation of testing data with testing data entered in for test
set.seed(17)
sample = sample(1:nrow(microbiome), nrow(microbiome) * .75)
test <- microbiome[-sample, c(1:98)]
train <- microbiome[sample, c(1:98)]
test_div <- microbiome[-sample,] #includes both diversity indices
train_div <- microbiome[sample,] #includes both diversity indices

microbiome_div <- microbiome
microbiome <- microbiome[, c(1:98)]


#----------------- Youden's Index ######

#This will be used to establish a cutoff/threshold for the logistic regression model predictions

YoudensIndex <- function(X){   # Input to the function is the model
  # Split validation set into MMR gotten and not gotten observations
  Test_LowFat <- test[which(test$Diet==0),]  
  Test_HighFat <- test[which(test$Diet==1),] 
  
  # Initialize the data frame to store results in
  Results <- as_data_frame(matrix(NA, ncol = 11))
  names(Results) <- c("CutoffProb","TP","TN","FP","FN","Sensitivity",
                      "Specificity","PPV","NPV","Accuracy","Youden's Index")
  
  # Initialize variables to store Youden's Index
  YI <- 0
  SSsum <- 0
  
  # Loop over various cutoff probabilities from 0.01 to 0.99
  sequenceprobs <- seq(0.01,0.99,0.01)
  for (i in sequenceprobs){
    # Get the predictions from the model for No MMR group
    probs <- as.vector(predict(X,newdata=Test_HighFat, type="response"))
    
    # Use cutoff to classify vaccination status
    preds <- rep(0,nrow(Test_HighFat))  # Initialize prediction vector
    preds[probs>i] <- 1
    
    # Add predictions to mortality group data frame
    Positives <- cbind(Test_HighFat,preds)
    
    # Identify true positives and false negatives
    TruePositives <- nrow(Positives[which(Positives$preds == 1),])
    FalseNegatives <- nrow(Positives[which(Positives$preds == 0),])
    
    # Get the predictions from the model for MMR Received group
    probs <- as.vector(predict(X,newdata=Test_LowFat, type="response"))
    
    # Use cutoff to classify mortality
    preds <- rep(0,nrow(Test_LowFat))  # Initialize prediction vector
    preds[probs>i] <- 1 
    
    # Add predictions to MMR received group data frame
    Negatives <- cbind(Test_LowFat,preds)
    
    # Identify false positives and true negatives
    FalsePositives <- nrow(Negatives[which(Negatives$preds == 1),])
    TrueNegatives <- nrow(Negatives[which(Negatives$preds == 0),])
    
    # Generate statistics
    Sensitivity <- TruePositives / (TruePositives + FalseNegatives)
    Specificity <- TrueNegatives / (TrueNegatives + FalsePositives)
    PPV <- TruePositives / (TruePositives + FalsePositives)
    NPV <- TrueNegatives / (TrueNegatives + FalseNegatives)
    Accuracy <- (TruePositives + TrueNegatives) / (TruePositives + TrueNegatives + FalsePositives + FalseNegatives)
    
    # Generate entry to append to the results dataframe
    Entry <- as_data_frame(matrix(NA, nrow = 1, ncol = 11))
    names(Entry) <- c("CutoffProb","TP","TN","FP","FN","Sensitivity",
                      "Specificity","PPV","NPV","Accuracy","Youden's Index")
    Entry[1,] <- c(i,TruePositives,TrueNegatives,FalsePositives,FalseNegatives,
                   Sensitivity,Specificity,PPV,NPV,Accuracy,0)
    
    # Append to the dataframe
    Results <- rbind(Results,Entry)
    
    # Update Youden's Index variables
    candidate <- Sensitivity + Specificity - 1
    if (candidate > SSsum){
      YI <- i
      SSsum <- candidate
    }
  }
  
  # Retrieve Youden's Index and put into dataframe
  Results$`Youden's Index` <- YI
  
  # Drop the first row (which is all NAs)
  Results <- Results[-1,]
  
  # Output the resulting dataframe
  return(Results)
}


#----------------- Youden's Index for Testing Set with Diversity Indices ######

YoudensIndexDiv <- function(X){   # Input to the function is the model
  # Split validation set into MMR gotten and not gotten observations
  Test_LowFat <- test_div[which(test_div$Diet==0),]  
  Test_HighFat <- test_div[which(test_div$Diet==1),] 
  
  # Initialize the data frame to store results in
  Results <- as_data_frame(matrix(NA, ncol = 11))
  names(Results) <- c("CutoffProb","TP","TN","FP","FN","Sensitivity",
                      "Specificity","PPV","NPV","Accuracy","Youden's Index")
  
  # Initialize variables to store Youden's Index
  YI <- 0
  SSsum <- 0
  
  # Loop over various cutoff probabilities from 0.01 to 0.99
  sequenceprobs <- seq(0.01,0.99,0.01)
  for (i in sequenceprobs){
    # Get the predictions from the model for No MMR group
    probs <- as.vector(predict(X,newdata=Test_HighFat, type="response"))
    
    # Use cutoff to classify vaccination status
    preds <- rep(0,nrow(Test_HighFat))  # Initialize prediction vector
    preds[probs>i] <- 1
    
    # Add predictions to mortality group data frame
    Positives <- cbind(Test_HighFat,preds)
    
    # Identify true positives and false negatives
    TruePositives <- nrow(Positives[which(Positives$preds == 1),])
    FalseNegatives <- nrow(Positives[which(Positives$preds == 0),])
    
    # Get the predictions from the model for MMR Received group
    probs <- as.vector(predict(X,newdata=Test_LowFat, type="response"))
    
    # Use cutoff to classify mortality
    preds <- rep(0,nrow(Test_LowFat))  # Initialize prediction vector
    preds[probs>i] <- 1 
    
    # Add predictions to MMR received group data frame
    Negatives <- cbind(Test_LowFat,preds)
    
    # Identify false positives and true negatives
    FalsePositives <- nrow(Negatives[which(Negatives$preds == 1),])
    TrueNegatives <- nrow(Negatives[which(Negatives$preds == 0),])
    
    # Generate statistics
    Sensitivity <- TruePositives / (TruePositives + FalseNegatives)
    Specificity <- TrueNegatives / (TrueNegatives + FalsePositives)
    PPV <- TruePositives / (TruePositives + FalsePositives)
    NPV <- TrueNegatives / (TrueNegatives + FalseNegatives)
    Accuracy <- (TruePositives + TrueNegatives) / (TruePositives + TrueNegatives + FalsePositives + FalseNegatives)
    
    # Generate entry to append to the results dataframe
    Entry <- as_data_frame(matrix(NA, nrow = 1, ncol = 11))
    names(Entry) <- c("CutoffProb","TP","TN","FP","FN","Sensitivity",
                      "Specificity","PPV","NPV","Accuracy","Youden's Index")
    Entry[1,] <- c(i,TruePositives,TrueNegatives,FalsePositives,FalseNegatives,
                   Sensitivity,Specificity,PPV,NPV,Accuracy,0)
    
    # Append to the dataframe
    Results <- rbind(Results,Entry)
    
    # Update Youden's Index variables
    candidate <- Sensitivity + Specificity - 1
    if (candidate > SSsum){
      YI <- i
      SSsum <- candidate
    }
  }
  
  # Retrieve Youden's Index and put into dataframe
  Results$`Youden's Index` <- YI
  
  # Drop the first row (which is all NAs)
  Results <- Results[-1,]
  
  # Output the resulting dataframe
  return(Results)
}


#----------------- Basic Binomial Logistic Regression ######

#NOTE: the first few logistic models were run WITHOUT the diversity indices

#Create a model with all possible predictors
model <- glm(Diet~.,family=binomial(link='logit'),data=train)
summary(model) #Every variable is significant...

#ANOVA 
anova(model, test="Chisq")
#A large p-value here indicates that the model without the variable explains 
#more or less the same amount of variation.

# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                     443      600.3              
# Source 12     59.3       431      541.0 2.998e-08 ***
# Sex     1      3.6       430      537.4 0.0591044 .  
# OTU0    1      2.2       429      535.2 0.1371864    
# OTU1    1      2.8       428      532.4 0.0926576 .  
# OTU2    1     37.5       427      494.9 9.225e-10 ***
# OTU3    1      6.5       426      488.3 0.0105465 *  
# OTU4    1      0.5       425      487.8 0.4627754    
# OTU5    1      0.6       424      487.2 0.4295599    
# OTU6    1      0.0       423      487.1 0.8667479    
# OTU7    1      0.5       422      486.6 0.4750674    
# OTU8    1      8.3       421      478.4 0.0040625 ** 
# OTU9    1     33.9       420      444.4 5.658e-09 ***
# OTU10   1      8.4       419      436.0 0.0037840 ** 
# OTU11   1      0.4       418      435.6 0.5108424    
# OTU12   1      0.0       417      435.6 0.8449512    
# OTU13   1      4.2       416      431.4 0.0406067 *  
# OTU14   1      0.6       415      430.8 0.4448434    
# OTU15   1      1.8       414      429.0 0.1752748    
# OTU16   1      1.1       413      427.9 0.2912993    
# OTU17   1      1.1       412      426.8 0.3026798    
# OTU18   1      9.5       411      417.3 0.0020388 ** 
# OTU19   1      0.6       410      416.7 0.4489548    
# OTU20   1      0.4       409      416.3 0.5055763    
# OTU21   1      0.5       408      415.7 0.4636021    
# OTU22   1      2.9       407      412.8 0.0884765 .  
# OTU23   1      0.8       406      412.0 0.3804481    
# OTU24   1      2.5       405      409.5 0.1124038    
# OTU25   1     15.8       404      393.8 7.213e-05 ***
# OTU26   1      2.9       403      390.9 0.0913290 .  
# OTU27   1      0.0       402      390.9 0.8350321    
# OTU28   1      0.0       401      390.9 0.9077022    
# OTU29   1      0.0       400      390.9 0.9689933    
# OTU30   1      0.6       399      390.3 0.4524332    
# OTU31   1      4.5       398      385.8 0.0337931 *  
# OTU32   1      0.6       397      385.2 0.4387198    
# OTU33   1      4.3       396      380.9 0.0385292 *  
# OTU34   1     11.8       395      369.1 0.0005996 ***
# OTU35   1      0.0       394     7497.1 1.0000000    
# OTU36   1   7136.6       393      360.5 < 2.2e-16 ***
# OTU37   1      0.9       392      359.6 0.3411374    
# OTU38   1      0.3       391      359.3 0.5668152    
# OTU39   1      3.8       390      355.4 0.0507401 .  
# OTU40   1      4.0       389      351.5 0.0457329 *  
# OTU41   1     50.9       388      300.6 9.756e-13 ***
# OTU42   1      2.8       387      297.8 0.0957792 .  
# OTU43   1      0.3       386      297.5 0.5623545    
# OTU44   1      0.2       385      297.3 0.6938438    
# OTU45   1      0.4       384      296.9 0.5262524    
# OTU46   1      0.0       383      296.9 0.9957600    
# OTU47   1      1.9       382      294.9 0.1629259    
# OTU48   1      2.4       381      292.6 0.1221490    
# OTU49   1      0.4       380      292.2 0.5274164    
# OTU50   1      0.0       379     6415.8 1.0000000    
# OTU51   1   6125.3       378      290.5 < 2.2e-16 ***
# OTU52   1      0.0       377      290.5 0.9988007    
# OTU53   1      0.5       376      290.0 0.4682265    
# OTU54   1      0.0       375     5767.0 1.0000000    
# OTU55   1      0.0       374     6055.3 1.0000000    
# OTU56   1      0.0       373     6199.5 1.0000000    
# OTU57   1   1081.3       372     5118.2 < 2.2e-16 ***
# OTU58   1      0.0       371     5694.9 1.0000000    
# OTU59   1      0.0       370     6848.3 1.0000000    
# OTU60   1    937.1       369     5911.2 < 2.2e-16 ***
# OTU61   1      0.0       368     6271.6 1.0000000    
# OTU62   1    648.8       367     5622.8 < 2.2e-16 ***
# OTU63   1      0.0       366     5983.2 1.0000000    
# OTU64   1    288.3       365     5694.9 < 2.2e-16 ***
# OTU65   1    288.3       364     5406.5 < 2.2e-16 ***
# OTU66   1      0.0       363     5911.2 1.0000000    
# OTU67   1      0.0       362     5983.2 1.0000000    
# OTU68   1    144.2       361     5839.1 < 2.2e-16 ***
# OTU69   1   5594.5       360      244.6 < 2.2e-16 ***
# OTU70   1      3.0       359      241.6 0.0823970 .  
# OTU71   1     11.2       358      230.4 0.0008157 ***
# OTU72   1      0.0       357      230.4 0.9259223    
# OTU73   1      0.0       356     7136.6 1.0000000    
# OTU74   1   6921.8       355      214.9 < 2.2e-16 ***
# OTU75   1      0.0       354     5983.2 1.0000000    
# OTU76   1      0.0       353     6704.1 1.0000000    
# OTU77   1      0.0       352     7136.6 1.0000000    
# OTU78   1   1513.8       351     5622.8 < 2.2e-16 ***
# OTU79   1      0.0       350     5839.1 1.0000000    
# OTU80   1      0.0       349     6055.3 1.0000000    
# OTU81   1      0.0       348     6127.4 1.0000000    
# OTU82   1    793.0       347     5334.5 < 2.2e-16 ***
# OTU83   1      0.0       346     5334.5 1.0000000    
# OTU84   1      0.0       345     5334.5 1.0000000    
# OTU85   1   5168.6       344      165.8 < 2.2e-16 ***
# OTU86   1      0.0       343     4397.3 1.0000000    
# OTU87   1      0.0       342     5190.3 1.0000000    
# OTU88   1   1369.7       341     3820.6 < 2.2e-16 ***
# OTU89   1    360.4       340     3460.2 < 2.2e-16 ***
# OTU90   1      0.0       339     4109.0 1.0000000    
# OTU91   1    360.4       338     3748.5 < 2.2e-16 ***
# OTU92   1    504.6       337     3243.9 < 2.2e-16 ***
# OTU93   1      0.0       336     4181.1 1.0000000    
# OTU94   1    288.3       335     3892.7 < 2.2e-16 ***

#Many of the OTU groups do not help lower the deviance

pR2(model)
#      llh      llhNull           G2        McFadden      r2ML         r2CU 
# -1946.357283  -300.141633 -3292.431301    -5.484796  -1660.350015 -2239.857196

fitted.results <- predict(model,newdata=subset(test,select=c(2:98)),type='response')

LogisticFunction <- YoudensIndex(model)
unique(LogisticFunction$`Youden's Index`) #0.01

fitted.results <- ifelse(fitted.results > 0.01,1,0)

misClasificError <- mean(fitted.results != test$Diet)
print(paste('Accuracy',1-misClasificError))
#Accuracy 0.751677852348993

p <- predict(model,newdata=subset(test,select=c(2:98)),type='response')
p
pr <- prediction(p, test$Diet)
prf1 <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf1)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7410085

#----------------- Stepwise Logistic Regression ######

model.null = glm(Diet ~ 1, data=train, family = binomial(link="logit"))

step(model.null,
     scope = list(upper=model),
     direction="both",
     test="Chisq",
     data=train)

#Call:  glm(formula = Diet ~ OTU41 + Source + OTU77 + OTU9 + OTU54 + 
# OTU71 + OTU88 + OTU25 + OTU18 + OTU74 + OTU1 + OTU10 + OTU86 + 
#   OTU85 + OTU59 + OTU35 + OTU89 + OTU58 + OTU44 + OTU51 + OTU22 + 
#   OTU80 + OTU73 + OTU87 + OTU7 + OTU26 + OTU31 + OTU75 + OTU4 + 
#   OTU42 + OTU83 + OTU0 + OTU8 + OTU40 + OTU70 + OTU30 + OTU39 + 
#   OTU43 + OTU24, family = binomial(link = "logit"), data = train)
# 
# Degrees of Freedom: 443 Total (i.e. Null);  393 Residual
# Null Deviance:	    600.3 
# Residual Deviance: 170.5 	AIC: 272.5

#Step model with 39 predictors (fewer than half of the full model's predictors)
step.model <- glm(formula = Diet ~ OTU41 + Source + OTU77 + OTU9 + OTU54 + 
                    OTU71 + OTU88 + OTU25 + OTU18 + OTU74 + OTU1 + OTU10 + OTU86 +
                    OTU85 + OTU59 + OTU35 + OTU89 + OTU58 + OTU44 + OTU51 + OTU22 +
                    OTU80 + OTU73 + OTU87 + OTU7 + OTU26 + OTU31 + OTU75 + OTU4 +
                    OTU42 + OTU83 + OTU0 + OTU8 + OTU40 + OTU70 + OTU30 + OTU39 +
                    OTU43 + OTU24, family = binomial(link = "logit"), data = train)

summary(step.model) 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -5.378e-01  2.525e+00  -0.213 0.831297    
# OTU41        1.209e+04  2.745e+03   4.404 1.06e-05 ***
# Source1     -1.852e+00  1.953e+00  -0.948 0.343167    
# Source2      3.956e-02  1.734e+00   0.023 0.981795    
# Source3     -3.047e+00  1.968e+00  -1.548 0.121583    
# Source4     -8.161e+00  1.646e+00  -4.958 7.14e-07 ***
# Source5     -5.117e+00  1.931e+00  -2.650 0.008049 ** 
# Source6     -1.654e+00  1.808e+00  -0.915 0.360086    
# Source7     -1.199e+00  2.099e+00  -0.571 0.567870    
# Source8     -3.932e+00  1.835e+00  -2.143 0.032092 *  
# Source9     -3.031e+00  1.740e+00  -1.742 0.081569 .  
# Source10    -5.685e-01  2.094e+00  -0.271 0.786010    
# Source11     2.244e+01  1.072e+05   0.000 0.999833    
# Source12    -1.172e+01  2.961e+00  -3.958 7.55e-05 ***
# OTU77        1.251e+03  3.376e+02   3.707 0.000210 ***
# OTU9        -3.743e+03  7.977e+02  -4.692 2.70e-06 ***
# OTU54        1.053e+10  8.131e+09   1.295 0.195317    
# OTU71        5.974e+04  6.837e+06   0.009 0.993029     -- note high p-value
# OTU88        1.201e+04  2.835e+03   4.238 2.25e-05 ***
# OTU25        4.876e+09  6.775e+09   0.720 0.471718    
# OTU18        5.134e+10  1.018e+10   5.043 4.59e-07 ***
# OTU74       -2.008e+10  8.236e+09  -2.438 0.014749 *  
# OTU1        -3.064e+10  9.198e+09  -3.332 0.000863 ***
# OTU10       -3.911e+04  1.047e+07  -0.004 0.997018     -- note high p-value    
# OTU86        6.309e+08  7.754e+09   0.081 0.935153     -- note high p-value    
# OTU85        6.257e+02  1.835e+02   3.409 0.000651 ***
# OTU59        5.432e+03  2.140e+03   2.539 0.011126 *  
# OTU35       -4.905e+04  5.146e+06  -0.010 0.992395     -- note high p-value    
# OTU89        1.798e+10  8.087e+09   2.223 0.026220 *  
# OTU58       -3.050e+09  9.560e+08  -3.190 0.001423 ** 
# OTU44        2.745e+10  8.604e+09   3.190 0.001423 ** 
# OTU51        1.125e+10  8.008e+09   1.405 0.160108    
# OTU22        1.885e+10  8.483e+09   2.222 0.026303 *  
# OTU80       -1.723e+10  8.292e+09  -2.078 0.037695 *  
# OTU73        2.671e+10  8.565e+09   3.119 0.001816 ** 
# OTU87        2.929e+09  8.316e+09   0.352 0.724667    
# OTU7        -2.197e+10  8.160e+09  -2.693 0.007081 ** 
# OTU26       -3.203e+10  9.625e+09  -3.328 0.000876 ***
# OTU31        1.081e+04  4.654e+03   2.323 0.020201 *  
# OTU75       -1.920e+10  7.862e+09  -2.442 0.014588 *  
# OTU4         2.147e+10  7.719e+09   2.781 0.005417 ** 
# OTU42        5.417e+03  2.578e+03   2.101 0.035601 *  
# OTU83       -1.373e+04  6.638e+03  -2.069 0.038563 *  
# OTU0        -1.151e+10  8.096e+09  -1.422 0.155000    
# OTU8         2.610e+03  1.085e+03   2.406 0.016143 *  
# OTU40        1.549e+10  8.002e+09   1.936 0.052841 .  
# OTU70        1.705e+10  8.031e+09   2.123 0.033786 *  
# OTU30        1.449e+10  7.912e+09   1.832 0.066954 .  
# OTU39        1.613e+10  8.523e+09   1.892 0.058495 .  
# OTU43        2.304e+03  1.293e+03   1.782 0.074741 .  
# OTU24       -1.221e+10  7.732e+09  -1.580 0.114185    

#ANOVA 
anova(step.model, test="Chisq")

#       Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                     443     600.28              
# OTU41   1   98.712       442     501.57 < 2.2e-16 ***
# Source 12   81.245       430     420.33 2.388e-12 ***
# OTU77   1   29.529       429     390.80 5.510e-08 ***
# OTU9    1   22.053       428     368.75 2.653e-06 ***
# OTU54   1   18.562       427     350.18 1.645e-05 ***
# OTU71   1   15.394       426     334.79 8.725e-05 ***
# OTU88   1   11.435       425     323.35 0.0007208 ***
# OTU25   1   10.817       424     312.54 0.0010055 ** 
# OTU18   1   12.477       423     300.06 0.0004119 ***
# OTU74   1    8.713       422     291.35 0.0031603 ** 
# OTU1    1    6.807       421     284.54 0.0090789 ** 
# OTU10   1    7.369       420     277.17 0.0066359 ** 
# OTU86   1    5.843       419     271.33 0.0156420 *  
# OTU85   1    5.156       418     266.17 0.0231625 *  
# OTU59   1    5.730       417     260.44 0.0166820 *  
# OTU35   1    4.469       416     255.97 0.0345116 *  
# OTU89   1    4.987       415     250.99 0.0255455 *  
# OTU58   1    3.642       414     247.34 0.0563481 .  
# OTU44   1    2.949       413     244.40 0.0859259 .  
# OTU51   1    4.597       412     239.80 0.0320361 *  
# OTU22   1    3.568       411     236.23 0.0588862 .  
# OTU80   1    5.263       410     230.97 0.0217838 *  
# OTU73   1    4.462       409     226.51 0.0346651 *  
# OTU87   1    4.251       408     222.25 0.0392235 *  
# OTU7    1    3.272       407     218.98 0.0704789 .  
# OTU26   1    6.924       406     212.06 0.0085040 ** 
# OTU31   1    4.097       405     207.96 0.0429509 *  
# OTU75   1    3.551       404     204.41 0.0595241 .  
# OTU4    1    3.689       403     200.72 0.0547698 .  
# OTU42   1    2.990       402     197.73 0.0837558 .  
# OTU83   1    3.279       401     194.45 0.0701671 .  
# OTU0    1    3.209       400     191.24 0.0732154 .  
# OTU8    1    3.196       399     188.05 0.0738122 .  
# OTU40   1    3.582       398     184.46 0.0584185 .  
# OTU70   1    3.470       397     181.00 0.0625068 .  
# OTU30   1    1.974       396     179.02 0.1600272    
# OTU39   1    3.656       395     175.37 0.0558664 .  
# OTU43   1    2.275       394     173.09 0.1314678    
# OTU24   1    2.593       393     170.50 0.1073404    

pR2(step.model)
#      llh      llhNull           G2        McFadden      r2ML         r2CU 
# -85.2486643 -300.1416326  429.7859366    0.7159719    0.6201528    0.836603

fitted.results <- predict(step.model,newdata=subset(test,select=c(2:98)),type='response')

LogisticFunction <- YoudensIndex(step.model) #Use Youden's Index to select cutoff
unique(LogisticFunction$`Youden's Index`) #0.83

fitted.results <- ifelse(fitted.results > 0.83,1,0)

misClasificError <- mean(fitted.results != test$Diet)
print(paste('Accuracy',1-misClasificError))
#Accuracy 0.785234899328859 - Better, compared to the model with all predictors

p <- predict(step.model,newdata=subset(test,select=c(2:98)),type='response')
p
pr <- prediction(p, test$Diet)
prf2 <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf2)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7875417 -- Better!

#----------------- Manually Edited Stepwise Model #####

#Manually Selected Model (based off of summary from stepwise model) with 26 predictors
#For this, I only selected variables with a p-value relatively close to 0.1 or below
#I deleted each of the insignicant variables one by one to make sure each one's 
#removal didn't significantly affect the p-values of the other variables being
#considered for removal

step.model.2 <- glm(formula = Diet ~ OTU41 + Source + OTU77 + OTU9 + OTU54 + 
                    OTU88 + OTU18 + OTU74 + OTU1 + OTU70 +
                    OTU85 + OTU59 + OTU89 + OTU58 + OTU44 +
                    OTU80 + OTU73 + OTU7 + OTU26 + OTU31 +  
                    OTU75 + OTU4 + OTU42 + OTU8 + OTU40 + OTU39, 
                  family = binomial(link = "logit"), data = train)

summary(step.model.2)

pR2(step.model.2)
#McFadden
#0.5641426

fitted.results <- predict(step.model.2,newdata=subset(test,select=c(2:98)),type='response')

LogisticFunction <- YoudensIndex(step.model.2) #Use Youden's Index to select cutoff
unique(LogisticFunction$`Youden's Index`) #0.46

fitted.results <- ifelse(fitted.results > 0.46,1,0)

misClasificError <- mean(fitted.results != test$Diet)
print(paste('Accuracy',1-misClasificError))
#Accuracy 0.785234899328859 - Unchanged

p <- predict(step.model.2,newdata=subset(test,select=c(2:98)),type='response')
p
pr <- prediction(p, test$Diet)
prf3 <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf3)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.8201706 -- Better than automatic/original stepwise model

#----------------- Comparing ROC Curves

par(mfrow=c(1,3))
plot(prf1)
plot(prf2)
plot(prf3)

#Reset
par(mfrow=c(1,1))

#----------------- Basic Logistic Regression + Diversity Indices ######

#Create a model with all possible predictors
model <- glm(Diet~.,family=binomial(link='logit'),data=train_div)
summary(model) #Every variable is significant...

#ANOVA 
anova(model, test_div="Chisq") #See output -- neither diversity index appears to lower residual deviance

pR2(model)
# McFadden
#-5.604885 (versus -5.484796 prior)  

fitted.results <- predict(model,newdata=subset(test_div,select=c(2:100)),type='response')

LogisticFunction <- YoudensIndexDiv(model)
unique(LogisticFunction$`Youden's Index`) #0.01

fitted.results <- ifelse(fitted.results > 0.01,1,0)

misClasificError <- mean(fitted.results != test_div$Diet)
print(paste('Accuracy',1-misClasificError))
#Accuracy 0.751677852348993 - unchanged

p <- predict(model,newdata=subset(test_div,select=c(2:100)),type='response')
p
pr <- prediction(p, test_div$Diet)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7456433 - basically unchanged


#----------------- Stepwise Model + Diversity Indices ######

model <- glm(Diet~.,family=binomial(link='logit'),data=train_div)
model.null = glm(Diet ~ 1, data=train_div, family = binomial(link="logit"))

step(model.null,
     scope = list(upper=model),
     direction="both",
     test="Chisq",
     data=train)

#This stepwise selection (for the training set with the diversity indices) does not select the
#diversity variables. To be sure that these variables do not add anything to the model,
#we ran the same stepwise model, but manually added the diversity indices...

#Same step model as before, but with the addition of ShannonIndex and Renyi predictor variables
step.model.div <- glm(formula = Diet ~ OTU41 + Source + OTU77 + OTU9 + OTU54 + 
                    OTU71 + OTU88 + OTU25 + OTU18 + OTU74 + OTU1 + OTU10 + OTU86 +
                    OTU85 + OTU59 + OTU35 + OTU89 + OTU58 + OTU44 + OTU51 + OTU22 +
                    OTU80 + OTU73 + OTU87 + OTU7 + OTU26 + OTU31 + OTU75 + OTU4 +
                    OTU42 + OTU83 + OTU0 + OTU8 + OTU40 + OTU70 + OTU30 + OTU39 +
                    OTU43 + OTU24+ShannonIndex+Renyi, family = binomial(link = "logit"), data = train_div)

summary(step.model.div) 
#Neither diversity index appears to be significant in this model

pR2(step.model.div)
#McFadden
#0.7164411 (slightly higher than before)

#Previous McFadden
#0.7159719

fitted.results <- predict(step.model.div,newdata=subset(test_div,select=c(2:100)),type='response')

LogisticFunction <- YoudensIndexDiv(step.model)
unique(LogisticFunction$`Youden's Index`) #0.83

fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError <- mean(fitted.results != test_div$Diet)
print(paste('Accuracy',1-misClasificError))
#Accuracy 0.74496644295302 - Unchanged from before

p <- predict(step.model.div,newdata=subset(test_div,select=c(2:100)),type='response')
p
pr <- prediction(p, test_div$Diet)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7886541 -- Within the range of error, so it is likely not worth including the extra predictors
#Previous AUC of 0.7875417 


#----------------- K-Nearest Neighbors (KNN) ######

#Remove Diet column for KNN testing and training sets
training <- train[,!colnames(train) %in% "Diet"]
testing <- test[,!colnames(test) %in% "Diet"]

set.seed(10)
#create model with training data, test data, and training set classifier
knn.pred1 <- knn(training, testing, train$Diet)
knn.pred5 <- knn(training, testing, train$Diet, k=5)
knn.pred10 <- knn(training, testing, train$Diet, k=10)

#create a prediction table to understand how classes of sentiment values are being defined
pred.mat1 <- table("Predictions" = knn.pred1, Actual = test$Diet); pred.mat1
#             Actual
# Predictions  0  1
#           0 74 18
#           1 13 44

pred.mat5 <- table("Predictions" = knn.pred5, Actual = test$Diet); pred.mat5
#             Actual
# Predictions  0  1
#           0 77 21
#           1 10 41

pred.mat10 <- table("Predictions" = knn.pred10, Actual = test$Diet); pred.mat10
#              Actual
# Predictions  0  1
#           0 79 26
#           1 8  36

#assess accuracy of model prediction
(accuracy <- sum(diag(pred.mat1))/nrow(test) * 100) #79.19463% -- best model tied with k=5 below
(accuracy <- sum(diag(pred.mat5))/nrow(test) * 100) #79.19463%
(accuracy <- sum(diag(pred.mat10))/nrow(test) * 100) #77.18121%


#----------------- K-Means Clustering #####
set.seed(11)
diet.cluster <- kmeans(microbiome[,4:98],2)

table(diet.cluster$cluster, microbiome$Diet)
#      0   1
# 1  349 239
# 2    1   4

#This does not look promising. Kmeans is likely not going to be helpful in creating an estimation model.

#----------------- Random Forest ######

# went by the rule of sqrt(p variables) when building a random forest of classification trees
# sqrt(6701) = 81.86
set.seed(100)
rf.biome= randomForest(Diet~.,data=train, mtry=80, importance =TRUE)
rf.biome
#No. of variables tried at each split: 80
#OOB estimate of  error rate: 15.99%
#Confusion matrix:
#  0   1 class.error
#0 236  27   0.1026616
#1  44 137   0.2430939

importance (rf.biome)
varImpPlot (rf.biome)
# source is showing as the most significant, sex is showing little signficance
#the OTUs 77, 85, 71,54, 41, 9 look very significant
#Renyi looks more predictive than the ShannonIndex


rf.p <- predict(rf.biome,newdata=subset(test,select=c(2:98)),type='response')
rf.p
table(rf.p, test$Diet)
#         0  1
#       0 77 15
#       1 10 47
25/149
#16.78 % misclassification rate

sapply(c(is.vector, is.matrix, is.list, is.data.frame), do.call, list(rf.p))
rf.pr <- prediction(as.numeric(rf.p), as.numeric(test$Diet))

#ROC
rf.prf <- performance(rf.pr, measure = "tpr", x.measure = "fpr")
plot(rf.prf)

#AUC
#rf.pred = predict(rf.biome,type="prob")[, 2]
#rf.roc = prediction(rf.pred, train$Diet)
rf.auc <- as.numeric(performance(rf.pr, "auc")@y.values)
rf.auc
#0.821561


#RF with variables with Mean Decrease Accuracy ~10 or higher, smaller number of predictors changed mtry to 3
set.seed(102)
rf.biome2 = randomForest(Diet~ OTU77+ OTU85 + OTU71 + OTU54 + OTU41 + OTU9 + Source, data=train, mtry=3, importance =TRUE)
rf.biome2
#Number of trees: 500, No. of variables tried at each split: 3
#OOB estimate of  error rate: 16.22%
#Confusion matrix:
#    0   1 class.error
#0 235  28   0.1064639
#1  44 137   0.2430939

importance (rf.biome2)
varImpPlot (rf.biome2)
#Mean Decrease accuracy and Mean Decrease Gini went down across the board. 
#            0         1           MeanDecreaseAccuracy MeanDecreaseGini
#OTU77  12.594688 23.538983            23.954547         30.08362
#OTU85   7.630045 19.021186            17.401198         26.73102
#OTU71   7.566652  2.699498             7.485637         21.28655
#OTU54   2.810747  6.697529             6.214892         19.83674
#OTU41  42.446493 52.255283            60.191907         55.55406
#OTU9    9.702282  9.090364            13.092183         22.52913
#Source 30.358108 49.766697            51.591274         37.58660

rf.p2 <- predict(rf.biome2,newdata=subset(test,select=c(2:98)),type='response')
rf.p2
table(rf.p2, test$Diet)
#  0  1
#0 78 15
#1  9 47
(9+15)/149
# 16.107% Misclassification rate (minute improvement from full model)

#coming as false for all the below, use as.numeric to get the right format
sapply(c(is.vector, is.matrix, is.list, is.data.frame), do.call, list(rf.p2))
rf.pr2 <- prediction(as.numeric(rf.p2), as.numeric(test$Diet))

#rf.pr <- prediction(rf.p, test$Diet)
rf.prf2 <- performance(rf.pr2, measure = "tpr", x.measure = "fpr")
plot(rf.prf2)

#rf.pred2 = predict(rf.biome2,type="prob")[, 2]
#rf.roc2 = prediction(rf.pred2, train$Diet)
rf.auc2 <- as.numeric(performance(rf.pr2 , "auc")@y.values)
rf.auc2
#0.8330552

#RF model remvoving an regressors with less than 1 
importance    <- importance(rf.biome)
varImportance <- data.frame(Variables = row.names(importance), 
                            Importance = round(importance[ ,'MeanDecreaseGini'],2))

#Create a rank variable based on importance
rankImportance <- varImportance %>%
  mutate(Rank = paste0('#',dense_rank(desc(Importance))))
rankImportance

#RF with mean decrease gini < 1 removed, mtry updated to 7 to reflect change in p 
set.seed(103)
rf.biome3 <- randomForest(Diet~.-Sex - OTU70 - OTU55- OTU59 - OTU5 - OTU27 - OTU92 - OTU84 - OTU37 
                          - OTU21- OTU50 - OTU0 - OTU66- OTU86 - OTU11 - OTU35 - OTU32 - OTU76 - OTU38
                          - OTU46 - OTU65 - OTU53 - OTU93 - OTU49 - OTU52 - OTU89- OTU24- OTU44- OTU36
                          - OTU1-OTU48 - OTU58 - OTU94 - OTU16 - OTU17 - OTU2 - OTU13 - OTU19 - OTU42
                          - OTU51 - OTU64 - OTU56 - OTU91 - OTU45 - OTU61 - OTU62- OTU80 - OTU4, data = train, mtry=7, importance =TRUE )

rf.biome3
#No. of variables tried at each split: 7
#OOB estimate of  error rate: 18.02%
#Confusion matrix:
#    0   1 class.error
#0 240  23  0.08745247
#1  57 124  0.31491713

importance (rf.biome3)
varImpPlot (rf.biome3)
#decrease in MDA and MDG for most regressors


#check against test set
rf.p3 <- predict(rf.biome3,newdata=subset(test,select=c(2:98)),type='response')
rf.p3
table(rf.p3, test$Diet)
#   0  1
#0 81 20
#1  6 42
(6+20)/149
# 17.45% misclassification rate - a bit worse than the two other models. 


sapply(c(is.vector, is.matrix, is.list, is.data.frame), do.call, list(rf.p3))
rf.pr3 <- prediction(as.numeric(rf.p3), as.numeric(test$Diet))

#ROC
rf.prf3 <- performance(rf.pr3, measure = "tpr", x.measure = "fpr")
par(mfrow=c(2,2))
plot(rf.prf3)

#AUC
#rf.pred3 = predict(rf.biome3,type="prob")[, 2]
#rf.roc3 = prediction(rf.pred3, test$Diet)
rf.auc3 <- as.numeric(performance(rf.pr3, "auc")@y.values)
rf.auc3
#0.8042269


#---------------- Boosting ######

set.seed(14)
#The bernoulli distribution only produced NAs - does not like the factor, so created a new train set with diet as a character
train.boost <- train
train.boost$Diet <- as.character(train.boost$Diet)
boost.biome=gbm(Diet~.,data=train.boost, distribution="bernoulli",n.trees=5000, interaction.depth=5, shrinkage = 0.01)

boost.biome
#5000 iterations were performed.
#There were 97 predictors of which 0 had non-zero influence.
summary(boost.biome)
importance(boost.biome)
varImp(boost.biome,numTrees = 50)

par(mfrow=c(1,2))
plot(boost.biome ,i="Source") #to plot against significant regressors
plot(boost.biome ,i="OTU85")
plot(boost.biome ,i="OTU41")
plot(boost.biome ,i="OTU14")
plot(boost.biome ,i="OTU53")
plot(boost.biome ,i="OTU73")
plot(boost.biome ,i="OTU75")

#predict against test set
boost.p=predict(boost.biome ,newdata =subset(test,select=c(2:98)),n.trees=5000, type = 'response')
boost.p

#With Boost, the Youden's Index function wasn't working well, so we used F1 to choose a threshold instead
#F1 score
my_reference <- as.factor(test$Diet) #Save response var as factor to be used in f1_score_func function
f1_score_func <- function(threshold){ #Function inputs threshold and outputs threshold with f1 pairing as named vector
  rounded_preds <- as.factor(as.integer(boost.p >= threshold)) #Round predictions according to threshold
  my_F1 <- confusionMatrix(data = rounded_preds, reference = my_reference, mode = "prec_recall")[[4]][7] #Extract F1
  attributes(my_F1) <- NULL #Remove attribute (name)
  c("threshold" = threshold, "F1" = my_F1) #Return threshold and F1 pairing
}
f1_scores <- sapply(seq(from = 0.01, to = 0.5, by = 0.01), f1_score_func) #Create list of f1scores paired to thresholds
f1_scores <- as_tibble(t(f1_scores)) #Convert to tibble
plot(f1_scores$threshold, f1_scores$F1) #Inspect pattern. Appears to level off around threshold of .45


boost.p2 <- ifelse(boost.p > 0.45,1,0)
table(boost.p2, test$Diet)
#   0  1
#0 75 21
#1 12 41
(12+21)/149
#22.15% misclassification rate

sapply(c(is.vector, is.matrix, is.list, is.data.frame), do.call, list(boost.p))
boost.pr <- prediction(as.numeric(boost.p), as.numeric(test$Diet))
boost.auc <- as.numeric(performance(boost.pr , "auc")@y.values)
boost.auc
#0.8073786

#ROC curve
plot(performance(boost.pr, measure = "tpr", x.measure = "fpr"))
#Looks okay

#----------- Boosting model with smaller shrinkage number

set.seed(14)
boost.biom2=gbm(Diet~.,data=train.boost, distribution="bernoulli",n.trees=5000, interaction.depth=5, shrinkage = 0.0005)

boost.biom2
#5000 iterations were performed.
#There were 97 predictors of which 0 had non-zero influence.
summary(boost.biom2)
varImp(boost.biom2,numTrees = 50)

par(mfrow=c(1,2))
plot(boost.biome ,i="Source") #to plot against significant regressors
plot(boost.biome ,i="OTU6")#
plot(boost.biome ,i="OTU9")
plot(boost.biome ,i="OTU25")
plot(boost.biome ,i="OTU41")
plot(boost.biome ,i="OTU77")
plot(boost.biome ,i="OTU82")

#predict against test set
boos2.p=predict(boost.biom2 ,newdata=subset(test,select=c(2:98)),n.trees=5000, type = 'response')
boos2.p

#F1 score
my_reference <- as.factor(test$Diet) #Save response var as factor to be used in f1_score_func function
f1_score_func <- function(threshold){ #Function inputs threshold and outputs threshold with f1 pairing as named vector
  rounded_preds <- as.factor(as.integer(boos2.p >= threshold)) #Round predictions according to threshold
  my_F1 <- confusionMatrix(data = rounded_preds, reference = my_reference, mode = "prec_recall")[[4]][7] #Extract F1
  attributes(my_F1) <- NULL #Remove attribute (name)
  c("threshold" = threshold, "F1" = my_F1) #Return threshold and F1 pairing
}
f1_scores <- sapply(seq(from = 0.01, to = 0.5, by = 0.01), f1_score_func) #Create list of f1scores paired to thresholds
f1_scores <- as_tibble(t(f1_scores)) #Convert to tibble
plot(f1_scores$threshold, f1_scores$F1) #Inspect pattern. Appears to level off around threshold of .3

boos2b.p <- ifelse(boos2.p > 0.3,1,0)
table(boos2b.p, test$Diet)
#   0  1
# 67 10
# 20 52
(10+20)/149
#20.13% misclassification rate

sapply(c(is.vector, is.matrix, is.list, is.data.frame), do.call, list(boos2.p))
boos2.pr <- prediction(as.numeric(boos2.p), as.numeric(test$Diet))
boos2.auc <- as.numeric(performance(boos2.pr , "auc")@y.values)
boos2.auc
#0.8516871

#ROC curve
plot(performance(boos2.pr, measure = "tpr", x.measure = "fpr"))
#This appears to be very good

#----------------------- Next Boosting Model with cross validation
fitControl <- trainControl(method = "cv", number = 10 ) #5folds)
tune_Grid <-  expand.grid(interaction.depth = 2, n.trees = 2500, shrinkage = 0.01, n.minobsinnode = 10)
set.seed(108)
fit <- train(Diet ~ ., data = train.boost, method = "gbm", trControl = fitControl, verbose = FALSE, tuneGrid = tune_Grid)

fit
#Stochastic Gradient Boosting 
#444 samples
#97 predictor
#2 classes: '0', '1' 

#No pre-processing
#Resampling: Cross-Validated (10 fold) 
#Summary of sample sizes: 400, 399, 400, 399, 400, 399, ... 
#Resampling results:
#Accuracy   Kappa   
#0.7725033  0.5113699

summary(fit)
varImp(fit,numTrees = 50)
#         Overall
#OTU41   100.0000
# Source4  57.7576
# OTU85     5.6704
# OTU6      4.7223
# OTU9      2.4721
# OTU8      2.4475
# OTU77     2.3635
# OTU54     2.1988
# OTU39     1.3119
# All others = 0.0000

#try against the test data
fit.p= predict(fit ,newdata =subset(test,select=c(2:98)),n.trees=5000, type = 'raw')
fit.p #Notes, this output is all 1s and 0s, not probabilities

table(fit.p, test$Diet)
#   0  1
#0 73 23
#1 14 39
(15+23)/149
#25.5% Misclassification rate

sapply(c(is.vector, is.matrix, is.list, is.data.frame), do.call, list(fit.p))
fit.pr <- prediction(as.numeric(fit.p), as.numeric(test$Diet))
fit.auc <- as.numeric(performance(fit.pr , "auc")@y.values)
fit.auc
#0.7340564

#ROC curve
plot(performance(fit.pr, measure = "tpr", x.measure = "fpr"))
#Very similar to original "full model" binomial logistic regression
