#Microbiome Project

#----------------- Packages
library(tidyverse) # Load the core tidyverse packages: ggplot2, tibble, 
# tidyr, readr, purrr, and dplyr
install.packages("vegan") 
library(vegan) #calculation of diversity metrics
library(psych) #descriptive statistics
library(class) #package for KNN model

#----------------- Working Directory
setwd("~/Documents/GitHub/DataMining_MicrobiomeProject")

#----------------- Read in Data
microbiome <- read.csv("MicrobiomeWithMetadata.csv", encoding = 'utf-8', stringsAsFactors = FALSE)
#View(microbiome)

#----------------- Grading Criteria
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
#This dataset was pre-curated from the original Science article and is thus relatively clean to begin with
library(psych)

nrow(microbiome) #675 observations

sum(is.na(microbiome)) #No missing data requiring imputation

colnames(microbiome[1:5])
#"Diet" "Source" "Donor" "CollectionMet" "Sex"  
#Beyond this, all column names are OTU groups

table(microbiome$Donor)
#0   1   2   3   4   5   6 
#313  64   6 151  92  46   3 

#   Donor	
# 0 - HMouseLFPP -- mice
# 1 - CONVR
# 2 - Human
# 3 - Fresh -- mice
# 4 - Frozen -- mice
# 5 - HMouseWestern -- mice
# 6 - CONVD

#We will be using only non-human (mouse) donors for our analysis
microbiome <- subset(microbiome, Donor == 0 | Donor == 3 | Donor == 4 | Donor == 5)

nrow(microbiome) #602 observations

#Further explore the distribution of these variables
table(microbiome$Diet)
# 0   1   4 
# 350 243 9 

#   Diet	
#0	LFPP: low fat, high-plant polysaccharide diet
#1	Western: high-fat, high-sugar Western diet
#4	Suckling

#Remove all rows with suckling diets (want to compare only high and low fat diets)
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

table(microbiome$Sex)
#  0   1 
# 543  50  -- Uneven distribution between males and females (see below for reference)

#Sex
#0 - Male
#1 - Female

plot(Diet~Sex, data=microbiome)
#Females appear to be more unevenly split between diets (favoring the low fat diet)

#   Diet	
#0	LFPP: low fat, high-plant polysaccharide diet
#1	Western: high-fat, high-sugar Western diet

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
#Shannon Diversity Index
microbiome$ShannonIndex <- NULL
microbiome$ShannonIndex <- diversity(microbiome[,4:98])
summary(microbiome$ShannonIndex)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000001 0.348860 0.693277 1.213849 1.168294 4.432951

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
# Min.  1st Qu.   Median   Mean  3rd Qu.     Max. 
# 0.0000  0.1216  0.5273  0.9932  0.9235  4.0791 

#Is there a significant difference between these indexes?
comp.index <- aov(microbiome$ShannonIndex ~ microbiome$Renyi, data=microbiome)
plot(comp.index)
summary(comp.index)
#                 Df Sum Sq Mean Sq F value Pr(>F)    
#microbiome$Renyi   1 1183.2    1183   37225 <2e-16 ***
#Residuals        591   18.8       0  

#----------------- Further Descriptive and ANOVA Exploration with Featured Variables 
#----------------- Shannon Index
describeBy(microbiome$ShannonIndex, microbiome$Diet)
#Descriptive statistics by group 
#group: 0
#   vars   n mean   sd median trimmed  mad min  max range skew kurtosis   se
# X1    1 350 1.26 1.47   0.69    1.03 0.62   0 4.43  4.43 1.44     0.53 0.08
# --------------------------------------------------------------------------------------------- 
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
# --------------------------------------------------------------------------------------------- 
#   group: 1
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

#----------------- Multinomial Logistic Regression ######




#----------------- K-Nearest Neighbors (KNN) ######
biome.mat <- microbiome #make new dataframe to use for KNN only

#Cross validation of testing data with testing data entered in for test
train <- sample(nrow(biome.mat), ceiling(nrow(biome.mat) * .50))
test <- (1:nrow(biome.mat))

cl <- biome.mat[, "Diet"] #define the classifier variable of diet
modeldata <- biome.mat[,!colnames(biome.mat) %in% "Diet"]

#create model with training data, test data, and training set classifier
knn.pred1 <- knn(modeldata[train, ], modeldata[test, ], cl[train])
knn.pred5 <- knn(modeldata[train, ], modeldata[test, ], cl[train], k=5)
knn.pred10 <- knn(modeldata[train, ], modeldata[test, ], cl[train], k=10)

#create a prediction table to understand how classes of sentiment values are being defined
pred.mat1 <- table("Predictions" = knn.pred1, Actual = cl[test]); pred.mat1
#             Actual
#Predictions   0   1
#           0 301  60
#           1  49 183

pred.mat5 <- table("Predictions" = knn.pred5, Actual = cl[test]); pred.mat5
#             Actual
#Predictions   0   1
#           0 299 127
#           1  51 116

pred.mat10 <- table("Predictions" = knn.pred10, Actual = cl[test]); pred.mat10
#               Actual
#Predictions   0   1
#           0 295 131
#           1  55 112

#assess accuracy of model prediction
(accuracy <- sum(diag(pred.mat1))/length(test) * 100) #81.61889%
(accuracy <- sum(diag(pred.mat5))/length(test) * 100) #69.98314%
(accuracy <- sum(diag(pred.mat10))/length(test) * 100) #68.63406%


#----------------- K-Means Clustering #####
set.seed(1)
diet.cluster <- kmeans(microbiome[,4:98], 6)
diet.cluster

table(diet.cluster$cluster, microbiome$Diet)
#    0   1   2   3   4   5
#0   2   6   0   0   0   0
#1  61  10   0   0   0   0
#2 321 247   1   1   9   6
#3   3   0   0   0   0   0
#4   2   4   0   0   0   0
#5   0   2   0   0   0   0

#visualize clusters with ggplot2
library(ggplot2)
diet.cluster$cluster <- as.factor(diet.cluster$cluster)
#ggplot(microbiome, color = diet.cluster$cluster)) + geom_point()


#----------------- Tree Based Methods ######
#----------------- Random Forest
library(randomForest)

#Cross validation of testing data with testing data entered in for test
set.seed(17)
train = sample(1:nrow(microbiome), nrow(microbiome) * .75)
test <- microbiome[-train, ]
train <- microbiome[train, ]


# went by the rule of ???p variables when building a random forest of classification trees
# sqrt(6701) = 81.86
rf.biome= randomForest(Diet~.,data=train, mtry=80, importance =TRUE)
yhat.rf = predict(rf.biome,newdata=test)
mean((yhat.rf-test)^2)
#[1] 0.3987529
importance (rf.biome)
# nothing is showing as significantly important



