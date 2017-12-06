#Microbiome Project

#----------------- Packages
library(tidyverse) # Load the core tidyverse packages: ggplot2, tibble, 
# tidyr, readr, purrr, and dplyr
install.packages("vegan") 
library(vegan) #calculation of diversity metrics
library(psych) #descriptive statistics

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
#This is likely coincidence

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
microbiome$ShannonIndex <- diversity(microbiome[,6:94])
summary(microbiome$ShannonIndex)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000001 0.279000 0.693300 1.218000 1.161000 4.365000 

#Reyyi Index
#specnumber fucntion finds the number of species in the sample
specnumber(microbiome[,6:94]) #89 species present in each sample
h <- diversity(microbiome[,6:94])
J <- h/log(specnumber(microbiome[,6:94])) #assessing Pielou's eveness

#now take a random subset of six sites to assess Renyi diversity index
k <- sample(nrow(microbiome[,6:94]), 6) 
R <- renyi(microbiome[,6:94][k,])
plot(R) #can visualize these six sites, labels sites chosen

#----------------- Further Descriptive and ANOVA Exploration with Featured Variable 
describeBy(microbiome$ShannonIndex, microbiome$Diet)
#Descriptive statistics by group 
#group: 0
#vars   n mean   sd median trimmed  mad min  max range skew kurtosis   se
#X1    1 350 1.25 1.45   0.69    1.02 0.64   0 4.36  4.36 1.43     0.51 0.08
#-------------------------------------------------------------------------- 
#  group: 1
#vars   n mean   sd median trimmed  mad min  max range skew kurtosis   se
#X1    1 243 1.14 1.34   0.71    0.89 0.75   0 4.36  4.36 1.63     1.41 0.09

#significance testing with ANOVA
fit <- aov(microbiome$ShannonIndex ~ microbiome$Diet, data=microbiome)
plot(fit)
summary(fit)
#               Df Sum Sq Mean Sq F value Pr(>F)
#microbiome$Diet   1    1.7   1.712   0.869  0.352
#Residuals       591 1164.2   1.970    
#No significant on the 

#posthoc TukeyHSD (Honestly Significant Differences)
TukeyHSD(fit) #where fit comes from aov()
#Tukey multiple comparisons of means
#95% family-wise confidence level
#Fit: aov(formula = microbiome$ShannonIndex ~ microbiome$Diet, data = microbiome)
#$`microbiome$Diet`
#       diff        lwr       upr     p adj
#1-0 -0.1092519 -0.3394191 0.1209153 0.3515974
#Tukey's posthoc testing tells us that there is truly no significance in diet's
#effect on the mean diversity index measure.

#significance testing with ANOVA using sex
fit.sex <- aov(microbiome$ShannonIndex ~ microbiome$Sex, data=microbiome)
plot(fit.sex)
summary(fit.sex)
#               Df Sum Sq Mean Sq F value Pr(>F)
#microbiome$Sex   1      0  0.0022   0.001  0.974
#Residuals      591   1166  1.9727 

#posthoc TukeyHSD 
TukeyHSD(fit.sex) 
#         diff        lwr       upr     p adj
#1-0 0.006854582 -0.4008207 0.4145299 0.9736681

#----------------- Multinomial Logistic Regression ######




#----------------- K-Nearest Neighbors (KNN) ######
library(class) #package for KNN model
biome.mat <- microbiome #make new dataframe to use for KNN only

#Cross validation of testing data with testing data entered in for test
train <- sample(nrow(biome.mat), ceiling(nrow(biome.mat) * .50))
test <- (1:nrow(biome.mat))

cl <- biome.mat[, "Diet"] #define the classifier variable of sentiment

#create model with training data, test data, and training set classifier
knn.pred1 <- knn(biome.mat[train, ], biome.mat[test, ], cl[train])
knn.pred5 <- knn(biome.mat[train, ], biome.mat[test, ], cl[train], k=5)
knn.pred10 <- knn(biome.mat[train, ], biome.mat[test, ], cl[train], k=10)

#create a prediction table to understand how classes of sentiment values are being defined
pred.mat1 <- table("Predictions" = knn.pred1, Actual = cl[test]); pred.mat1
#               Actual
#Predictions       0   1   2   3   4   5
#              0 389   2   0   0   0   0
#              1   0 267   0   0   0   0
#              3   0   0   1   1   0   0
#              4   0   0   0   0   9   0
#              5   0   0   0   0   0   6

pred.mat5 <- table("Predictions" = knn.pred5, Actual = cl[test]); pred.mat5
#               Actual
#Predictions       0   1   2   3   4   5
#              0 370   0   0   0   0   0
#              1  19 269   1   1   0   0
#              3   0   0   0   0   0   0
#              4   0   0   0   0   9   0
#              5   0   0   0   0   0   6

pred.mat10 <- table("Predictions" = knn.pred10, Actual = cl[test]); pred.mat10
#               Actual
#Predictions       0   1   2   3   4   5
#              0 358   0   0   0   0   0
#              1  31 269   1   1   0   0
#              3   0   0   0   0   0   0
#              4   0   0   0   0   5   5
#              5   0   0   0   0   4   1

#assess accuracy of model prediction
(accuracy <- sum(diag(pred.mat1))/length(test) * 100) #97.33333%
(accuracy <- sum(diag(pred.mat5))/length(test) * 100) #94.66667%
(accuracy <- sum(diag(pred.mat10))/length(test) * 100) #93.48148%

#----------------- K-Means Clustering #####
set.seed(1)
diet.cluster <- kmeans(microbiome[,6:94], 6)
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

#Recently performed studies in rodents have indicated that Akkermansia muciniphila in 
#the intestinal tract may reduce obesity, diabetes, and inflammation.

#----------------- K-Means Clustering from Scratch
#need to make a function for euclidean distance from scratch
eu.dist <- function(x1, x2) {
  dist.m <- matrix(NA, nrow=dim(x1)[1], ncol=dim(x2)[1])
  for(i in 1:nrow(x2)) { dist.m[,i] <- sqrt(rowSums(t(t(x1)-x2[i,])^2))}
  dist.m
}

cluster.kmeans <- function(x, centroids, distance, num) {
  cluster.list <- vector(num, mode="list"); center.list <- vector(num, mode="list")
  for(i in 1:num) { #start the loop through the number of specified iterations from num
    dist <- distance(x, centroids) #take the distance between the point and the calculated centroids
    clusters <- apply(dist, 1, which.min) #apply the which.min function to the distance
    centroids <- apply(x, 2, tapply, clusters, mean) #apply the mean function to create the clusters
    cluster.list[[i]] <- clusters #fill the list with the calculated clusters
    center.list[[i]] <- centroids } 
  list(clusters=cluster.list, centroids=center.list) 
}

test = microbiome[, 6:94] #make a test set with OTU columns 
test.m = as.matrix(test)
centroid <- test.m[sample(nrow(test.m), 10),] #take the centroids
results <- cluster.kmeans(test.m, centroid, eu.dist, 10); results


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



