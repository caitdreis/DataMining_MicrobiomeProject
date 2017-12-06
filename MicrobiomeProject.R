#Microbiome Project

#----------------- Packages
library(tidyverse) # Load the core tidyverse packages: ggplot2, tibble, 
# tidyr, readr, purrr, and dplyr
library(vegan) #calculation of diversity metrics
library(psych) #descriptive statistics
library(class) #package for KNN model
library(randomForest) #package for Random Forest model
library(pscl) #For logistic regression R^2
library(ROCR) #For ROC curves

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


#----------------- Training and Testing Sets ######

#Cross validation of testing data with testing data entered in for test
set.seed(17)
train = sample(1:nrow(microbiome), nrow(microbiome) * .75)
test <- microbiome[-train, c(1:98)]
train <- microbiome[train, c(1:98)]

microbiome <- microbiome[, c(1:98)]

#----------------- Multinomial Logistic Regression ######

model <- glm(Diet~.,family=binomial(link='logit'),data=train)
summary(model) #Every variable is significant...
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
fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError <- mean(fitted.results != test$Diet)
print(paste('Accuracy',1-misClasificError))
#"Accuracy 0.751677852348993"

p <- predict(model,newdata=subset(test,select=c(2:98)),type='response')
p
pr <- prediction(p, test$Diet)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc #0.7410085

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

#This does not look promising. Kmeans is likely not going to be helpful in creating a prediction model.


#----------------- Tree Based Methods ######
#----------------- Random Forest

# went by the rule of ???p variables when building a random forest of classification trees
# sqrt(6701) = 81.86
rf.biome= randomForest(Diet~.,data=train, mtry=80, importance =TRUE)
yhat.rf = predict(rf.biome,newdata=test)
mean((yhat.rf-test)^2)
#[1] 0.3987529
importance (rf.biome)
# nothing is showing as significantly important



