#Microbiome Project

library(tidyverse) # Load the core tidyverse packages: ggplot2, tibble, 
# tidyr, readr, purrr, and dplyr

setwd("~/Documents/GitHub/DataMining_MicrobiomeProject")

microbiome <- read.csv("MicrobiomeWithMetadata.csv", encoding = 'utf-8', stringsAsFactors = FALSE)

#----------------- Grading Criteria
#Data exploration
#Data cleaning (missing values, outliers, etc.)
#Rationale for the selected statistical modeling methods
#Correct implementation and use of statistical modeling methods
#Appropriate model selection approach (train/test, cross-validation, etc.)
#Thoroughly documented code (R comments)
#Appropriate use of functions to reduce code complexity and redundancy
#Writing quality for final report, evaluated in terms of conformance to process outline, level of detail, and correctness.

#----------------- K-Nearest Neighbors (KNN)
library(class) #package for KNN model
biome.mat <- as.data.frame(data.matrix(microbiome)) #make matrix a data frame

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
#Predictions   0   1   2   3   4   5
#              0 389   2   0   0   0   0
#              1   0 267   0   0   0   0
#              3   0   0   1   1   0   0
#              4   0   0   0   0   9   0
#              5   0   0   0   0   0   6

pred.mat5 <- table("Predictions" = knn.pred5, Actual = cl[test]); pred.mat5
#               Actual
#Predictions   0   1   2   3   4   5
#              0 370   0   0   0   0   0
#              1  19 269   1   1   0   0
#              3   0   0   0   0   0   0
#              4   0   0   0   0   9   0
#              5   0   0   0   0   0   6

pred.mat10 <- table("Predictions" = knn.pred10, Actual = cl[test]); pred.mat10
#               Actual
#Predictions   0   1   2   3   4   5
#              0 358   0   0   0   0   0
#              1  31 269   1   1   0   0
#              3   0   0   0   0   0   0
#              4   0   0   0   0   5   5
#              5   0   0   0   0   4   1

#assess accuracy of model prediction
(accuracy <- sum(diag(pred.mat1))/length(test) * 100) #97.33333%
(accuracy <- sum(diag(pred.mat5))/length(test) * 100) #94.66667%
(accuracy <- sum(diag(pred.mat10))/length(test) * 100) #93.48148%

#----------------- K-Means Clustering
set.seed(1)
diet.cluster <- kmeans(microbiome[,6:94], 6)
diet.cluster

table(diet.cluster$cluster, microbiome$Diet)
#    0   1   2   3   4   5
#1   2   6   0   0   0   0
#2  61  10   0   0   0   0
#3 321 247   1   1   9   6
#4   3   0   0   0   0   0
#5   2   4   0   0   0   0
#6   0   2   0   0   0   0

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

