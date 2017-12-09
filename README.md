# DataMining_MicrobiomeProject

## Team Members

Elizabeth Homan | eih2nn@virginia.edu 
Jennifer Cruser | jc4pg@virginia.edu
Caitlin Dreisbach | cnd2y@virginia.edu

## Delegation of Project Tasks

Elizabeth Homan | Logistic regression, ROC curves/evaluation metrics, data preprocessing, approach, and analysis.

Jennifer Cruser | Random Forest, boosted models, analysis, and interpretation.

Caitlin Dreisbach | Github owner, descriptive t-tests, KNN, background research, and OTU interpretation.

## Comparison of Models and Evaluation Metrics

# Type of Model: Logistic Regression

Best Evaluation Metrics: AUC: 0.82

Accuracy: 79%

Details on Best Model Manually selected variables, 26 predictors selected by p-value.

Model Rank in Performance: 3

# K-Nearest Neighbors

AUC: 0.78

Accuracy: 79%

KNN with k=1 (noted to have the same accuracy as k=5)

Model Rank in Performance: 4


# Random Forest

AUC: 0.83

Accuracy: 85%

Fewer predictors (OTU: 9, 41, 54, 71, 77, 85)

Model Rank in Performance: 1


# Boosted 

AUC: 0.85

Accuracy: 80%

Bernoulli distributed boosted model with shrinkage tuning, previously attempted with 5-fold CV

Model Rank in Performance: 2
 
