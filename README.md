# Introduction
Gene Regulatory Networks (GRNs) are crucial to understanding biological systems since they show the regulatory relationships between transcription factors and target genes. A conceptual expansion on the GRN is the Predictive Expression Network (PEN). Algorithms that produce PENs utilize all genes in the input data, creating an all-to-all genes network.

# RF as a feature selection method
Some feature selection methods:

* Pearson Correlation
* Mutual Information (MI)
* Sequential Feature Selection (SFS)
* Lasso
* Ridge Regression
* Random Forest (RF): not a classical feature selection method

An RF is an ensemble of decision trees (DT). The number of trees in a forest is a parameter that is chosen by the user. During the training phase, DTs are built, where a random subset of features for a random subset of samples are examined at each decision point and the one that best divides the data is chosen. Once a forest has been generated, the importance of each feature can be calculated from node impurity using one of the following methods:

* Gini index: for classification
* Variance explained: for regression
* Permutation importance

Because of the nature of DTs, the importance of any chosen feature is inherently conditional on the previously chosen features. In this way, RF can account for some of the interconnected dependencies that occur in biological systems.

# Iterative RF (iRF)
It is an algorithmic advancement of RF, which produces a more accurate model by iteratively creating weighted forests. First, an RF is created where features are unweighted and have an equal chance of being randomly sampled at any given node. The resulting importance scores for the features are then used for weighting the features in the next forest, thus increasing the chance that important features are evaluated at any given node. This process is repeated *i* times.

![Figure 1](https://github.com/thisisamirv/PEN/blob/main/Figure%201.jpeg)

**Note:** Due to the ability to easily follow the decisions that these models make, they have been named explainable AI (X-AI).

## Implementation of iRF in R
**iRF** is an R package for iRF implementation. Unfortunately, this implementation is not very effective with big data.

## Implementation of iRF in C++
**Ranger-Based iRF (RB-iRF)** is an implementation of iRF in C++. **Ranger** is an open-source RF implementation in C++. RB-iRF has used *ranger* as the core of its iRF implementation.

# iRF-LOOP
iRF Leave One out Prediction (iRF-LOOP) is a method for the creation of PENs on the order of 40,000 genes or more. The method is as follows:

1. Given a gene expression matrix of *m* samples (rows) and *n* genes (columns), iRF-LOOP starts by treating each gene as the dependent variable (*Y*) and the remaining *n - 1* genes as predictor matrix (*X*) of size *m x (n - 1)*.
2. Using an iRF model, the importance score of each gene in *X*, for predicting *Y*, is calculated. The result is a vector, of size *n*, of importance scores (the importance score of *Y*, for predicting itself, has been set to zero).
3. This process is repeated for each of the *n* genes, requiring *n* iRF runs.
4. The *n* vectors of importance scores are concatenated into an *n x n* importance matrix.
5. To keep importance scores on the same scale across the importance matrix, each column is normalized relative to the sum of the column.
6. Finally, from the importance matrix we will generate a network. Generally, the scores (weights) are thresholded at some value and only edges with large enough weights are included in the final network. For example, we may produce four thresholded networks, keeping the top 10%, 5%, 1%, and 0.1% of edge scores, respectively.

![Figure 2](https://github.com/thisisamirv/PEN/blob/main/Figure%202.jpeg)

**Note:** A common setting for the number of trees is the square root of the number of genes.

# Test
Run the following code in shell to test the PEN-iRF-LOOP:

  ./test.sh

# References
1- Cliff, A., Romero, J., Kainer, D., Walker, A., Furches, A., & Jacobson, D. (2019). A high-performance computing implementation of iterative random forest for the creation of predictive expression networks. Genes, 10(12), 996.
2- Walker, A. M., Cliff, A., Romero, J., Shah, M. B., Jones, P., Gazolla, J. G. F. M., ... & Kainer, D. (2022). Evaluating the performance of random forest and iterative random forest based methods when applied to gene expression data. Computational and Structural Biotechnology Journal, 20, 3372-3386.
