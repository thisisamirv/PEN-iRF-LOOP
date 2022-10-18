# The following file contains RNA-seq data for 22 samples of different human tissues and 52580 genes:
# Import data
library(Biobase)
load("Example.RData")
exprs_data <- exprs(Example)
rm(Example)

# Sample names include: SRX003935, SRX003921, ...
# Gene names include: ENSG00000000003, ENSG00000000005, ...
# Let's take a peak:
dim(exprs_data)
exprs_data[1:5, 1:3]

# m = 22
# n = 52580

# Now, we will create an iRF with the expression values for the first gene (ENSG00000000003) as Y and the expression values for the rest of the 52579 genes as predictor matrix X of size 22 * 52579:
Y = matrix(exprs_data[1, ],
  ncol = 1,
  dimnames = list(colnames(exprs_data), rownames(exprs_data)[1]))
X = t(exprs_data[c(-1), ])

# Take a peak
X[1:3, 1:3]
Y[1:3, ]

# Now we perform iRF:
imp_score <- iRF(X, Y, n.iter = 5, ntree = 229, verbose = F, n.core = -1)

# Set the importance score of Y for predicting itself as 0
imp_score <- rbind(0, imp_score[[2]])
rownames(imp_score)[1] <- rownames(exprs_data)[1]

# Check the highest value
imp_score[which(imp_score == max(imp_score)), ]

# Let's do min-max feature scaling:
imp_score <- (imp_score - min(imp_score)) / (max(imp_score) - min(imp_score))

# Now, we will define a loop to do the same for all genes:
imp_matrix <- c()
n_trees <- round(sqrt(dim(exprs_data)[1]) , digits = 0)

for (i in 1:dim(exprs_data)[1]) {
  Y = matrix(exprs_data[i, ],
    ncol = 1,
    dimnames = list(colnames(exprs_data), rownames(exprs_data)[i]))
  X = t(exprs_data)
  imp_score <- iRF(X, Y, n.iter = 5, ntree = n_trees, verbose = F, n.core = -1)
  imp_score <- imp_score[[2]]
  colnames(imp_score) <- rownames(exprs_data)[i]
  imp_score[i, ] <- 0
  imp_score <- (imp_score - min(imp_score)) / (max(imp_score) - min(imp_score))
  imp_matrix <- cbind(imp_matrix, imp_score)
  if (i %% 1000 == 0) {
    print(paste(i, "of", dim(exprs_data)[1], "Completed"))
  }
}

# Keep the top *n* percent of each column:
n = 10
for (i in 1:dim(exprs_data)[1]) {
  imp_matrix[, i] <- sort(imp_matrix[, i], decreasing = T)
  thresh <- imp_matrix[round(length(imp_matrix[, i]) / n, digits = 0), i]
  imp_matrix[which(imp_matrix[, i] <= thresh), i] <- NA
}

# Save the network:
save(imp_matrix, "PEN.RData")
