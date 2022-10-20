###############################################################################
# The process of iRF-LOOP for one feature column is presented in here.
# The aim of this file is for educational purposes.
###############################################################################

# The following file contains RNA-seq data for 22 samples of different human tissues and 52580 genes.
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

# Set the importance score of Y for predicting itself as 0:
imp_score <- rbind(0, imp_score[[2]])
rownames(imp_score)[1] <- rownames(exprs_data)[1]

# Check the highest value:
imp_score[which(imp_score == max(imp_score)), ]

# Let's do min-max feature scaling:
imp_score <- (imp_score - min(imp_score)) / (max(imp_score) - min(imp_score))
