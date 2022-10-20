###############################################################################
# PEN-iRF-LOOP
# Create Predictive Expression Network using Iterative Random Forest Leave One
# out Prediction
# Author: Amir Valizadeh
###############################################################################

parse_arguments <- function() {
  suppressPackageStartupMessages(library(optparse))
  
  option_list = list(
    make_option(c("-d", "--data"),
                action = "store",
                default = NULL,
                type = 'character',
                help = "Path to the .Rdata file containing an ExpressionSet object with a count matrix for RNA-Seq expression values with features as rows and samples as columns."),
    make_option(c("-i", "--iter"),
                action = "store",
                default = "5",
                help = "The number of iterations to run iRF."),
    make_option(c("-k", "--keep"),
                action = "store",
                default = "10",
                help = "The top 'k' percent of each column of network to keep.")
  )
  
  desc <- "PEN-iRF-LOOP.R"
  opt <- parse_args(OptionParser(option_list = option_list, description = desc),
                    convert_hyphens_to_underscores = TRUE)
  
  errors <- 0
  # Check whether data have been set by the user
  if(is.null(opt$data)) {
    message("ERROR:: --data is required but is not set.")
    errors <- errors + 1
  } else if(!file.exists(opt$data)) {
    message("ERROR:: --data must be an existing Rdata file.")
    errors <- errors + 1
  }
  
  if(errors > 0) {
    quit()
  } else {
    return(opt)
  }
}

########################################################################
# Functions
########################################################################

load_data = function(opt) {
  suppressPackageStartupMessages(library(Biobase))
  suppressPackageStartupMessages(library(iRF))
  
  if (is.null(opt$data)) {
    return(NULL)
  } else {
    assign('data', get(load(opt$data)))
    exprs_data <- exprs(data)
  }
    
  rm(data)
  return(exprs_data)
}

get_or_set_iter = function(opt) {
  # If arguments not passed, set iter to 5
  if (is.null(opt$iter)) {
    iter <- 5
    return(iter)
  }
  # If arguments is passed, set iter to the argument
  if (!is.null(opt$iter)) {
    iter <- opt$iter
    return(iter)
  }
  return(iter)
}

get_or_set_k = function(opt) {
  # If arguments not passed, set k to 10
  if (is.null(opt$keep)) {
    k <- 10
    return(k)
  }
  # If arguments is passed, set k to the argument
  if (!is.null(opt$keep)) {
    k <- opt$keep
    return(k)
  }
  return(k)
}

########################################################################
# iRF-LOOP
########################################################################

main <- function() {
  # Parse command line arguments
  opt <- parse_arguments()
  
  # Load arguments
  exprs_data <- load_data(opt)
  iter <- as.numeric(get_or_set_iter(opt))
  k <- as.numeric(get_or_set_k(opt))
  
  # Set other arguments
  imp_matrix <- c()
  n_trees <- round(sqrt(dim(exprs_data)[1]), digits = 0)
  
  # iRF-LOOP
  for (i in 1:dim(exprs_data)[1]) {
    Y = matrix(exprs_data[i, ],
      ncol = 1,
      dimnames = list(colnames(exprs_data), rownames(exprs_data)[i]))
    X = t(exprs_data)
    imp_score <- iRF(X, Y, n.iter = iter, ntree = n_trees, verbose = F, n.core = -1)
    imp_score <- imp_score[[2]]
    colnames(imp_score) <- rownames(exprs_data)[i]
    imp_score[i, ] <- 0
    imp_score <- (imp_score - min(imp_score)) / (max(imp_score) - min(imp_score))
    imp_matrix <- cbind(imp_matrix, imp_score)
    if (i %% round((dim(exprs_data)[1] / 10), digits = 0) == 0) {
      print(paste(i, "of", dim(exprs_data)[1], "features done!"))
    }
  }
  
  # Keep the top 'k' percent
  for (i in 1:dim(exprs_data)[1]) {
    imp_matrix[, i] <- sort(imp_matrix[, i], decreasing = T)
    thresh <- imp_matrix[round(length(imp_matrix[, i]) / k, digits = 0), i]
    imp_matrix[which(imp_matrix[, i] <= thresh), i] <- NA
  }
  
  return(imp_matrix)
}

########################################################################
# Save the network and finish
########################################################################

main()
save(imp_matrix, "PEN.RData")
print("PEN Created")
quit(save = 'no', status = 0)