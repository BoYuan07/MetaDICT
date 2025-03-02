#' @title Data Check
#' 
#' @description Check the format of inputs.
#' 
#' @param count The integrated count table of taxa by samples. The
#' \code{count} parameter should be provided as either a \code{matrix} or a
#' \code{data.frame}.
#' @param meta The integrated meta table \code{meta} contains sample information
#' and batch id.
#' @param covariates The covariates used in data integration. Default is all.
#' @param distance_matrix The distance matrix that measures sequence dissimilarity.
#' @param tree The phylogenetic tree (optional if distance matrix or taxonomy is provided).
#' @param taxonomy The taxonomy table (optional if distance matrix or phylogenetic tree is provided).
#' @param tax_level The taxonomic level of count table.
#'
#' @return a \code{list} contains count list, meta table list, sequencing distance
#' matrix and parameters.
#' 
#' @export 
#' 
#' 
data_check <- function(count, meta, covariates = "all", distance_matrix = NULL, tree = NULL, taxonomy = NULL, tax_level = NULL, verbose = TRUE){
    #=================check input data type===================#
    if (verbose){
        message("Checking data input...")
    }
    
    if (!inherits(count, c("data.frame","matrix"))){
        stop("The integrated count table should be either a data frame
        or a matrix object.")
    }

    # check non-numeric value
    non_numeric_columns = vapply(count, function(col) !is.numeric(col), logical(1))
    if (any(non_numeric_columns)){
        stop("All columns must be numeric.")
    }
    taxa_name = rownames(count)
    if(is.null(taxa_name)){
        stop("Please set the rownames of count table as taxa names.")
    }
    
    # check taxa
    if (!is.null(distance_matrix)){
        if (nrow(distance_matrix)!=length(taxa_name)){
            stop("Taxa number does not match with distance matrix.")
        }
        dist_mat <- distance_matrix
        beta <- 0.1
        neighbor <- 5
        sigma <- median(as.vector(dist_mat))
    }else if (!is.null(tree)){
        message("A phylogenetic tree is provided. Sequencing similarity will be 
        estimated using phylogenetic tree.")
        if (! requireNamespace("ape", quietly = TRUE)){
            stop("The 'ape' package is needed to process the phylogenetic tree object.")
        }
        if (! all(taxa_name %in% tree$tip.label)){
            stop("Not all taxa are founded on the tree.")
        }
        drop_tips <- tree$tip.label[!(tree$tip.label %in% taxa_name)]
        tree <- ape::drop.tip(tree,drop_tips) %>% ape::as.phylo()
        dist_mat <- as.matrix(ape::cophenetic.phylo(tree))
        beta <- 0.1
        neighbor <- 5
        sigma <- median(as.vector(dist_mat))
    }else if(!is.null(taxonomy)){
        "A taxonomy table is provided. Sequencing similarity will be 
        estimated using taxonomic information."
        if (!all(taxa_name %in% rownames(taxonomy))){
            stop("Taxa names do not match between count table and taxonomy table.")
        }
        if (!(tax_level %in% colnames(taxonomy))){
            stop("The provided taxonomy level does not match with colnames of taxonomy.")
        }
        # Find the column index of taxonomic level
        col_index <- which(colnames(taxonomy) == tax_level)

        # generate adjacency matrix
        adj_mat <- matrix(0, nrow = length(taxa_name), ncol = length(taxa_name))
        taxonomy_subset <- taxonomy[taxa_name, 1:(col_index-1)]  # Extract columns before col_index
        adj_mat <- as.matrix(
            outer(1:length(taxa_name), 1:length(taxa_name), 
            Vectorize(function(i, j) as.numeric(all(taxonomy_subset[i, ] == taxonomy_subset[j, ])))
            )
        )
        dist_mat <- 1-adj_mat
        beta <- 0.01
        neighbor <- max(rowSums(dist_mat == 0))
        sigma <- 1
    }else{
        stop("Neither a phylogenetic tree nor a taxonomy table is proveded.
        Please supply at least one to proceed.")
    }
    
    if (!all(colnames(count) %in% rownames(meta))){
        stop("Sample names do not match between count table and meta table.")
    }
    meta <- meta[colnames(count),]
    if (!"batch" %in% colnames(meta)) {
        stop("The metadata table does not contain 'batch' column.")
    }
    batchid <- meta[['batch']]
    if (is.character(covariates) && covariates == "all") {
    meta_filtered <- meta[, !colnames(meta) %in% "batch", drop = FALSE]
    } else if (is.vector(covariates)) {
        # Ensure all covariates exist in meta
        if (!all(covariates %in% colnames(meta))) {
            stop("Not all covariates are included in the meta table.")
        } else {
            meta_filtered <- meta[, covariates, drop = FALSE]  # Select only the specified columns
        }
    } else {
        meta_filtered <- NULL  # Ensure meta_filtered is defined
    }


    if (!is.null(meta_filtered) && ncol(meta_filtered)>=3){
        alpha <- 0.01
    }else{
        alpha <- 0.1
    }

    batch_name = unique(batchid)

    # Create a list of sample IDs for each batch
    sample_list <- lapply(batch_name, function(x) {
    rownames(meta)[batchid == x]
    })

    # For each study, subset the count.genus matrix using the corresponding sample IDs
    O.list <- lapply(seq_along(batch_name), function(i) {
    as.matrix(count[, sample_list[[i]]])
    })
    
    # Create a list of meta data for each study based on sample_list
    meta.list <- lapply(seq_along(batch_name), function(i) {
        meta_filtered[sample_list[[i]], , drop = FALSE]  # Keep as data frame
    })

    r = sum(svd(as.matrix(count))$d>1e-3)

    if (verbose) {
        message("PASS")
    }
    output = list(count_list = O.list,
                  meta_list = meta.list,
                  dist_mat = dist_mat,
                  controls = list(alpha = alpha, beta = beta, neighbor = neighbor, sigma = sigma, 
                  r = r))

}

   