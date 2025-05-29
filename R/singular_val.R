#' @title Generate Singular Value Plots for Each Dataset.
#' 
#' @import ggplot2
#' @import ggpubr
#'
#' @description This function produces singular value plots for each input dataset to assess
#' the validity of the low-rank assumption. A rapid decay in the singular values
#' indicates that the dataset can be effectively approximated by matrix factorization.
#' 
#' @param count The integrated count table of taxa by samples. The
#' \code{count} parameter should be provided as either a \code{matrix} or a
#' \code{data.frame}.
#' @param meta The integrated meta table \code{meta} contains a column named ``batch'' 
#' which stores batch id.
#'  
#' @return A list of ggplot objects displaying the singular values for each dataset.
#' 
#' @export 
#' 
plot_singular_values <- function(count, meta){
    if (!all(colnames(count) %in% rownames(meta))){
        stop("Sample names do not match between count table and meta table.")
    }
    meta <- meta[colnames(count),]
    if (!"batch" %in% colnames(meta)) {
        stop("The metadata table does not contain 'batch' column.")
    }
    batchid <- meta[['batch']]
    batch_name = unique(batchid)

    # Create a list of sample IDs for each batch
    sample_list <- lapply(batch_name, function(x) {
    rownames(meta)[batchid == x]
    })

    # For each study, subset the count.genus matrix using the corresponding sample IDs
    O.list <- lapply(seq_along(batch_name), function(i) {
    as.matrix(count[, sample_list[[i]]])
    })

    svd_results <- list()

    # Loop through each dataset in O.list
    for (i in seq_along(O.list)) {
    svd_d <- svd(t(O.list[[i]]))$d  # Get singular values
    
    # Create a temporary data frame for the current dataset
    temp_df <- data.frame(
        value = svd_d,
        taxon = seq_along(svd_d),  # Correct indexing
        dataset = batch_name[i]
    )
    
    # Store in list
    svd_results[[i]] <- temp_df
    }
    # Combine all results into a single data frame
    svd_value <- do.call(rbind, svd_results)

    # plot
    plot_list <- lapply(svd_results, function(df) {
        ggplot(df, aes(x = taxon, y = value)) +
        geom_line(size = 1) +
        facet_wrap(~dataset, ncol = length(svd_results)) +
        theme_bw() +
        theme(
            plot.title  = element_blank(),
            axis.title  = element_blank(),
            axis.text   = element_text(size = 10),
            legend.text = element_text(size = 10)
        )
    })
    
    # Combine plots into one figure
    combined_plot <- ggarrange(plotlist = plot_list, nrow = 1)

    # Add a title
    final_plot <- annotate_figure(combined_plot,
                                    top = text_grob("Plot of Singular Values", face = "bold", size = 14))

  return(final_plot)
}
