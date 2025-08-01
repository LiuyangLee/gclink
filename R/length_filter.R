#' @title Remove Length Outliers from BLAST Results
#'
#' @description Filters BLAST hits by removing ORFs whose gene (protein) length
#'              is an outlier within the corresponding gene group, as defined
#'              by the inter-quartile range (IQR).  Hits whose length falls
#'              outside the interval
#'              \code{[Q1 - down_IQR * IQR, Q3 + up_IQR * IQR]} are discarded.
#'
#' @param Data A data frame containing BLAST results.  Must include the columns
#'             \code{gene} (gene symbol) and \code{length} (ORF length in amino
#'             acids).
#' @param down_IQR Numeric multiplier applied to the IQR for the **lower**
#'                 bound (default: 1.5).
#' @param up_IQR Numeric multiplier applied to the IQR for the **upper** bound
#'               (default: 1.5).
#'
#' @return The input data frame with outlier rows removed.  The returned object
#'         is **ungrouped** regardless of the input grouping.
#' @export
#' @details
#'   - Filtering is performed **within each gene group**; outliers are determined
#'     independently for every gene symbol.
#'   - Progress messages report the number of rows before and after filtering.
#'   - Missing values in \code{length} are ignored when computing quantiles.
#'
#' @examples
#' \dontrun{
#' # Remove ORFs whose lengths are > 2 IQR beyond Q1 and Q3
#' clean <- length_filter(
#'   Data      = blast_hits,
#'   down_IQR  = 2,
#'   up_IQR    = 2
#' )
#' head(clean)
#' }
length_filter = function(Data = bin_genes,
                         down_IQR = 1.5,
                         up_IQR = 1.5
){
  # filter the blastp result based on the gene length.
  # the gene column is gerenated by orf_extract(), the filter threshold is that
  # A qaccver will be retained only when its length distribution is not an outliers,
  # e.g., length should distributed at least between the range of 1.5 * IQR
  # this would be done for all blast qaccvers within corresponding gene catogory.
  # length, the colnames for the gene length column
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"The qaccver number before filter is:",dim(Data)[1],'\n\n')
  Data = Data %>%
    dplyr::mutate(., length = as.numeric(.$length)) %>%
    dplyr::group_by(gene) %>%
    dplyr::filter(!
                    (
                      (length <
                         (stats::quantile(length, probs=c(.25), na.rm = T)
                          - down_IQR * stats::IQR(length, na.rm = T))
                      ) |
                        (length >
                           (stats::quantile(length, probs=c(.75), na.rm = T)
                            + up_IQR * stats::IQR(length, na.rm = T))
                        )
                    )
    ) %>%
    dplyr::ungroup()
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"The qaccver number after filter is:",dim(Data)[1],'\n\n')
  return(Data)
}
