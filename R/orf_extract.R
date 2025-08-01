#' @title Extract ORF and Genome Information from BLAST or BLASTP Results
#'
#' @description This function parses BLASTP result tables to extract structured genome,
#'              contig, ORF, and gene information from the query and subject identifiers.
#'              It is designed for downstream analyses requiring explicit separation of
#'              genome, contig, and ORF identifiers from concatenated BLAST headers.
#' @param bin_genes A data frame containing BLASTP results with at least 2 standard columns:
#'                  \code{qaccver}, \code{saccver}.
#' the column of qaccver should include both of the genome name and predicted contig name, which is concatenated by a separator "---".
#' for example, for the qaccver "p__Myxococcota--c__Kuafubacteria--o__Kuafubacteriales--f__Kuafubacteriaceae--GCA_016703535.1---JADJBV010000001.1_150",
#' the genome name is "p__Myxococcota--c__Kuafubacteria--o__Kuafubacteriales--f__Kuafubacteriaceae--GCA_016703535.1",
#' the contig name is "JADJBV010000001.1", the orf name is "JADJBV010000001.1_150", and the orf_position is "150".
#' the column of saccver must include the gene name and may include the gene information, which are concatenated by a separator "_".
#' for example, for the saccver "bchC_Methyloversatilis_sp_RAC08_BSY238_2447_METR",
#' the gene name is "bchC", the gene information is Methyloversatilis_sp_RAC08_BSY238_2447_METR that can help you understand the source of gene.
#' @return The original data frame with six additional columns:
#'         \describe{
#'           \item{\code{genome}}{Genome identifier extracted from \code{qaccver}.}
#'           \item{\code{contig}}{Contig identifier extracted from \code{qaccver}.}
#'           \item{\code{orf}}{Full ORF identifier extracted from \code{qaccver}.}
#'           \item{\code{genome_contig}}{Concatenated genome and contig IDs (genome---contig).}
#'           \item{\code{gene}}{Gene symbol extracted from \code{saccver}.}
#'           \item{\code{orf_position}}{Numeric ORF position extracted from the ORF identifier.}
#'         }
#' @export
#' @examples
#' \dontrun{
#' # Example BLASTP input
#' data(blastp_df)
#' parsed_df <- orf_extract(blastp_df)
#' head(parsed_df)
#' }
orf_extract = function(bin_genes = blastp_df
){
  bin_genes$genome = sapply(as.character(bin_genes$qaccver),
                            function(x) unlist(strsplit(x, "---"))[1])
  bin_genes$orf = sapply(as.character(bin_genes$qaccver),
                         function(x) utils::tail(unlist(strsplit(x, "---")), 1))
  bin_genes$contig = sapply(as.character(bin_genes$orf),
                            function(x) unlist(strsplit(x, "_\\d*$"))[1])
  bin_genes$genome_contig = bin_genes %>%
    {paste0(.$genome, "---", .$contig)}
  bin_genes$gene = sapply(as.character(bin_genes$saccver),
                          function(x) unlist(strsplit(x, "_"))[1])
  bin_genes$orf_position = as.numeric(sapply(as.character(bin_genes$orf),
                                             function(x) (utils::tail((unlist(strsplit(x, "_"))),n=1))))
  return(bin_genes)
}
