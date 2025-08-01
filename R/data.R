#' BLASTP Results for test Proteins
#'
#' A dataset containing BLASTP alignment results for proteins from test genomes of bacteria,
#' including alignment metrics.
#'
#' @format A data frame with multiple rows and 12 variables:
#' \describe{
#'   \item{qaccver}{Character. Protein query identifier.}
#'   \item{saccver}{Character. Subject protein identifier from reference databases.}
#'   \item{pident}{Numeric. Percentage identity between query and subject sequences.}
#'   \item{length}{Numeric. Length of the aligned sequence region.}
#'   \item{mismatch}{Numeric. Number of mismatches in the alignment.}
#'   \item{gapopen}{Numeric. Number of gap openings in the alignment.}
#'   \item{qstart}{Numeric. Start position in query sequence.}
#'   \item{qend}{Numeric. End position in query sequence.}
#'   \item{sstart}{Numeric. Start position in subject sequence.}
#'   \item{send}{Numeric. End position in subject sequence.}
#'   \item{evalue}{Numeric. Expect value for the alignment.}
#'   \item{bitscore}{Numeric. Bit score for the alignment.}
#' }
#' @export
"blastp_df"

#' Photosynthesis Gene Classification Groups
#'
#' A dataset mapping photosynthesis-related genes to their functional groups and subunits.
#'
#' @format A data frame with multiple rows and 3 variables:
#' \describe{
#'   \item{gene}{Character. Gene identifier (e.g., "bchB").}
#'   \item{gene_group}{Character. Functional group (e.g., "bch").}
#'   \item{gene_label}{Character. Subunit designation (e.g., "B").}
#' }
#' @export
"PGC_group"

#' Photosynthesis Gene List
#'
#' A comprehensive list of genes involved in photosynthetic processes.
#'
#' @format A character vector containing gene identifiers:
#' \describe{
#'   Each element represents a photosynthesis-related gene symbol (e.g., "acsF", "bchB").
#' }
#' @export
"photosynthesis_gene_list"

#' Genomic Sequence Data with Annotations
#'
#' A dataset containing DNA sequences from test bacteria with detailed annotation metadata.
#' The first column combines multiple annotation elements separated by semicolons.
#'
#' @format A data frame with multiple rows and 2 variables:
#' \describe{
#'   \item{SeqName}{Character. Combined annotation fields separated by semicolons, containing:
#'     \itemize{
#'       \item \code{ID}: Sequence identifier (e.g., "1_7")
#'       \item \code{partial}: Completion status ("00" for complete, "01" for partial)
#'       \item \code{start_type}: Translation initiation codon (e.g., "GTG", "ATG")
#'       \item \code{rbs_motif}: Ribosome binding site motif (e.g., "GGAG/GAGG")
#'       \item \code{rbs_spacer}: RBS spacer length (e.g., "5-10bp")
#'       \item \code{gc_cont}: GC content (e.g., "0.673")
#'     }
#'   }
#'   \item{Sequence}{Character. DNA sequence (when available) in FASTA format}
#' }
#' @export
"seq_data"

.reference_gc_data <- function() {
  list(
    ref1 = blastp_df,
    ref2 = PGC_group,
    ref3 = photosynthesis_gene_list,
    ref4 = seq_data
  )
  invisible(NULL)
}
