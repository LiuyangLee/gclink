#' @title Parse ORF Coordinates from Prodigal FASTA Headers
#'
#' @description Extracts ORF identifiers, start/end positions and strand
#'              orientation directly from the FASTA headers produced by
#'              Prodigal.  The resulting table is ready for downstream
#'              gene-cluster analyses.
#'
#' @param in_seq_data A data frame with two columns:
#'         \describe{
#'           \item{\code{SeqName}}{ORF Sequence Name.} \code{SeqName} includes the information of \code{qaccver},
#'           \code{start}, \code{end}, and \code{direction}.
#'           \item{\code{Sequence}}{ORF Sequence.}. e.g., "Kuafubacteriaceae--GCA_016703535.1---JADJBV010000001.1_1
#'           # 74 # 1018 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.685"
#'         }.
#' This data frame can be scaled from **Prodigal** FASTA output by the command:
#' seq_data = "Prodigal.fasta" %>%
#'    {Biostrings::readBStringSet(., format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)} %>%
#'    {data.frame(Sequence = .)} %>%
#'    {tibble::rownames_to_column(., "SeqName")}
#'
#'  **Prodigal** FASTA output: Each row must correspond to a header line (starting with \code{>}).
#'  The header is expected to follow the standard Prodigal format:
#'  \code{>contig_id # start # end # strand # ...}.
#' @return A data frame
#' @export
#' @examples
#' \dontrun{
#' loc <- orf_locate(data)
#' }
orf_locate = function(in_seq_data = seq_data
                      ){
  # obtain the positions of orf. this function will generete five columns
  # "SeqName", the name of "orf", and the position of orf including "start", "end", and "direction".
  # It is note that using awk would be faster:
  GC.location = in_seq_data %>%
    {dplyr::mutate(.,
                   qaccver = sapply(as.character(.$SeqName),
                                function(x) unlist(strsplit(x, " # "))[1]),
                   start = sapply(as.character(.$SeqName),
                                  function(x) unlist(strsplit(x, " # "))[2]),
                   end = sapply(as.character(.$SeqName),
                                function(x) unlist(strsplit(x, " # "))[3]),
                   direction = sapply(as.character(.$SeqName),
                                      function(x) unlist(strsplit(x, " # "))[4])
    )
    }
  return(GC.location)
}
