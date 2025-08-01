#' @title Parse ORF Coordinates from Prodigal FASTA Headers
#'
#' @description Extracts ORF identifiers, start/end positions and strand
#'              orientation directly from the FASTA headers produced by
#'              Prodigal.  The resulting table is ready for downstream
#'              gene-cluster analyses.
#'
#' @param in_seq_data     A data frame with two columns:
#'   \describe{
#'     \item{`SeqName`}{ORF identifier (Prodigal format: `>ORF_id # start # end # strand # ...`).}
#'     \item{`Sequence`}{ORF sequence.}
#'   }
#'   Example:
#'   `"Kuafubacteriaceae--GCA_016703535.1---JADJBV010000001.1_1 # 74 # 1018 # 1 # ..."`
#'   Can be imported from **Prodigal** FASTA using:
#'   ```r
#'   seq_data <- Biostrings::readBStringSet("Prodigal.fasta",format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE) %>%
#'     data.frame(Sequence = .) %>%
#'     tibble::rownames_to_column("SeqName")
#'   ```
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
