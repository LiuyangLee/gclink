#' @title Determine ORF Range for a Candidate Gene Cluster
#'
#' @description Internal helper used by \code{\link{gc_cal}}.
#'              After \code{\link{gc_position}} has isolated the ORF positions
#'              belonging to a single cluster, this function **validates** and
#'              **trims** that range so that the final span (distance between the
#'              first and last retained ORF) does not exceed \code{AllGeneNum}.
#'              The goal is to retain the **largest contiguous block** that still
#'              satisfies the user-defined size limit.
#'
#' @param Norf_position Numeric vector of ORF positions (ascending) that belong
#'                      to the current candidate cluster (output from
#'                      \code{\link{gc_position}}).
#' @param AllGeneNum Integer.  Maximum allowed genomic span (in ORF count) for
#'                   the final cluster.
#' @param MinConSeq Integer.  Minimum number of consecutive reference genes
#'                  required for the cluster.
#'
#' @return A numeric vector containing the **final ORF positions** that define
#'         the validated gene cluster.  If no valid block can be produced, the
#'         vector will be empty.
#' @export
#' @details
#'   - For every reference gene in \code{Norf_position}, the function evaluates
#'     whether a window of at least \code{MinConSeq} consecutive reference genes
#'     centred on that gene can fit within \code{AllGeneNum} consecutive ORFs.
#'   - Genes that pass the test are collected in \code{retain.site}.
#'   - The **minimal** and **maximal** positions in \code{retain.site} are then
#'     used to slice the full ORF range, guaranteeing that the final cluster
#'     length ≤ \code{AllGeneNum}.
#'
#' @examples
#' \dontrun{
#' gc_range(
#'   Norf_position = c(10, 12, 14, 16, 18, 30, 60),
#'   AllGeneNum    = 20,
#'   MinConSeq     = 5
#' )
#' }
gc_range =  function(Norf_position = Norf_position,
                     AllGeneNum = 20,
                     MinConSeq = 10){
  #MinConSeq = floor(as.numeric(AllGeneNum) * as.numeric(Proportion)) # ceiling
  retain.site = vector()
  GC.site = vector()
  for (Nsite in Norf_position) {
    Nordination = which(Norf_position %in% Nsite)
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('Nsite: ', Nsite),'\n\n')
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('Ntarget gene ordination: ', Nordination),'\n\n')
    if ((Nordination < MinConSeq) & ((length(Norf_position)-Nordination) < MinConSeq)){
      clusterdif.tmp = Norf_position %>%
        {max(c(abs((.[1]) - (.[Nordination])),
               abs((.[length(Norf_position)]) - (.[Nordination]))))}
    } else if (Nordination < MinConSeq){
      clusterdif.tmp = Norf_position %>%
        {max(c(abs((.[1]) - (.[Nordination])),
               abs((.[Nordination+(MinConSeq-1)]) - (.[Nordination]))))}
    } else if ((length(Norf_position)-Nordination) < MinConSeq){
      clusterdif.tmp = Norf_position %>%
        {max(c(abs((.[Nordination-(MinConSeq-1)]) - (.[Nordination])),
               abs((.[length(Norf_position)]) - (.[Nordination]))))}
    } else {
      clusterdif.tmp = Norf_position %>%
        {min(c(abs((.[Nordination-(MinConSeq-1)]) - (.[Nordination])),
               abs((.[Nordination+(MinConSeq-1)]) - (.[Nordination]))))}
    }
    if (clusterdif.tmp <= AllGeneNum){
      #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0(' [ Retain this site because Max difference between start/end and here is ', clusterdif.tmp, ' ]'),'\n\n')
      retain.site = c(retain.site, Nsite)
    } else {
      #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0(' [ Discard this site because Max difference between start/end and here is ', clusterdif.tmp, ' ]'),'\n\n')
    }
  }
  if (length(retain.site) >= 2){
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0(' [ min (retain.site) is ', min(retain.site),' ]'),'\n\n')
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0(' [ max (retain.site) is ', max(retain.site),' ]'),'\n\n')
    GC.site = Norf_position  %>% {.[. %in% c(min(retain.site):max(retain.site))]}
  } else if (length(retain.site) == 1){
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0(' [ the (retain.site) is ', retain.site,' ]'),'\n\n')
    GC.site = Norf_position  %>% {.[. %in% c(min(retain.site):max(retain.site))]}
  }
  return(GC.site)
}
