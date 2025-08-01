#' @title Identify Breakpoints of Gene Clusters within a Contig
#'
#' @description Internal helper used by \code{\link{gc_cal}}.
#'              Given the ordered positions of *reference* genes on a contig,
#'              this function returns the genomic coordinates that mark the
#'              boundary of each candidate cluster.
#'              A boundary is declared whenever the gap between two successive
#'              reference genes exceeds the maximum spacing allowed by the
#'              cluster definition (\code{AllGeneNum - MinConSeq}).
#'
#' @param Data A numeric vector (ascending order) of ORF positions that carry
#'             one of the reference genes of interest.  Usually the vector
#'             \code{orf_position.tmp} created inside \code{gc_cal}.
#' @param AllGeneNum Integer.  Maximum genomic span (in ORF count) that the
#'                   algorithm is allowed to cover when defining a single
#'                   cluster.  Passed unchanged from \code{gc_cal}.
#' @param MinConSeq Integer.  Minimum number of consecutive reference genes
#'                  required to form a cluster.  Passed unchanged from
#'                  \code{gc_cal}.
#'
#' @return A numeric vector containing every position that marks the **end**
#'         of a putative gene cluster.  These values are subsequently used as
#'         breakpoints to slice the contig into candidate regions in the
#'         downstream functions \code{gc_position()} and \code{gc_range()}.
#' @export
#' @details
#'   - If the gap between two consecutive reference genes is larger than
#'     \code{AllGeneNum - MinConSeq}, the left-hand gene is recorded as the
#'     **last member** of the preceding cluster.
#'   - The final reference gene is always appended to the vector as the last
#'     boundary, ensuring the rightmost cluster is not lost.
#'
#' @examples
#' \dontrun{
#' # Example: reference genes at positions 5, 7, 9, 25, 27, 29
#' breakpoints <- gc_cluster(
#'   Data       = c(5, 7, 9, 25, 27, 29),
#'   AllGeneNum = 20,
#'   MinConSeq  = 10
#' )
#' cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),breakpoints)
#' }
gc_cluster =  function(Data = orf_position.tmp,
                       AllGeneNum = 30,
                       MinConSeq = 15){
  # obtain the sites of potential GCs within one contig
  # MinConSeq = floor(as.numeric(AllGeneNum) * as.numeric(Proportion)) # ceiling
  GeneNum=length(Data)
  if (GeneNum >= MinConSeq){
    cluster.tmp = vector()
    for (site in Data) {
      ordination = which(Data %in% site)
      # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('orf site: ', site))
      # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('target gene ordination: ', ordination))
      if (ordination < GeneNum){
        eachdif.tmp = Data %>%
          {abs((.[ordination+1]) - (.[ordination]))}
        # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('eachdif.tmp: ', eachdif.tmp))
        # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('[Difference between up and down is ', eachdif.tmp))
        # the max site-difference between two orfs is AllGeneNum-MinConSeq.
        # 两两orf之间相差至多AllGeneNum-MinConSeq个位置
        if (eachdif.tmp > (AllGeneNum-MinConSeq)){
          cluster.tmp = c(cluster.tmp, site)
          # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('cluster.tmp: ', cluster.tmp))
        }
      } else {
        cluster.tmp = c(cluster.tmp, site)
        # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('cluster.tmp: ', cluster.tmp))
      }
    }
  }
  return(cluster.tmp)
}
