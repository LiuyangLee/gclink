#' @title Extract ORF Positions for One Specific Gene Cluster
#'
#' @description Internal helper used by \code{\link{gc_cal}}.
#'              Given the ordered positions of **all** reference genes on a contig
#'              (\code{Data}) and the vector of cluster breakpoints produced by
#'              \code{\link{gc_cluster}} (\code{Cluster}), this function slices
#'              \code{Data} to return the positions that belong to the
#'              \code{EachCluster}-th cluster.
#'
#' @param Data A numeric vector (ascending order) of ORF positions that carry
#'             one of the reference genes of interest.  Usually the vector
#'             \code{orf_position.tmp} created inside \code{gc_cal}.
#' @param Cluster Numeric vector returned by \code{\link{gc_cluster}}; each
#'                element marks the **last position** of a candidate cluster.
#' @param EachCluster Integer index specifying which cluster to extract
#'                    (1 = first cluster, 2 = second cluster, …).
#' @export
#' @return A numeric vector containing the ORF positions that constitute the
#'         requested gene cluster.
#' @examples
#' \dontrun{
#' gc_position(Data = orf_position.tmp,
#'             Cluster = cluster.tmp,
#'             EachCluster = eachcluster)
#' }
gc_position =  function(Data = orf_position.tmp,
                        Cluster = cluster.tmp,
                        EachCluster = eachcluster){
  # obtain the new position of each GC within a contig
  # 分割基因cluster
  # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('EachCluster: ', EachCluster))
  breakpoint = which(Data %in% Cluster[EachCluster])
  if (EachCluster == 1){
    # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('EachCluster: ', EachCluster))
    Norf_position = Data[1:breakpoint]
    # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('Norf_position: ', Norf_position))
  } else {
    Norf_position = Data[(which(Data==Cluster[EachCluster-1])+1):breakpoint]
    # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('Norf_position: ', Norf_position))
  }
  return(Norf_position)
}
