#' @title Identify and Extract Gene Clusters from Scaled BLAST Data
#'
#' @description This function screens contigs for regions that contain a
#'              pre-defined set of “reference” genes (e.g., photosynthetic genes, viral genes)
#'              arranged in a continuous block.  Contigs are
#'              first coarsely filtered by the minimum number of reference genes
#'              they carry, then finely scanned for clusters that satisfy user-
#'              defined density and contiguity criteria.  Each detected cluster
#'              is returned with a unique \code{gene_cluster} identifier.
#'
#' @param Data A data frame produced by \code{\link{orf_extract}} (i.e., a scaled
#'             BLAST table).  Must include the columns \code{genome_contig},
#'             \code{gene}, and \code{orf_position}.
#' @param in_gene_list A character vector of “reference” gene symbols (e.g.,
#'                   \code{photosynthesis_gene_list}) that are expected
#'                   to appear in the target cluster(s).
#' @param AllGeneNum Integer.  Maximum total ORF count (annotated plus hypothetical)
#'                   that the algorithm is allowed to span when defining a cluster
#'                   (default: 30).
#' @param MinConSeq Integer.  Minimum number of **reference genes** that must be
#'                  present **and consecutive** within the candidate cluster
#'                  (default: 15).  Must satisfy \code{1 <= MinConSeq <= AllGeneNum}.
#'
#' @return A data frame identical in structure to \code{Data} but filtered to
#'         contain only those rows that belong to valid clusters.  An extra
#'         column \code{gene_cluster} (format: \code{genome_contig---N}) is added
#'         to uniquely label every cluster.
#' @export
#' @details
#'   1. **Coarse filter**:  Contigs with fewer than \code{MinConSeq} reference
#'      genes are discarded.
#'   2. **Fine scan**:  For each remaining contig, the algorithm slides a
#'      window that can encompass up to \code{AllGeneNum} consecutive ORFs
#'      and retains windows that contain at least \code{MinConSeq} reference
#'      genes in uninterrupted order.
#'   3. **Cluster labelling**:  Each valid cluster receives a unique ID
#'      (\code{genome_contig---1}, \code{genome_contig---2}, …).
#'
#' @examples
#' \dontrun{
#' # Detect photosynthetic gene clusters (photosynthesis_gene_list) up to 30 ORFs wide
#' # requiring at least 15 reference genes to be consecutive
#' clusters <- gc_cal(
#'   Data       = scaled_blast,
#'   in_gene_list = photosynthesis_gene_list,
#'   AllGeneNum = 30,
#'   MinConSeq  = 15
#' )
#' head(clusters)
#' }
gc_cal = function(Data = bin_genes,
                  in_gene_list = photosynthesis_gene_list,
                  AllGeneNum = 30,
                  MinConSeq = 15){
  # Data is the scaled blast table.
  # gene_lists is defined as a vectore include any genes involved in the "gene" column in the scaled blast table.
  # default is photosynthesis_gene_list,
  # AllGeneNum is the expected gene number in the gene cluster.
  # The AllGeneNum is generally set larger than the number of "reference genes" to include unknown (hypothetical) genes in the cluster.
  # The Proportion is an expected minimum proportion that "reference genes" should account for in the gene cluster.
  # MinConSeq: the min number of target gene in the gene cluster
  # MinConSeq = floor(as.numeric(AllGeneNum) * as.numeric(Proportion)) # ceiling
  if (MinConSeq < 1) {
    stop("Be sure that threshold number of connective sequences (MinConSeq) must >= 1 !!!")
  }
  if (MinConSeq > AllGeneNum) {
    stop("Be sure that threshold number of connective sequences (MinConSeq) must <= AllGeneNum !!!")
  }
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"The gene lists include: ")
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),in_gene_list)
  result=data.frame()
  # Firstly roughly filter
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0("The qaccver number before roughly filter is: ", dim(Data)[1]))
  Data = Data %>%
    {.[.$genome_contig %in% {{.} %>%
        {stats::aggregate(
          x=list(length = (.$qaccver)), #处理对象
          by=list(genome = (.$genome),
                  genome_contig = (.$genome_contig)
          ), #条件
          FUN=length
        )} %>%
        {.[.$length>=MinConSeq,]}
    }$genome_contig,]} %>%
    {.[.$gene %in% in_gene_list, ]}
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0("The qaccver number after roughly filter is: ", dim(Data)[1]))
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"Starting precisely filter!")
  # Then precisely filter
  for (variable in unique(Data$genome_contig)) {
    # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('Genome-Contig ', variable, ' is running!'))
    Data.tmp = Data[Data$genome_contig %in% variable, ]
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('orf_position: ', Data.tmp$orf_position))
    order.tmp = base::order(Data.tmp$orf_position, decreasing = F)
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('order.tmp: ', order.tmp))
    Data.order.tmp = Data.tmp[order.tmp, ]
    orf_position.tmp = Data.order.tmp$orf_position
    # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('orf_position.tmp: ', Data.order.tmp$orf_position))
    if (length(orf_position.tmp) >= MinConSeq){
      # obtain the sites of potential GCs within one contig
      cluster.tmp = orf_position.tmp %>%
        {gc_cluster(Data = .,
                    AllGeneNum = AllGeneNum,
                    MinConSeq = MinConSeq)}
      num = 1 # GC计数器
      for (eachcluster in 1:(length(cluster.tmp))) {
        # 分割基因cluster
        Norf_position = orf_position.tmp %>%
          {gc_position(Data = .,
                       Cluster = cluster.tmp,
                       EachCluster = eachcluster)}
        # 每隔MinConSeq个基因，orf位置差异<=AllGeneNum
        if (length(Norf_position) >= MinConSeq){
          GC.site = Norf_position %>%
            {gc_range(Norf_position = .,
                      AllGeneNum = AllGeneNum,
                      MinConSeq = MinConSeq)}
          if (length(GC.site) >= MinConSeq){
            result.tmp = Data.order.tmp %>%
              {.[.$orf_position %in% GC.site, ]}
            # add the gene_cluster column. This is useful especially when one genome or contig has multiple clusters.
            result.tmp$gene_cluster = result.tmp %>%
              {paste0((.$genome_contig), '---', num)}
            result = rbind(result, result.tmp)
            if (length(result.tmp$gene_cluster)>0){
              cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('Length of gene cluster *', unique(result.tmp$gene_cluster), '* is ', length(GC.site), '!'))
              # cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('[ GC site is ', (GC.site), ' ]'))
            }
            num = num + 1
          }
        }
      }
    }
  }
  if (length(result$gene_cluster) >= 1){
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0("The qaccver number after precisely filter is: ", dim(result)[1]))
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"The gene clusters have been successfully filtered~~~")
  } else {
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"No gene clusters can be found based on current settings.")
  }
  result = result %>%
    {.[base::order(.$gene_cluster, .$orf_position), ]}
  return(result)
}
