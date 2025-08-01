#' @title Scale Gene-Cluster Coordinates for Visualization
#'
#' @description Prepares a gene‐cluster annotation table for downstream plotting by
#'              converting absolute genomic coordinates into relative positions,
#'              ensuring that every cluster starts at 0 and is oriented consistently.
#'              Hypothetical ORFs (originally labeled with \code{NA}) and missing
#'              labels are replaced with placeholders, and factor levels are set
#'              as requested.
#'
#' @param GC_meta A data frame containing gene-cluster information.  Must include
#'                the columns \code{gene_cluster}, \code{gene_group}, \code{gene_label},
#'                \code{start}, \code{end}, and \code{direction} (numeric: \code{1}
#'                for forward, \code{-1} for reverse).
#' @param levels_gene_group Character vector specifying the desired factor levels
#'                          for \code{gene_group}.  Group names should
#'                          appear in the order required for plotting legends.
#'
#' @return The input data frame with the following **new or overwritten** columns:
#'         \describe{
#'           \item{\code{gene_label}}{Empty string (\code{""}) if originally \code{NA}.}
#'           \item{\code{gene_group}}{Set to \code{"hypothetical ORF"} if originally
#'                                    \code{NA}, then coerced to a factor using
#'                                    \code{levels_gene_group}.}
#'           \item{\code{Pgenome}}{Factor version of \code{gene_cluster}; levels
#'                                 follow the order of appearance in the data.}
#'           \item{\code{Pstart}, \code{Pend}}{Relative start and end coordinates
#'                                              (numeric) within each cluster,
#'                                              scaled so that the left-most gene
#'                                              starts at 0.}
#'           \item{\code{Pdirection}}{Logical vector: \code{TRUE} for forward,
#'                                     \code{FALSE} for reverse.}
#'         }
#' @export
#' @details
#'   - Absolute \code{start}/\code{end} values are **not** modified; scaled
#'     values are stored in new columns (\code{Pstart}, \code{Pend}).
#'   - \code{Pgenome} can be swapped for any unique identifier (e.g., \code{Genome})
#'     downstream if each genome contains only one cluster.
#'
#' @examples
#' \dontrun{
#' # Prepare levels so "hypothetical ORF" appears last in legends
#' scaled <- gc_scale(
#'   GC_meta          = GC_annotation,
#'   levels_gene_group = c('bch','puh','puf','crt', 'acsF', 'assembly','regulator','hypothetical ORF')
#' )
#' head(scaled)
#' }
gc_scale = function(GC_meta = GC_meta,
                    levels_gene_group = levels_gene_group
){
  # This function is aimed to scale the data for further arrangement plot.
  # The input of "levels_gene_group" should include the group of 'hypothetical ORF'.
  # Meanwhile, the "gene_group" column with NA will be substituted when 'hypothetical ORF'.
  # The "gene_label" column with NA will be substituted when ''.
  # Columns "Pgenome", "gene_group", "Pstart", "Pend", will generate and scaled.
  # Column "Pgenome" can be changed to other label when they are unique to plot,
  # e.g., when one "Genome" only have one cluster, you can specify the "Genome" as the "Pgenome".
  GC_meta[is.na(GC_meta$gene_label),'gene_label']=''
  GC_meta[is.na(GC_meta$gene_group),'gene_group']='hypothetical ORF'
  # specify the unit for plot: Pgenome
  GC_meta$Pgenome=factor(GC_meta$gene_cluster,
                         levels = unique(GC_meta$gene_cluster))
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"Find Gene clusters: ",'\n\n', paste0(unique(GC_meta$gene_cluster),'\n'))
  GC_meta$gene_group=factor(GC_meta$gene_group,
                            levels = levels_gene_group)
  GC_meta$Pstart = as.numeric(GC_meta$start)
  GC_meta$Pend = as.numeric(GC_meta$end)
  for (variable in unique(GC_meta$gene_cluster)) {
    Start = (min(GC_meta[GC_meta$gene_cluster %in% variable, 'Pstart']))
    GC_meta[GC_meta$gene_cluster %in% variable, 'Pstart']=GC_meta[GC_meta$gene_cluster %in% variable, 'Pstart'] - Start
    GC_meta[GC_meta$gene_cluster %in% variable, 'Pend']=GC_meta[GC_meta$gene_cluster %in% variable, 'Pend'] - Start
  }
  GC_meta$Pdirection = GC_meta$direction
  GC_meta$Pdirection[(GC_meta$Pdirection==1)] <- TRUE
  GC_meta$Pdirection[(GC_meta$Pdirection==-1)] <- FALSE
  return(GC_meta)
}
