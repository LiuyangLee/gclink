#' @title Complete Gene Clusters by Adding Missing ORFs
#'
#' @description
#' Expands gene cluster tables to include **all ORFs** (annotated and hypothetical) within contigs,
#' normalizing cluster representations for downstream analysis and plotting. Ensures consistent
#' ORF spacing/length across clusters by inserting missing rows.
#'
#' @param Data A `data.frame` of annotated ORFs with required columns:
#'   \itemize{
#'     \item `qaccver`: ORF identifier (genome---contig_orf_position format).
#'     \item `genome`, `contig`: Genome and contig names.
#'     \item `gene`: Gene symbol (NA for hypothetical ORFs).
#'     \item `orf_position`: Absolute ORF position on the contig.
#'     \item `gene_cluster`: Cluster identifier.
#'   }
#' @param Annotation A `data.frame` of full ORF annotations (e.g., from \code{\link{orf_extract}}).
#'   Must include `qaccver` and `orf_position`.
#' @param orf_before_first Integer. Hypothetical ORFs to insert **before** the first annotated ORF
#'                         in each cluster (bounded by contig start). Default: `0`.
#' @param orf_after_last Integer. Hypothetical ORFs to append **after** the last annotated ORF
#'                       (bounded by contig end). Default: `0`.
#' @param orf_range Character. Controls ORF inclusion and annotation merging:
#'   \itemize{
#'     \item `"All"`: Include every ORF in the contig range and merge all annotations (default).
#'     \item `"OnlyAnnotated"`: Keep only ORFs present in `Annotation` and merge their annotations.
#'     \item `"IgnoreAnnotated"`: Include all ORFs but skip merging with `Annotation`.
#'   }
#'
#' @return A `data.frame` with one row per ORF (real/hypothetical), sorted by `gene_cluster` and
#'         `orf_position`. Added columns:
#'   \describe{
#'     \item{`GC_orf_position`}{Relative position within cluster (1-indexed).}
#'     \item{`GC_present_length`}{Count of annotated ORFs in the cluster.}
#'     \item{`GC_absent_length`}{Count of inserted hypothetical ORFs.}
#'     \item{`GC_length`}{Total ORFs (`GC_present_length + GC_absent_length`).}
#'   }
#'
#' @details
#' * **Hypothetical ORFs** are inserted as rows with `gene = NA`.
#' * Output is always sorted by `gene_cluster` and `orf_position`.
#' * Progress messages are printed to console with timestamps.
#' * Contig bounds are respected—insertions never exceed actual ORF positions in `Annotation`.
#'
#' @export
#' @examples
#' \dontrun{
#' # Expand clusters with 2 ORFs upstream and downstream of each annotated block
#' sbgc_add <- gc_add(
#'   Data        = sbgc,
#'   Annotation  = bin_genes,
#'   orf_before_first = 2,
#'   orf_after_last   = 2,
#'   orf_range       = "All"
#' )
#' head(sbgc_add)
#' }
gc_add = function(Data = sbgc, Annotation = bin_genes,
                  orf_before_first = 0, orf_after_last = 0, orf_range = "All"){
  # This script is used to scale the data for further plotting the arrangement of the gene cluster.
  # the function include adding the information of orf_position for hypothetical or unknown genes.
  # and adding the numbers of reference gene and hypothetical (unknown) gene to the final table.
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"Starting add rows for unannotated ORFs!",'\n\n')
  Data = dplyr::select(Data, c("qaccver","genome","orf","contig",
                        "genome_contig","gene", "orf_position", "gene_cluster"))

  result = data.frame()
  for (variable in unique(Data$gene_cluster)) {
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0('Gene cluster ', variable, ' is running!'))
    result.GC = Data %>%
      {.[.$gene_cluster %in% variable, ]}
    GC.site = result.GC$orf_position

    # get the max site based on the orf site information from Annotation file
    Annotation.GC = Annotation %>%
      {.[.$contig %in% unique(result.GC$contig), ]}
    max.orf_site.Annotation.GC = max(Annotation.GC$orf_position)

    start_site = max((min(GC.site)-orf_before_first), 1)
    end_site = max(GC.site)+orf_after_last

    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0("start_site: ", start_site),'\n\n')
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),paste0("end_site: ", end_site),'\n\n')

    if (end_site > max.orf_site.Annotation.GC) {
      end_site = max.orf_site.Annotation.GC
    }

    all_sites = c(start_site:end_site)
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'), paste0("Gene cluster sites of ",variable, ':\n\n'), all_sites,'\n\n')
    absent_site = all_sites[!(all_sites %in% result.GC$orf_position)]

    result.GC$GC_orf_position = GC.site - min(all_sites) + 1

    result.unGC = data.frame()
    for (site in absent_site){
      result.unGC.tmp = data.frame(
        qaccver = result.GC %>%
          {paste0(unique(.$genome),'---', unique(.$contig),'_', site)},
        genome = unique(result.GC$genome),
        orf = paste0(unique(result.GC$contig), '_', site),
        contig = unique(result.GC$contig),
        genome_contig = unique(result.GC$genome_contig),
        orf_position = site,
        gene_cluster = unique(result.GC$gene_cluster),
        GC_orf_position = site - min(all_sites) + 1
      )
      result.unGC = dplyr::bind_rows(result.unGC, result.unGC.tmp)
    }
    result.tmp = dplyr::bind_rows(result.GC, result.unGC) %>%
      {.[order(.$GC_orf_position, decreasing = F), ]}
    result.tmp$GC_present_length = length(GC.site)
    result.tmp$GC_absent_length = length(absent_site)
    result = dplyr::bind_rows(result, result.tmp)
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"result: ")
    #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),result[1:16,])
  }
  if (length(result$gene_cluster) >= 1){
    result$GC_orf_position = as.numeric(result$GC_orf_position)
    result$orf_position = as.numeric(as.character(result$orf_position))
    #result = result %>%
    #  {dplyr::arrange(.,gene_cluster, orf_position)}
    result = result %>% {.[order(.$gene_cluster, .$orf_position), ]}
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"Rows for absent genes have be added~~~",'\n\n')
  } else {
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"No gene clusters can be found based on current settings, no rows can be added.",'\n\n')
  }
  result$GC_length = result$GC_present_length + result$GC_absent_length
  Annotation_to_merge = dplyr::select(Annotation, -c("genome","orf","contig",
                                              "genome_contig","gene", "orf_position"))
  if (orf_range == "All"){
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"All ORFs are used~~~",'\n\n')
    output = result %>%
      #{.[.$qaccver %in% unique(Annotation_to_merge$qaccver), ]} %>%
      merge(Annotation_to_merge, ., all.y = T, by = 'qaccver') %>%
      {.[order(.$gene_cluster, .$orf_position), ]}
  } else if (orf_range == "OnlyAnnotated"){
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"Only Annotated ORFs are used~~~",'\n\n')
    output = result %>%
      {.[.$qaccver %in% unique(Annotation_to_merge$qaccver), ]} %>%
      merge(Annotation_to_merge, ., all.y = T, by = 'qaccver') %>%
      {.[order(.$gene_cluster, .$orf_position), ]}
  } else if (orf_range == "IgnoreAnnotated"){
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),"Warning: Annotation informantion of ORFs are not used!",'\n\n')
    output = result %>%
      {.[order(.$gene_cluster, .$orf_position), ]}
  }
  return(output)
}
