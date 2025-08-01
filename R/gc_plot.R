#' @title Plot Scaled Gene Clusters with Arrows
#'
#' @description Generates a compact, publication-ready arrow plot of gene
#'              clusters previously processed with \code{\link{gc_scale}}.
#'              Each cluster is displayed on its own horizontal track, with
#'              gene directionality, labels, and user-defined colour schemes.
#'
#' @param data A data frame produced by \code{\link{gc_scale}}, must include
#'             the columns \code{Pstart}, \code{Pend}, \code{Pgenome},
#'             \code{gene_group}, \code{Pdirection}, and
#'             \code{gene_label}.
#' @param color_theme Character vector of colours, **in the same order** as the
#'                    factor levels of \code{gene_group}.  Defaults to an
#'                    8-colour palette; supply fewer or more colours as needed.
#' @export
#' @return A \code{ggplot} object that can be further customised or printed.
#'         The plot uses \code{gggenes::geom_gene_arrow()} and
#'         \code{gggenes::geom_gene_label()} with the following elements:
#'         \itemize{
#'           \item One facet per \code{Pgenome} (cluster).
#'           \item Genes coloured by \code{gene_group}.
#'           \item Arrow direction encoded by \code{Pdirection}.
#'           \item Optional gene labels inside arrows.
#'         }
#' @examples
#' \dontrun{
#' library(gggenes)
#' p <- gc_plot(
#'   data        = GC_meta,
#'   color_theme = c("bch" = "#3BAA51", "crt" = "#6495ED", "hypothetical ORF" = "grey")
#' )
#' print(p)
#' }
gc_plot = function(data = GC_meta,
                   color_theme = c('#3BAA51',
                                   '#6495ED',
                                   '#DD2421',
                                   '#EF9320',
                                   '#F8EB00',
                                   '#FF0683',
                                   '#956548',
                                   'grey')
){
  # data will be plotted in this function based on the geom_gene_arrow().
  # The numbers of Colors should be same with the gene_group.
  # Colors can be arbitrarily defined, e.g., color_theme_test= c('#EF9320','#EF5920', 'grey')
  color_theme = color_theme[1:length(unique(data$gene_group))]
  p=data %>%
    # {.[.$Pgenome %in% unique(.$Pgenome)[1:200], ]} %>%
    ggplot2::ggplot(.,
                    ggplot2::aes(xmin = Pstart, xmax = Pend,
                                 y = Pgenome, fill = gene_group,
                                 forward = Pdirection, label = gene_label)) +
    gggenes::geom_gene_arrow() +
    gggenes::geom_gene_label(fontface = 'bold') +
    ggplot2::scale_fill_manual(values=color_theme) +
    ggplot2::facet_wrap(~ Pgenome, scales = "free_y", ncol = 1) +
    gggenes::theme_genes() +
    ggplot2::theme(legend.position="bottom",
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank()
    )
  return(p)
}
