#' @importFrom dplyr %>%
NULL
#' @importFrom utils packageVersion data
#' @export blastp_df
#' @export seq_data
#' @export photosynthesis_gene_list
#' @export PGC_group
NULL  # 仅用于生成NAMESPACE的占位符
#' @title Gene-Cluster Discovery, Annotation and Visualization
#'
#' @description
#' `gclink()` performs **all steps** from raw BLASTp output to a publication-ready
#' gene-cluster plot **in one call**. It dynamically adapts its workflow:
#' * If `in_seq_data` is provided: Extracts coordinates/sequences, merges data, and generates plots.
#' * If `in_seq_data` is `NULL`: Skips sequence-dependent steps, returning only the cluster table.
#' * If `in_GC_group` is `NULL`: Omits functional-group merging and plotting.
#' * Supports custom gene lists (e.g., photosynthesis, viral genes) via `in_gene_list` for universal
#'   cluster detection in bacteria, archaea, or phages.
#'
#' Gene clusters are identified via a density threshold (`AllGeneNum` and `MinConSeq`) due to
#' incomplete gene annotation coverage.
#'
#' @param in_blastp_df    A `data.frame` of Diamond BLASTp output with columns:
#'   \itemize{
#'     \item `qaccver`: Genome + contig name (separated by `"---"`), e.g.,
#'       `"Kuafubacteriaceae--GCA_016703535.1---JADJBV010000001.1_150"`:
#'       \itemize{
#'         \item Genome: `"Kuafubacteriaceae--GCA_016703535.1"` (`"--"` separator),
#'         \item Contig: `"JADJBV010000001.1"`,
#'         \item ORF: `"JADJBV010000001.1_150"` (`"_"` separator),
#'         \item Position: `"150"`.
#'       }
#'     \item `saccver`: Gene name + metadata (separated by `"_"`), e.g.,
#'       `"bchC_Methyloversatilis_sp_RAC08_BSY238_2447_METR"`:
#'       \itemize{
#'         \item Gene: `"bchC"`,
#'         \item Metadata: `"Methyloversatilis_sp_RAC08_BSY238_2447_METR"`.
#'       }
#'   }
#'   EggNOG results are supported by renaming annotation columns (e.g., `"GOs"`) to `saccver`.
#'   Default: `blastp_df`.
#'
#' @param in_seq_data     `NULL` or a `data.frame` with:
#'   \describe{
#'     \item{`SeqName`}{ORF identifier (Prodigal format: `>ORF_id # start # end # strand # ...`).}
#'     \item{`Sequence`}{ORF sequence.}
#'   }
#'   Example:
#'   `"Kuafubacteriaceae--GCA_016703535.1---JADJBV010000001.1_1 # 74 # 1018 # 1 # ..."`
#'   Can be imported from **Prodigal** FASTA using:
#'   ```r
#'   seq_data <- Biostrings::readBStringSet("Prodigal.fasta") %>%
#'     data.frame(Sequence = .) %>%
#'     tibble::rownames_to_column("SeqName")
#'   ```
#'   `NULL` skips coordinate extraction and plotting. Default: `seq_data`.
#'
#' @param in_gene_list    Character vector of reference genes for cluster detection.
#'                        Default: `photosynthesis_gene_list`.
#' @param in_GC_group     `NULL` or a `data.frame` mapping genes to functional groups
#'                        (columns: `gene`, `gene_group`). `NULL` skips plotting.
#'                        Default: `PGC_group`.
#' @param AllGeneNum      Integer; max ORFs per cluster. Default: `50`.
#' @param MinConSeq       Integer; min consecutive reference genes per cluster. Default: `25`.

#' @param down_IQR        Numeric; lower-bound scale factor for IQR length
#'                        filtering (see `length_filter`). Default: `10`.
#' @param up_IQR          Numeric; upper-bound scale factor for IQR length
#'                        filtering (see `length_filter`). Default: `10`.
#' @param orf_before_first Integer; number of hypothetical ORFs to insert
#'                        **before the first gene** of each detected cluster.
#'                        Useful when the upstream annotation is incomplete or
#'                        low-confidence; insertion is bounded by the first ORF
#'                        of the contig. Default: `0`.
#' @param orf_after_last  Integer; number of hypothetical ORFs to append
#'                        **after the last gene** of each detected cluster.
#'                        Helpful when the downstream annotation is incomplete
#'                        or low-confidence; insertion is bounded by the last
#'                        ORF of the contig. Default: `0`.
#' @param orf_range       Character. ORF inclusion mode:
#'                        \itemize{
#'                          \item `"All"`: Include all ORFs + annotations (default)
#'                          \item `"OnlyAnnotated"`: Keep only annotated ORFs
#'                          \item `"IgnoreAnnotated"`: Include all ORFs without annotation merging
#'                        }
#' @param levels_gene_group Character vector; factor levels for gene-group
#'                        legends (must include "hypothetical ORF" in case
#'                        some genes remain unclassified).
#'                        Ignored if `in_GC_group` is `NULL`.
#'                        Default: `c('bch','puh','puf','crt','acsF','assembly','regulator','hypothetical ORF')`.

#' @param color_theme     Character vector; colours for `gc_plot`, matched to
#'                        the length and order of `levels_gene_group`.
#'                        Ignored if plotting is skipped.
#'                        Default: `c('#3BAA51','#6495ED','#DD2421','#EF9320','#F8EB00','#FF0683','#956548','grey')`.
#' @param genome_subset   Character vector or `NULL`; genomes to retain.  If `NULL`, all genomes
#'                        are retained. Default: `NULL`.
#' @return A named list with:
#'   \describe{
#'     \item{`GC_meta`}{Annotated cluster table (`data.frame`).}
#'     \item{`GC_seq`}{FASTA sequences (if `in_seq_data` provided).}
#'     \item{`GC_plot`}{ggplot object (if `in_seq_data` and `in_GC_group` provided).}
#'   }
#' @export
#' @examples
#'    # Full pipeline (Find Cluster + Extract FASTA + Plot Cluster)
#'    data(blastp_df)
#'    data(seq_data)
#'    data(photosynthesis_gene_list)
#'    data(PGC_group)
#'    gc_list <- gclink(in_blastp_df = blastp_df,
#'                      in_seq_data = seq_data,
#'                      in_gene_list = photosynthesis_gene_list,
#'                      in_GC_group  = PGC_group,
#'                      AllGeneNum = 50,
#'                      MinConSeq  = 25,
#'                      down_IQR   = 10,
#'                      up_IQR     = 10,
#'                      orf_before_first = 0,
#'                      orf_after_last = 0,
#'                      levels_gene_group = c('bch','puh','puf','crt','acsF','assembly','regulator',
#'                                            'hypothetical ORF'),
#'                      color_theme = c('#3BAA51','#6495ED','#DD2421','#EF9320','#F8EB00',
#'                                      '#FF0683','#956548','grey'),
#'                      genome_subset = NULL)
#'    gc_meta = gc_list[["GC_meta"]]
#'    gc_seq = gc_list[["GC_seq"]]
#'    gc_plot = gc_list[["GC_plot"]]
#'    head(gc_meta)
#'    head(gc_seq)
#'    plot(gc_plot)
gclink <- function(in_blastp_df = blastp_df,
                   in_seq_data = seq_data,
                   in_gene_list = photosynthesis_gene_list,
                   in_GC_group  = PGC_group,
                   AllGeneNum = 50,
                   MinConSeq = 25,
                   down_IQR = 10,
                   up_IQR = 10,
                   orf_before_first = 0,
                   orf_after_last = 0,
                   orf_range = "All",
                   levels_gene_group = c('bch','puh','puf','crt','acsF','assembly','regulator',
                                         'hypothetical ORF'),
                   color_theme = c('#3BAA51','#6495ED','#DD2421','#EF9320','#F8EB00',
                                   '#FF0683','#956548','grey'),
                   genome_subset = NULL
){
  # 1 filter based on gene length from blastp result
  bin_genes <- in_blastp_df %>%
    orf_extract() %>%
    length_filter(down_IQR = down_IQR, up_IQR = up_IQR)

  # 2 find gene cluster
  sbgc <- gc_cal(Data = bin_genes,
                 in_gene_list = in_gene_list,
                 AllGeneNum = AllGeneNum,
                 MinConSeq  = MinConSeq)

  # 3 supplement ORF
  sbgc_add <- gc_add(Data = sbgc, Annotation = bin_genes, orf_before_first = orf_before_first, orf_after_last = orf_after_last, orf_range = orf_range)

  if (!is.null(in_seq_data)) {
    # 4 get ORF loci from CDS FASTA
    GC_location <- orf_locate(in_seq_data)

    # 5 merge orf_position and annotion
    GC_gene = merge(sbgc_add, GC_location, by='qaccver', all.x = TRUE) %>%
      {.[order(.$gene_cluster, .$orf_position), ]}
    GC_gene$gene[is.na(GC_gene$gene)] <- 'unrelated'

    # 6 scale FASTA
    GC_seq <- GC_gene %>%
      dplyr::mutate(gene_name_seq = paste0(">", gene, "_", gene_cluster,
                                           "\n", Sequence)) %>%
      {.[!(is.na(.$Sequence)),]} %>%
      {.$gene_name_seq}

    if (!is.null(in_GC_group)) {
      # 7 merge GC_group
      GC_meta  <- merge(GC_gene, in_GC_group, by = 'gene', all.x = TRUE)

      # 8 filter based on genome_subset
      if (!is.null(genome_subset)) {
        GC_meta <- base::subset(GC_meta, genome %in% genome_subset)
      }

      # 9 scale & plot
      GC_meta <- gc_scale(GC_meta, levels_gene_group = levels_gene_group)
      GC_plot  <- gc_plot(GC_meta, color_theme = color_theme)

      # 10 reture list
      list(
        GC_meta = GC_meta,
        GC_seq  = GC_seq,
        GC_plot = GC_plot
      )
    } else {
      list(
        GC_meta = GC_gene,
        GC_seq  = GC_seq,
        GC_plot = NULL
      )
    }
  } else {
    list(
      GC_meta = sbgc_add,
      GC_seq  = NULL,
      GC_plot = NULL
    )
  }
}
