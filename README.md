```markdown
# gclink: Gene-Cluster Discovery, Annotation and Visualization

## Overview

`gclink` performs end-to-end analysis of gene clusters (e.g., photosynthesis,<br>
carbon/nitrogen/sulfur cycling, carotenoid, antibiotic, or viral genes)<br>
from (meta)genomes or (meta)transcriptomes. It provides:

- Parsing of Diamond/BLASTp-style results (including EggNOG output)
- Contiguous cluster detection
- Publication-ready visualization

## Installation

```r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("LiuyangLee/gclink")
```

## Quick Start

```r
library(gclink)

# Full pipeline example
data(blastp_df)
data(seq_data)
data(photosynthesis_gene_list)
data(PGC_group)

results <- gclink(
  in_blastp_df = blastp_df,
  in_seq_data = seq_data,
  in_gene_list = photosynthesis_gene_list,
  in_GC_group = PGC_group
)

# Access results
head(results$GC_meta)  # Cluster metadata
head(results$GC_seq)   # FASTA sequences
plot(results$GC_plot)  # Visualization
```

## Key Features

### Adaptive Workflow
- Works with or without CDS FASTA input
- Skips plotting when functional grouping is absent
- Supports custom gene lists for universal cluster detection

### Cluster Detection
- Density-based identification via `AllGeneNum` and `MinConSeq` parameters
- Handles incomplete gene annotation coverage
- Optional insertion of hypothetical ORFs at cluster boundaries

### Visualization
- Publication-ready arrow plots with customizable based on `gggenes`:
  - Color themes
  - Functional group levels
  - Genome subsets

## Documentation

Full function reference:
```r
?gclink::gclink
```

## Citation

If you use `gclink` in your research, please cite:

> Li, L., Huang, D., Hu, Y., Rudling, N. M., Canniffe, D. P., Wang, F., & Wang, Y.
> "Globally distributed Myxococcota with photosynthesis gene clusters illuminate the origin and evolution of a potentially chimeric lifestyle."
> *Nature Communications* (2023), 14, 6450.
> https://doi.org/10.1038/s41467-023-42193-7

## Dependencies

- R (≥ 3.5)
- dplyr (≥ 1.1.4)
- gggenes (≥ 0.5.1)
- ggplot2 (≥ 3.5.2)

## License

GPL-3 © [Liuyang Li](https://orcid.org/0000-0001-6004-9437)

## Contact

- Maintainer: Liuyang Li <cyanobacteria@yeah.net>
- Bug reports: https://github.com/LiuyangLee/gclink/issues
```
