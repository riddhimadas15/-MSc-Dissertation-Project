# Run once to install the R packages used in this project
install.packages(c("BiocManager", "dplyr", "readr", "readxl", "igraph", "ggraph", "ggplot2", "pheatmap", "Cairo"))
BiocManager::install(c("oligo", "limma", "multtest",
                       "hgu95av2.db", "pd.hg.u95av2",      # GSE10500 (HG-U95Av2)
                       "hgu133plus2.db", "pd.hgu133plus2")) # GSE93272 (HG-U133 Plus 2.0)
