##---- to setup the R environment

CRANpackages <-
  c(
    "BBmisc",
    "Hmisc",
    "circlize",
    "cluster",
    "forestmodel",
    "ggplot2",
    "ggpubr",
    "igraph",
    "pheatmap",
    "snow",
    "survival",
    "survminer",
    "BiocManager"
  )

cpac <- setdiff(CRANpackages, rownames(installed.packages()))
if (length(cpac) > 0)
{
  cat(
    sprintf(
      "These packages are not installed:\n%s\n\nInstalling the packages\n\n",
      paste0(cpac, collapse = ",\n")
    )
  )
  install.packages(cpac)
}


bioCpackages <- c(
  "Biobase",
  "ChIPpeakAnno",
  "ComplexHeatmap",
  "ConsensusClusterPlus",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "piano",
  "survcomp"
)

bpac <- setdiff(bioCpackages, rownames(installed.packages()))
if (length(bpac) > 0)
{
  cat(
    sprintf(
      "These packages are not installed:\n%s\n\nInstalling the packages\n\n",
      paste0(bpac, collapse = ",\n")
    )
  )
  BiocManager::install(bpac)
}
