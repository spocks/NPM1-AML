library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

GenomicDist <- function(gr1) {
  if (require(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
    aCR <- assignChromosomeRegion(
      gr1,
      nucleotideLevel = FALSE,
      precedence = c(
        "Promoters",
        "immediateDownstream",
        "fiveUTRs",
        "threeUTRs",
        "Exons",
        "Introns"
      ),
      # "Enhancer.Silencer"),
      TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
      proximal.promoter.cutoff = 1000,
      immediate.downstream.cutoff = 1000
    )#TxDb.Hsapiens.UCSC.hg19.knownGene)
    return(aCR$percentage)
  }
}

######################################
path.result <- "../results/"
path.data <- "../data/ATACseq/"
ClustDir <- "../data/ATACseq/COREs/"

ClusFiles <-
  list.files(path = ClustDir,
             pattern = "_CREAM.bed",
             full.names = T)

SamNames <-
  unlist(lapply(ClusFiles, function(X) {
    strsplit(X, "_")[[1]][2]
  }))

StemSams = c("X110661", "X120170", "X120170", "X120541", "X120899", "X130023", 
             "X140659", "X8107", "X90376", "X90669")

GenomicDistMat <- c()
for (ClustFile in ClusFiles) {
  ClustTable <- read.table(ClustFile)
  bed <- as.data.frame(ClustTable[, 1:3])
  colnames(bed) <- c("seqnames", "start", "end")
  gr1 <- toGRanges(bed, format = "BED", header = FALSE)
  GenomicDistMat <- cbind(GenomicDistMat, GenomicDist(gr1))
}

colnames(GenomicDistMat) <- paste("X", SamNames, sep = "")

GenomicDistMat <-
  cbind(GenomicDistMat[, which(colnames(GenomicDistMat) %in% StemSams)],
        GenomicDistMat[, -which(colnames(GenomicDistMat) %in% StemSams)])

PvalVec <- c()
for (DistType in 1:nrow(GenomicDistMat)) {
  PvalVec <-
    c(PvalVec,
      wilcox.test(GenomicDistMat[DistType, which(colnames(GenomicDistMat) %in% StemSams)],
                  GenomicDistMat[DistType, -which(colnames(GenomicDistMat) %in% StemSams)])$p.value)
}

names(PvalVec) <- rownames(GenomicDistMat)
FDRVec <- p.adjust(PvalVec)
names(FDRVec) <- rownames(GenomicDistMat)
print(FDRVec)
##################################
##################################
setwd(path.result)
pdf(
  paste(
    "figure-3B_COREs_Promoters.pdf",
    sep = "",
    collapse = ""
  ),
  onefile = F,
  width = 5,
  height = 8
)
par(mar = c(6, 7, 3, 3), mgp = c(3.5, 1, 0))
barplot(
  sort(GenomicDistMat["Promoters", ] - 60, decreasing = F),
  col =
    c(rep('#d7191c', length(
      which(colnames(GenomicDistMat) %in% StemSams)
    )),
    rep('#2c7bb6', (
      ncol(GenomicDistMat) - length(which(colnames(GenomicDistMat) %in% StemSams))
    )))[sort(GenomicDistMat["Promoters", ] - 60,
             decreasing = F,
             index.return = T)[[2]]],
  xlim = c(0, 30),
  xaxt = "n",
  xlab = "Percentage of COREs\n at gene promoters",
  las = 2,
  cex.names = 1.2,
  cex.lab = 1.5,
  horiz = T
)
axis(
  side = 1,
  at = c(0, 15, 30),
  labels = 60 + c(0, 15, 30),
  cex.axis = 1.5
)
dev.off()


pdf(
  paste(
    "figure-3B_COREs_Intergenic_Region.pdf",
    sep = "",
    collapse = ""
  ),
  onefile = F,
  width = 5,
  height = 8
)
par(mar = c(6, 7, 3, 3), mgp = c(3.5, 1, 0))
barplot(
  sort(GenomicDistMat["Introns", ], decreasing = F),
  col =
    c(rep('#d7191c', length(
      which(colnames(GenomicDistMat) %in% StemSams)
    )),
    rep('#2c7bb6', (
      ncol(GenomicDistMat) - length(which(colnames(GenomicDistMat) %in% StemSams))
    )))[sort(GenomicDistMat["Introns", ],
             decreasing = F,
             index.return = T)[[2]]],
  xlim = c(0, 15),
  xaxt = "n",
  xlab = "Percentage of COREs\n in intergenic regions",
  las = 2,
  cex.names = 1.2,
  cex.lab = 1.5,
  horiz = T
)
axis(
  side = 1,
  at = c(0, 7.5, 15),
  labels = c(0, 7.5, 15),
  cex.axis = 1.5
)
dev.off()


pdf(
  paste(
    "Supple_Fig-3_COREs_GenomicDist_Intergenic.Region.pdf",
    sep = "",
    collapse = ""
  ),
  onefile = F,
  width = 5,
  height = 8
)
par(mar = c(6, 7, 3, 3), mgp = c(3.5, 1, 0))
barplot(
  sort(GenomicDistMat["Intergenic.Region", ], decreasing = F),
  col =
    c(rep('#d7191c', length(
      which(colnames(GenomicDistMat) %in% StemSams)
    )),
    rep('#2c7bb6', (
      ncol(GenomicDistMat) - length(which(colnames(GenomicDistMat) %in% StemSams))
    )))[sort(GenomicDistMat["Intergenic.Region", ],
             decreasing = F,
             index.return = T)[[2]]],
  xlim = c(0, 18),
  xaxt = "n",
  xlab = "Percentage of COREs\n in intergenic regions",
  las = 2,
  cex.names = 1.2,
  cex.lab = 1.5,
  horiz = T
)
axis(
  side = 1,
  at = c(0, 9, 18),
  labels = c(0, 9, 18),
  cex.axis = 1.5
)
dev.off()


pdf(
  paste(
    "Supple_Fig-3_COREs_GenomicDist_immediateDownstream.pdf",
    sep = "",
    collapse = ""
  ),
  onefile = F,
  width = 5,
  height = 8
)
par(mar = c(6, 7, 3, 3), mgp = c(3.5, 1, 0))
barplot(
  sort(GenomicDistMat["immediateDownstream", ], decreasing = F),
  col =
    c(rep('#d7191c', length(
      which(colnames(GenomicDistMat) %in% StemSams)
    )),
    rep('#2c7bb6', (
      ncol(GenomicDistMat) - length(which(colnames(GenomicDistMat) %in% StemSams))
    )))[sort(GenomicDistMat["immediateDownstream", ],
             decreasing = F,
             index.return = T)[[2]]],
  xlim = c(0, 5),
  xaxt = "n",
  xlab = "Percentage of COREs\n in intergenic regions",
  las = 2,
  cex.names = 1.2,
  cex.lab = 1.5,
  horiz = T
)
axis(
  side = 1,
  at = c(0, 2.5, 5),
  labels = c(0, 2.5, 5),
  cex.axis = 1.5
)
dev.off()



pdf(
  paste(
    "Supple_Fig-3_COREs_GenomicDist_fiveUTRs.pdf",
    sep = "",
    collapse = ""
  ),
  onefile = F,
  width = 5,
  height = 8
)
par(mar = c(6, 7, 3, 3), mgp = c(3.5, 1, 0))
barplot(
  sort(GenomicDistMat["fiveUTRs", ], decreasing = F),
  col =
    c(rep('#d7191c', length(
      which(colnames(GenomicDistMat) %in% StemSams)
    )),
    rep('#2c7bb6', (
      ncol(GenomicDistMat) - length(which(colnames(GenomicDistMat) %in% StemSams))
    )))[sort(GenomicDistMat["fiveUTRs", ],
             decreasing = F,
             index.return = T)[[2]]],
  xlim = c(0, 3),
  xaxt = "n",
  xlab = "Percentage of COREs\n in intergenic regions",
  las = 2,
  cex.names = 1.2,
  cex.lab = 1.5,
  horiz = T
)
axis(
  side = 1,
  at = c(0, 1.5, 3),
  labels = c(0, 1.5, 3),
  cex.axis = 1.5
)
dev.off()

pdf(
  paste(
    "Supple_Fig-3_COREs_GenomicDist_threeUTRs.pdf",
    sep = "",
    collapse = ""
  ),
  onefile = F,
  width = 5,
  height = 8
)
par(mar = c(6, 7, 3, 3), mgp = c(3.5, 1, 0))
barplot(
  sort(GenomicDistMat["threeUTRs", ], decreasing = F),
  col =
    c(rep('#d7191c', length(
      which(colnames(GenomicDistMat) %in% StemSams)
    )),
    rep('#2c7bb6', (
      ncol(GenomicDistMat) - length(which(colnames(GenomicDistMat) %in% StemSams))
    )))[sort(GenomicDistMat["threeUTRs", ],
             decreasing = F,
             index.return = T)[[2]]],
  xlim = c(0, 1),
  xaxt = "n",
  xlab = "Percentage of COREs\n in intergenic regions",
  las = 2,
  cex.names = 1.2,
  cex.lab = 1.5,
  horiz = T
)
axis(
  side = 1,
  at = c(0, 0.5, 1),
  labels = c(0, 0.5, 1),
  cex.axis = 1.5
)
dev.off()

pdf(
  paste(
    "Supple_Fig-3_COREs_GenomicDist_Exons.pdf",
    sep = "",
    collapse = ""
  ),
  onefile = F,
  width = 5,
  height = 8
)
par(mar = c(6, 7, 3, 3), mgp = c(3.5, 1, 0))
barplot(
  sort(GenomicDistMat["Exons", ], decreasing = F),
  col =
    c(rep('#d7191c', length(
      which(colnames(GenomicDistMat) %in% StemSams)
    )),
    rep('#2c7bb6', (
      ncol(GenomicDistMat) - length(which(colnames(GenomicDistMat) %in% StemSams))
    )))[sort(GenomicDistMat["Exons", ],
             decreasing = F,
             index.return = T)[[2]]],
  xlim = c(0, 4),
  xaxt = "n",
  xlab = "Percentage of COREs\n in intergenic regions",
  las = 2,
  cex.names = 1.2,
  cex.lab = 1.5,
  horiz = T
)
axis(
  side = 1,
  at = c(0, 2, 4),
  labels = c(0, 2, 4),
  cex.axis = 1.5
)
dev.off()
