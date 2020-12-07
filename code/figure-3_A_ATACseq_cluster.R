library(pheatmap)
path.result <- "../results/"
path.data <- "../data/ATACseq/"

BinaryMatrix <-
  readRDS(sprintf("%sCOREs/COREs_BinaryMatrix_OverlapBP0.rds", path.data))

rownames(BinaryMatrix) <-
  unlist(lapply(rownames(BinaryMatrix), function(X) {
    strsplit(X, "_")[[1]][2]
  }))
rownames(BinaryMatrix) <-
  paste("X", rownames(BinaryMatrix), sep = "")

JaccardMat <- c()
for (SamIter1 in 1:nrow(BinaryMatrix)) {
  JaccardVec <- c()
  for (SamIter2 in 1:nrow(BinaryMatrix)) {
    JaccardVec <-
      c(JaccardVec, length(which(BinaryMatrix[SamIter1,] == 1 &
                                   BinaryMatrix[SamIter2,] == 1)) / length(which(BinaryMatrix[SamIter1,] == 1 |
                                                                                   BinaryMatrix[SamIter2,] == 1)))
  }
  JaccardMat <- rbind(JaccardMat, JaccardVec)
}

JaccardMat_Z <- JaccardMat
JaccardMat_Z[which(JaccardMat_Z == 1)] <-
  max(JaccardMat_Z[which(JaccardMat_Z < 1)]) * 1.2
JaccardMat_Z <-
  (JaccardMat_Z - median(JaccardMat_Z)) / mad(JaccardMat_Z)

rownames(JaccardMat_Z) <- rownames(BinaryMatrix)
colnames(JaccardMat_Z) <- rownames(BinaryMatrix)

Categories <- readRDS(sprintf("%ssample_annotation.Rds", path.data))
Categories = Categories[colnames(JaccardMat_Z),]

CategoriesColor <- list(
  Subtype = c(Primtive = '#d7191c',
              Commited = '#2c7bb6'),
  'FLT3-ITD' = c(mutation = '#c51b7d',
                 wild = '#1a9641')
)

setwd(path.result)
pdf(
  paste(
    "figure-3_A_ATACseq_cluster.pdf",
    sep = "",
    collapse = ""
  ),
  onefile = F,
  width = 6,
  height = 5
)
par(mar = c(7, 3, 3, 3), mgp = c(4, 1, 0))
pheatmap::pheatmap(
  JaccardMat_Z,
  clustering_method = "ward.D2",
  annotation_col = Categories,
  annotation_colors = CategoriesColor,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  breaks = seq((min(JaccardMat_Z) - 0.01),
               (max(JaccardMat_Z) + 0.01), by = 0.01),
  color = colorRampPalette(c(
    rep("#009999", 6),
    rep("green", 1),
    rep("yellow", 3),
    rep("gold", 4),
    rep("orange", 2)
  ))(length(seq((min(JaccardMat_Z) - 0.01),
                (max(JaccardMat_Z) + 0.01), by = 0.01
  )) - 1)
)
dev.off()
