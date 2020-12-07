library(Biobase)
library(ggplot2)

source("important_fun.R")
ccRepet = 10

outpath = "../results/uhn_cluster"
dir.create(outpath)
uhn = readRDS("../data/UHN_geneExpression.rds")
UHN_silhouette = heirarchical_consensClust(
  t(exprs(uhn))[, 1:5000],
  outFileName = sprintf("%s/uhn_cluster", outpath),
  distance = "spearman",
  pItem = 0.75,
  pFeature = 0.75,
  reps = ccRepet,
  maxK = 8,
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2"
)


outpath = "../results/TCGA_cluster"
dir.create(outpath)
tcga = readRDS("../data/TCGA_geneExpression.rds")
TCGA_silhouette = heirarchical_consensClust(
  t(exprs(tcga))[, 1:5000],
  outFileName = sprintf("%s/TCGA_cluster", outpath),
  distance = "spearman",
  pItem = 0.75,
  pFeature = 0.75,
  reps = ccRepet,
  maxK = 8,
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2"
)

outpath = "../results/BeatAML_cluster"
dir.create(outpath)
beatAML = readRDS("../data/BeatAML_geneExpression.rds")
BeatAML_silhouette = heirarchical_consensClust(
  t(exprs(beatAML))[, 1:5000],
  outFileName = sprintf("%s/BeatAML_cluster", outpath),
  distance = "spearman",
  pItem = 0.75,
  pFeature = 0.75,
  reps = ccRepet,
  maxK = 8,
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2"
)


outpath = "../results/KI_cluster"
dir.create(outpath)
KI = readRDS("../data/KI_geneExpression.rds")
KI_silhouette = heirarchical_consensClust(
  t(exprs(KI))[, 1:5000],
  outFileName = sprintf("%s/KI_cluster", outpath),
  distance = "spearman",
  pItem = 0.75,
  pFeature = 0.75,
  reps = ccRepet,
  maxK = 8,
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2"
)

outpath = "../results/Leucegene_cluster"
dir.create(outpath)
Leucegene = readRDS("../data/Leucegene_geneExpression.rds")
Leucegene_silhouette = heirarchical_consensClust(
  t(exprs(Leucegene))[, 1:5000],
  outFileName = sprintf("%s/Leucegene_cluster", outpath),
  distance = "spearman",
  pItem = 0.75,
  pFeature = 0.75,
  reps = ccRepet,
  maxK = 8,
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2"
)


##---------- plot silhouette for clusters --------------------------------------
df = data.frame(
  dataset = c(
    rep("UHN", 7),
    rep("TCGA", 7),
    rep("BeatAML", 7),
    rep("KI", 7),
    rep("Leucegene", 7)
  ),
  clusters = c(
    names(UHN_silhouette),
    names(TCGA_silhouette),
    names(BeatAML_silhouette),
    names(KI_silhouette),
    names(Leucegene_silhouette)
  ),
  silhouette = c(
    UHN_silhouette,
    TCGA_silhouette,
    BeatAML_silhouette,
    KI_silhouette,
    Leucegene_silhouette
  ),
  stringsAsFactors = F
)

df$data <-
  factor(df$data, level = c("UHN", "TCGA", "BeatAML",  "KI",  "Leucegene"))
df$clusters <- as.numeric(gsub("K", "", df$clusters))
plt <-
  ggplot(data = df, aes(x = clusters, y = silhouette, colour = dataset))
plt <- plt + geom_point(shape = 21, fill = "white", size = 3)
plt <- plt + geom_line(size = 1)
plt <- plt + theme_bw()
plt <- plt + scale_x_continuous(breaks = unique(df$clusters),
                                labels = paste0("K=", unique(df$clusters)))
plt <- plt + scale_y_continuous(name = "average silhouette width")


pdf("../results/Supplementary_Figure-1_subtype.pdf", 
    width = 8,
    height = 5)
print(plt)
dev.off()
