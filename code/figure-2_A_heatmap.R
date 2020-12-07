library(BBmisc)
library(Biobase)
library(ComplexHeatmap)
library(circlize)
source("important_fun.R")
getData <- function(expressionFile, diffExpFile)
{
  dffExp = readRDS(diffExpFile)
  dffExp = sortByCol(dffExp, c("log2FoldChange"))
  topGenes = c(rownames(dffExp)[1:150], rev(rownames(dffExp))[1:150])
  df = readRDS(expressionFile)
  df = df[topGenes,]
  return(df)
}

plotHeatmap <- function(dt, geneList, title = "")
{
  exprs(dt) = t(scale(t(exprs(dt)))[, ])
  colSc = c(
    "#440154FF",
    "#482576FF",
    "#414487FF",
    "#35608DFF",
    "#2A788EFF",
    "#21908CFF",
    "#22A884FF",
    "#43BF71FF",
    "#7AD151FF",
    "#BBDF27FF",
    "#FDE725FF"
  )
  colPN <- c("positive" = "#d01c8b", "negative" = "#CCCCCC")
  colClust <- c("primitive" = "#e41a1c", "committed" = "#377eb8")
  colList <- list("subtype" = colClust,
                  "FLT3-ITD" = colPN,
                  "FLT3-TKD" = colPN)
  legTextPara <-
    list(
      title_gp = gpar(fontsize = 10, fontface = "plain"),
      labels_gp = gpar(fontsize = 9)
    )
  
  topAnn <- HeatmapAnnotation(
    df = pData(dt)[, names(colList)],
    show_annotation_name = T,
    col = colList,
    gp = gpar(col = "white", lwd = 0.01),
    annotation_name_gp = gpar(fontsize = 9, fontface = "plain"),
    annotation_legend_param = legTextPara
  )
  
  rng = c(-1.27,  1.27)
  colPal = colorRamp2(seq(rng[1], rng[2], length.out = length(colSc)), colSc)
  
  hlpara <- list(
    title = "Expression",
    color_bar = "continuous",
    at = range(attr(colPal, "breaks")),
    labels = c("low", "high")
  )
  
  hlpara = c(hlpara, legTextPara)
  
  ht <-
    Heatmap(
      exprs(dt),
      column_title = title,
      column_title_gp = gpar(fontsize = 11),
      col = colPal,
      cluster_rows = FALSE,
      cluster_columns = F,
      show_column_names = F,
      show_row_names = F,
      top_annotation = topAnn,
      heatmap_legend_param = hlpara,
      rect_gp = gpar(col = "white", lwd = 0.05)
    )
  
  cg.df = data.frame(present = ifelse(geneList[featureNames(dt)] == 5, 
                                      "TRUE", "FALSE"))
  rowAnn <- HeatmapAnnotation(
    df = cg.df,
    col = list(present = c("FALSE" = "white", "TRUE" = "#737373")),
    which = "row",
    width = unit(0.25, "cm"),
    show_legend = F,
    show_annotation_name = F
  )
  ht <- rowAnn + ht
  return(ht)
}


UHN = getData(expressionFile = "../data/UHN_geneExpression.rds",
              diffExpFile = "../data/differentialGene_UHN.rds")

TCGA = getData(expressionFile = "../data/TCGA_geneExpression.rds",
               diffExpFile = "../data/differentialGene_TCGA.rds")

KI = getData(expressionFile = "../data/KI_geneExpression.rds",
             diffExpFile = "../data/differentialGene_KI.rds")

BeatAML = getData(expressionFile = "../data/BeatAML_geneExpression.rds",
                  diffExpFile = "../data/differentialGene_BeatAML.rds")

Leucegene = getData(expressionFile = "../data/Leucegene_geneExpression.rds",
                    diffExpFile = "../data/differentialGene_Leucegene.rds")

geneList <-
  table(unlist(lapply(list(UHN, TCGA, KI, BeatAML, Leucegene),
                      function(x)
                        rownames(x))))


png("../results/figure-2_A_heatmap_UHN.png", width = 5, height = 6, units = "in", res = 300)
plotHeatmap(UHN, geneList, title = "UHN")
dev.off()

png("../results/figure-2_A_heatmap_TCGA.png", width = 5, height = 6, units = "in", res = 300)
plotHeatmap(TCGA, geneList, title = "TCGA")
dev.off()

png("../results/figure-2_A_heatmap_KI.png", width = 5, height = 6, units = "in", res = 300)
plotHeatmap(KI, geneList, title = "KI")
dev.off()

png("../results/figure-2_A_heatmap_BeatAML.png", width = 5, height = 6, units = "in", res = 300)
plotHeatmap(BeatAML, geneList, title = "BeatAML")
dev.off()

png("../results/figure-2_A_heatmap_Leucegene.png", width = 5, height = 6, units = "in", res = 300)
plotHeatmap(Leucegene, geneList, title = "Leucegene")
dev.off()
