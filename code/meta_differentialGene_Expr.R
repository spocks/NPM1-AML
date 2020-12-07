library(survcomp)
library(ggplot2)

UHN = readRDS("../data/differentialGene_UHN.rds")
TCGA = readRDS("../data/differentialGene_TCGA.rds")
KI = readRDS("../data/differentialGene_KI.rds")
BeatAML = readRDS("../data/differentialGene_BeatAML.rds")
Leucegene = readRDS("../data/differentialGene_Leucegene.rds")

mde <- data.frame(
  GENENAME = UHN$GENENAME,
  estimate = NA,
  se = NA,
  p = NA,
  stringsAsFactors = F
)
rownames(mde) = mde$GENENAME

for (g in rownames(mde))
{
  getValue <-
    function(r, c) {
      c(UHN[r, c], TCGA[r, c], KI[r, c], BeatAML[r, c], Leucegene[r, c])
    }
  
  x <- getValue(g, "log2FoldChange")
  x.se <- getValue(g, "lfcSE")
  p <- getValue(g, "pvalue")
  v <- survcomp::combine.est(x, x.se, hetero = TRUE, na.rm = TRUE)
  pc <- survcomp::combine.test(p = p)
  mde[g, c("estimate", "se")] <- v[c("estimate", "se")]
  mde[g, "p"] <- pc
}

mde$fdr <- p.adjust(mde$p, "fdr")
write.csv(mde, file = "../results/meta_differentialGene.csv")

df = read.csv("../results/meta_differentialGene.csv", row.names = 1, 
              stringsAsFactors = F)

df$threshold = (abs(df$estimate) > 1 & df$fdr < 0.01)
col.txt = sprintf("FDR<%0.2f & fold change>%1.1f", 0.01, 1)
df$type = ifelse(df$threshold == T, col.txt, "Other")
df$type = factor(as.character(df$type), levels = c(col.txt, "Other"))
colx = c("#9E9D24FF", "gray")
names(colx) = c(col.txt, "Other")
g = ggplot(data = df, aes(
  x = estimate,
  y = -log10(fdr),
  colour = type
))
g = g + geom_point(alpha = 0.4, size = 1.75) + xlab("log2 fold change") + ylab("-log10 (FDR)")
g = g + scale_colour_manual(values = colx)
rng = max(abs(range(df$estimate)))
g = g + xlim(-rng, rng)
g = g + scale_y_continuous(trans = "log1p")
g = g + theme_light()
g = g + theme(
  legend.title = element_blank(),
  legend.background = element_rect(linetype = "solid", colour =
                                     "black")
)

pdf("../results/Supplementary_Figure-2_Volcano_plot.pdf",
    width = 8.27,
    height = 8.27)
print(g + theme(legend.position = "bottom", legend.box = "horizontal"))
dev.off()
