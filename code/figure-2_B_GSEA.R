library(piano)
library(snow)
library(BBmisc)
library(ggplot2)

GSAsummaryDF <- function(nes, ncpus = 1)
{
  nperm <- nes$nPerm
  gseares <- piano::GSAsummaryTable(nes)
  if ("p (dist.dir.up)" %in% colnames(gseares) &
      "p (dist.dir.dn)" %in% colnames(gseares))
  {
    gseares$p <- sapply(1:nrow(gseares), function(i)
    {
      v <- gseares[i , c("p (dist.dir.up)", "p (dist.dir.dn)")]
      min(v[!is.na(v)])
    })
    
    gseares[!is.na(gseares$p) &
              gseares$p == 0, "p"] <- 1 / (nperm + ncpus - (nperm %% ncpus) + 1)
    gseares$fdr <- p.adjust(gseares$p, "fdr")
    
    gseares <-
      BBmisc::sortByCol(gseares, c("p", "Stat (dist.dir)"), asc = c(T, F))
    rownames(gseares) <- NULL
    
    gseares <-
      gseares[, c("Name",
                  "Genes (tot)",
                  "Stat (dist.dir)",
                  "Genes (up)",
                  "Genes (down)",
                  "p",
                  "fdr")]
    
  }
  return(gseares)
}
run_GSEA <- function(gst,
                     gmt.file = NULL,
                     return.piano = FALSE,
                     ...)
{
  library(snow)
  gsa.opt <- list(...)
  
  if (is.null(gmt.file))
  {
    if ("gsa" %in% names(gsa.opt))
    {
      gsa <- gsa.opt[["gsa"]]
      gsa.opt[["gsa"]] <- NULL
    } else
    {
      stop("specify gmt.file or gsa")
    }
  } else
  {
    if ("gsa" %in% names(gsa.opt))
    {
      warning("gsa will be ignored as gmt.file is specified")
    }
    gsc <- piano::loadGSC(file = gmt.file)
  }
  
  if ("gsSizeLim" %in% names(gsa.opt))
  {
    gscSize <- sapply(gsc$gsc, length)
    gscSize <- gscSize[gscSize >= min(gsa.opt[["gsSizeLim"]])]
    gscSize <- gscSize[gscSize <= max(gsa.opt[["gsSizeLim"]])]
    totalGene <- unique(unlist(gsc$gsc[names(gscSize)]))
  } else
  {
    totalGene <- unique(unlist(gsc$gsc))
  }
  
  commanGene <- intersect(totalGene, names(gst))
  cat(
    sprintf(
      "total gene is gmt %d\ntotal gene in input data %d\ncomman gene %d\n",
      length(totalGene),
      length(gst),
      length(commanGene)
    )
  )
  
  nes <- piano::runGSA(geneLevelStats = gst, gsc = gsc, ...)
  
  if ("ncpus" %in% names(gsa.opt)) {
    ncpus <- gsa.opt[["ncpus"]]
  } else{
    ncpus <- 1
  }
  gseares <- GSAsummaryDF(nes, ncpus = ncpus)
  
  if (return.piano == TRUE)
  {
    return(list(data = gseares, piano = nes))
  }
  return(gseares)
}


mde = read.csv("../results/meta_differentialGene.csv", row.names = 1, 
               stringsAsFactors = F)
gene.st <- mde$estimate
names(gene.st) <- rownames(mde)

gmt.file = "../data/c2.cp.reactome.v6.0.symbols.gmt"
gsea.res = run_GSEA(
  gst = gene.st,
  gmt.file = gmt.file,
  geneSetStat = "gsea",
  nPerm = 1000,
  gseaParam = 1,
  ncpus = 1
)

write.csv(gsea.res, file = "../results/gene_set_enrichment_analysis.csv")


pdf("../results/figure-2_B_pathway.pdf",
    width = 7,
    height = 8)
df = gsea.res
df$Stat = df$`Stat (dist.dir)`
fdrc = 0.15
th = 0.3
df$threshold = (abs(df$`Stat (dist.dir)`) > 0.3 & df$fdr < fdrc)
col.txt = sprintf("FDR<%0.2f & Stat>%1.1f", fdrc, th)
df$type = ifelse(df$threshold == T, col.txt, "Other")
df$type = factor(as.character(df$type), levels = c(col.txt, "Other"))
colx = c("#e41a1c", "gray")
names(colx) = c(col.txt, "Other")
g = ggplot(data = df, aes(
  x = Stat,
  y = -log10(fdr),
  colour = type
))
g = g + geom_point(alpha = 0.4, size = 1.75) + xlab("enrichment score") + ylab("-log10 (FDR)")
g = g + scale_colour_manual(values = colx)
rng = max(abs(range(df$Stat)))
g = g + xlim(-rng, rng)
g = g + scale_y_continuous(trans = "log1p")
g = g + theme_light()
g = g + theme(
  legend.title = element_blank(),
  legend.background = element_rect(linetype = "solid", colour =
                                     "black")
)

print(g + theme(legend.position = "bottom", legend.box = "horizontal"))

dev.off()