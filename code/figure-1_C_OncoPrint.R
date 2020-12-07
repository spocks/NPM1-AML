library(ComplexHeatmap)


make_side_ann <- function(mat, mat.info)
{
  rtx <-
    data.frame(
      name = rownames(mat),
      primitive = NA,
      committed = NA,
      stringsAsFactors = F
    )
  rownames(rtx) <- rtx$name
  for (gx in rownames(rtx))
  {
    v <- table(mat[gx,], mat.info[, "subtype"])
    rtx[gx, "primitive"] <- v["mutation", "primitive"]
    rtx[gx, "committed"] <- v["mutation", "committed"]
  }
  rty <- as.matrix(rtx[, c("primitive", "committed")])
  
  clustCol <-
    c('#b2182b', '#2166ac')
  names(clustCol) <- c("primitive", "committed")
  rowbar <-
    ComplexHeatmap::anno_barplot(
      rty,
      which = "row",
      axis = TRUE,
      border = FALSE,
      axis_param = list(side = "top"),
      gp = gpar(fill = clustCol, col =
                  NA)
    )
  row_ha <- rowAnnotation(
    foo = rowbar,
    show_annotation_name = F,
    width = unit(1.5, "cm")
  )
  return(row_ha)
}

make_top_ann <- function(topBarmat, lp)
{
  clustCol <- c('#b2182b', '#2166ac')
  names(clustCol) <- c("primitive", "committed")
  
  dataCol <- c("#CCDD99", "#98BA58",  "#76a21e", "#446520")
  names(dataCol) <- c("TCGA", "KI", "BeatAML", "Leucegene")
  topAnn <-
    HeatmapAnnotation(
      Subtype = topBarmat$subtype,
      Data = topBarmat$data,
      col = list("Subtype" = clustCol,  "Data" =
                   dataCol),
      show_annotation_name = T,
      annotation_name_gp = gpar(fontsize = 9),
      annotation_legend_param = lp,
      height = unit(0.7, "cm")
    )
  return(topAnn)
}

alter_fun_Onco = list(
  background = function(x, y, w, h)
  {
    grid.rect(x, y, w * 0.70, h * 0.95, gp = gpar(fill = bcCol, col = lineCol))
  },
  mutation = function(x, y, w, h) {
    grid.rect(x, y, w * 0.70, h * 0.95, gp = gpar(fill = mutCol, col = lineCol))
  }
)


mutationData = readRDS("../data/mutation.rds")
mat = mutationData$mutation
mat.info = mutationData$meta

row_barplot <- make_side_ann(mat, mat.info)
lp <-
  list(title_gp = gpar(fontsize = 9),
       labels_gp = gpar(fontsize = 8))
topAnn <- make_top_ann(mat.info[, c("subtype", "data")], lp)

mutCol = "#854E4B"
names(mutCol) = c("mutation")
bcCol = "#D3D3D3"
lineCol = NA

ht <- oncoPrint(
  mat,
  get_type = function(x)
    strsplit(x, ";")[[1]],
  alter_fun = alter_fun_Onco,
  col = mutCol,
  row_order = NULL,
  column_order = colnames(mat),
  show_pct = T,
  show_row_names = T,
  right_annotation = NULL,
  top_annotation = topAnn,
  pct_gp = gpar(fontsize = 7, fontface = "plain"),
  row_names_gp = gpar(fontsize = 7, fontface = "plain"),
  heatmap_legend_param = append(list(title = "Alterations"), lp)
)


png(
  "../results/figure-1_C_Oncoprint.png",
  width = 10.27,
  height = 3.3,
  units = 'in',
  res = 400
)
draw(ht + row_barplot,
     merge_legends = TRUE,
     heatmap_legend_side = "right")
dev.off()
