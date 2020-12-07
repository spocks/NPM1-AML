library(ggplot2)
library(ggpubr)

uhn.auc = readRDS("../data/UHN_drug_response.rds")

all.drugs = c("Sorafenib",
              "Ruxolitinib",
              "Sunitinib",
              "Quizartinib",
              "Imatinib",
              "Dasatinib")
colx = c("primitive" = '#d7191c', "committed" = '#2c7bb6')

pdf("../results/figure-6_UHN_drug_response.pdf",
    height = 6,
    width = 5)

for (drug in all.drugs)
{
  df = uhn.auc[uhn.auc$drug == drug,]
  pl = ggplot(df, aes(
    x = subtype,
    y = AUC,
    color = subtype,
    fill = subtype
  )) + geom_boxplot()
  pl = pl + scale_fill_manual(values = alpha(colx, 0.4))
  pl = pl + scale_color_manual(values = colx)
  pl = pl + stat_compare_means(comparisons = list(c("primitive", "committed")),
                               aes(label = "p={p.format}"),
                               method = "wilcox")
  pl = pl + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(drug)
  print(pl)
}

dev.off()


pdf(
  "../results/figure-6_BeatAML_drug_response.pdf",
  height = 6,
  width = 5
)

beatAML.auc = readRDS("../data/BeatAML_drug_response.rds")
for (drug in all.drugs)
{
  df = beatAML.auc[beatAML.auc$drug == drug,]
  pl = ggplot(df, aes(
    x = subtype,
    y = AUC,
    color = subtype,
    fill = subtype
  )) + geom_boxplot()
  pl = pl + scale_fill_manual(values = alpha(colx, 0.4))
  pl = pl + scale_color_manual(values = colx)
  pl = pl + stat_compare_means(comparisons = list(c("primitive", "committed")),
                               aes(label = "p={p.format}"),
                               method = "wilcox")
  pl = pl + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(drug)
  print(pl)
}

dev.off()
