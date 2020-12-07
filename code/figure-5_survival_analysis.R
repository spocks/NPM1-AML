library(survminer)
library(survival)
library(forestmodel)

dfx = readRDS("../data/survival_data.rds")

colx = c("#e41a1c", "#377eb8")
names(colx) =  c("primitive", "committed")

pdf("../results/figure-5A_survival_plots.pdf",
    height = 5,
    width = 5)

##----- for complete dataset ------------------------------
fit = survfit(Surv(OS, OS.status) ~ subtype, data = dfx)
plt = ggsurvplot(
  fit,
  data = dfx,
  size = 1,
  palette = c("#e41a1c", "#377eb8"),
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  fontsize = 2.75,
  risk.table.col = "strata",
  legend.labs = c("primitive", "committed"),
  risk.table.height = 0.22,
  xlab = "Time (in months)",
  ggtheme = theme_classic(),
  legend.title = "subtype",
  legend = "right",
  title = ""
)

plt$plot <-
  plt$plot + scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
plt$plot <-
  plt$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
plt$table <- plt$table + theme(axis.title.y = element_blank())
print(plt)
dev.off()



##----- for individual dataset ------------------------------
pdf("../results/Supplementary_Figure-5_survival_plots.pdf",
    height = 5,
    width = 5)
for (dataset in c("UHN", "TCGA", "Leucegene", "KI", "BeatAML"))
{
  fit <-
    survfit(Surv(OS, OS.status) ~ subtype, data = dfx[dfx$dataset == dataset,])
  plt <-
    ggsurvplot(
      fit,
      data = dfx[dfx$dataset == dataset,],
      size = 1,
      palette = c("#e41a1c", "#377eb8"),
      conf.int = TRUE,
      pval = TRUE,
      risk.table = TRUE,
      fontsize = 2.75,
      risk.table.col = "strata",
      legend.labs = c("primitive", "committed"),
      risk.table.height = 0.22,
      xlab = "Time (in months)",
      ggtheme = theme_classic(),
      legend.title = "subtype",
      legend = "right",
      title = dataset
    )
  
  plt$plot <-
    plt$plot + scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
  plt$plot <-
    plt$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  plt$table <- plt$table + theme(axis.title.y = element_blank())
  print(plt)
}

dev.off()

##------- do Cox analysis ----
res.cox = coxph(
  Surv(OS, OS.status) ~ subtype + Sex + WBC + Age + Karyotype +
    Transplant + FLT3.ITD + FLT3.TKD + DNMT3A + KRAS + NRAS +
    subtype * FLT3.ITD + strata(dataset),
  data = dfx
)

pdf("../results/figure-5B_CoxPH.pdf",
    height = 7,
    width = 10)
print(forest_model(res.cox))
dev.off()
