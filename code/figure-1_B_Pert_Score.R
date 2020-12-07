library(ggplot2)
library(ggpubr)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

colClust = c("primitive"='#d7191c', "committed"='#2c7bb6')


pdf("../results/figure-1_B_Pert_score.pdf", height = 7, width = 5)

uhn = readRDS("../data/PERT_Score_UHN.rds")
plt.uhn = ggplot(uhn, aes(x=cluster, y=score, fill=cluster))
plt.uhn = plt.uhn + geom_violin(trim=F, color=NA, scale="width")
plt.uhn = plt.uhn + stat_summary(fun.data=data_summary, geom="pointrange", color="white")
plt.uhn = plt.uhn + stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x.npc='center', label.y.npc = "top",
                     comparisons = list(c("primitive", "committed")))
plt.uhn = plt.uhn + scale_fill_manual(values=colClust) + theme_classic()
print(plt.uhn+theme(legend.position="none")+xlab("subtype")+ ggtitle("UHN"))

TCGA = readRDS("../data/PERT_Score_TCGA.rds")
plot.tcga = ggplot(TCGA, aes(x=cluster, y=score, fill=cluster))
plot.tcga = plot.tcga + geom_violin(trim=F, color=NA, scale="width")
plot.tcga = plot.tcga + stat_summary(fun.data=data_summary, geom="pointrange", color="white")
plot.tcga = plot.tcga + scale_fill_manual(values=colClust) + theme_classic()
plot.tcga = plot.tcga + stat_compare_means(method = "wilcox.test", label = "p.signif",
                                       label.x.npc='center', label.y.npc = "top",
                                       comparisons = list(c("primitive", "committed")))
print(plot.tcga+theme(legend.position="none")+xlab("subtype")+ ggtitle("TCGA"))


KI = readRDS("../data/PERT_Score_KI.rds")
plot.ki = ggplot(KI, aes(x=cluster, y=score, fill=cluster))
plot.ki = plot.ki + geom_violin(trim=F, color=NA, scale="width")
plot.ki = plot.ki + stat_summary(fun.data=data_summary, geom="pointrange", color="white")
plot.ki = plot.ki + scale_fill_manual(values=colClust) + theme_classic()
plot.ki = plot.ki + stat_compare_means(method = "wilcox.test", label = "p.signif",
                                           label.x.npc='center', label.y.npc = "top",
                                           comparisons = list(c("primitive", "committed")))
print(plot.ki+theme(legend.position="none")+xlab("subtype")+ ggtitle("KI"))



BeatAML = readRDS("../data/PERT_Score_BeatAML.rds")
plot.beatAML = ggplot(BeatAML, aes(x=cluster, y=score, fill=cluster))
plot.beatAML = plot.beatAML + geom_violin(trim=F, color=NA, scale="width")
plot.beatAML = plot.beatAML + stat_summary(fun.data=data_summary, geom="pointrange", color="white")
plot.beatAML = plot.beatAML + scale_fill_manual(values=colClust) + theme_classic()
plot.beatAML = plot.beatAML + stat_compare_means(method = "wilcox.test", label = "p.signif",
                                       label.x.npc='center', label.y.npc = "top",
                                       comparisons = list(c("primitive", "committed")))
print(plot.beatAML+theme(legend.position="none")+xlab("subtype")+ ggtitle("BeatAML"))


Leucegene = readRDS("../data/PERT_Score_Leucegene.rds")
plot.Leucegene = ggplot(Leucegene, aes(x=cluster, y=score, fill=cluster))
plot.Leucegene = plot.Leucegene + geom_violin(trim=F, color=NA, scale="width")
plot.Leucegene = plot.Leucegene + stat_summary(fun.data=data_summary, geom="pointrange", color="white")
plot.Leucegene = plot.Leucegene + scale_fill_manual(values=colClust) + theme_classic()
plot.Leucegene = plot.Leucegene + stat_compare_means(method = "wilcox.test", label = "p.signif",
                                                 label.x.npc='center', label.y.npc = "top",
                                                 comparisons = list(c("primitive", "committed")))
print(plot.Leucegene+theme(legend.position="none")+xlab("subtype")+ ggtitle("Leucegene"))


dev.off()
