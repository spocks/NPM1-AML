library(igraph)
set.seed(7)

network_data = readRDS("../data/network_data.rds")
anet = graph_from_data_frame(d=network_data$link, vertices=network_data$node, 
                              directed=F)

##--- find communities using the label propagating method 
clp = cluster_label_prop(anet)
l = layout_with_fr(anet, niter = 5000)

##--- plot data -----

pdf("../results/figure-1_A_network.pdf", width = 6, height = 7)
plot(clp, anet, layout=l, edge.color=E(anet)$color, col=V(anet)$color,
     vertex.frame.color	="gray50",
     vertex.label = V(anet)$dataset, vertex.label.color="white")

dev.off()

