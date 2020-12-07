library(ConsensusClusterPlus)
SEED = 1262118388.71279

get_silhouette <- function(df, cc)
{
  cord <- Hmisc::rcorr(t(df), type = c("pearson", "spearman")[2])
  cord <- as.dist(1 - cord$r)
  samplOrd <- rownames(as.matrix(cord))
  SS <- list()
  for (i in 2:8)
  {
    cl <- cc$results[[i]]$consensusClass[samplOrd]
    SS[[paste0("K", i)]] <- cluster::silhouette(cl, cord)
  }
  sw <- sapply(SS, function(x)
    mean(x[, 3]))
  return(sw)
}

ConsensusClusterPlus_modified <-
  function(data,
           maxK,
           reps,
           pItem,
           pFeature,
           clusterAlg,
           distance,
           outFileName,
           seed = SEED,
           innerLinkage = "average",
           finalLinkage = "average")
  {
    createRandString<- function() {
      v = c(sample(LETTERS, 3, replace = TRUE),
            sample(0:9, 3, replace = TRUE),
            sample(letters, 3, replace = TRUE))
      return(paste0(sample(v),collapse = ""))
    }
    px1 = createRandString()
    pathX = sprintf("ConClu_TMP_%s", px1)
    dir.create(pathX,
               showWarnings = TRUE,
               recursive = FALSE,
               mode = "0777")
    
    results = ConsensusClusterPlus(
      data,
      maxK = maxK,
      reps = reps,
      pItem = pItem,
      pFeature = pFeature,
      title = pathX,
      clusterAlg = clusterAlg,
      distance = distance,
      seed = seed,
      plot = "pdf",
      innerLinkage = innerLinkage,
      finalLinkage = finalLinkage
    )
    
    icl = calcICL(results, title = pathX, plot = "pdf")
    
    cmd1 = sprintf("cp %s/consensus.pdf %s_consensus.pdf", pathX, outFileName)
    system(cmd1)
    cmd2 = sprintf("cp %s/icl.pdf %s_icl.pdf", pathX, outFileName)
    system(cmd2)
    cmd3 = sprintf("rm -rf %s", pathX)
    system(cmd3)
    return(list(results = results, icl = icl))
  }

heirarchical_consensClust <-
  function(dfDt,
           outFileName,
           distance = "spearman",
           pItem = 0.8,
           pFeature = 0.8,
           reps = 1000,
           maxK = 4,
           innerLinkage = "ward.D2",
           finalLinkage = "average")
  {
    dataX = as.matrix(t(dfDt))
    consRes = ConsensusClusterPlus_modified(
      data = dataX,
      maxK = maxK,
      reps = reps,
      pItem = pItem,
      pFeature = pFeature,
      clusterAlg = "hc",
      distance = distance,
      outFileName = outFileName,
      seed = SEED,
      innerLinkage = innerLinkage,
      finalLinkage = finalLinkage
    )
    outFl = sprintf("%s_ClObj.Rda", outFileName)
    saveRDS(consRes, file = outFl)
    get_silhouette(dfDt, consRes)
  }
