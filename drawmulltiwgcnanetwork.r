getLevel <- function (level, design) 
{
  if (level == 1) {
    return("combined")
  }
  return(unique(design[, level]))
}


drawMultiWGCNAnetwork2 <- function (WGCNAlist, comparisonList, moduleOfInterest, design, 
                                    overlapCutoff = 0, padjCutoff = 1, removeOutliers = TRUE, 
                                    alpha = 0.00000000000000000000000000000000000000000000000001, 
                                    layout = NULL, hjust = 0.4, vjust = 0.3, width = 0.5, colors = NULL) 
{
  overlapList = lapply(comparisonList, function(x) x$overlap)
  overlapList = do.call(rbind, overlapList)
  filteredOverlapList = overlapList
  if (removeOutliers) {
    for (WGCNA in WGCNAlist) {
      filteredOverlapList = filteredOverlapList[!filteredOverlapList$mod1 %in% 
                                                  WGCNA@outlierModules, ]
      filteredOverlapList = filteredOverlapList[!filteredOverlapList$mod2 %in% 
                                                  WGCNA@outlierModules, ]
    }
  }
  admittedModules = unique(c(filteredOverlapList$mod1, filteredOverlapList$mod2))
  if (is.null(layout)) {
    myCoords = list()
    for (level in 1:3) {
      WGCNAs = getLevel(level, design)
      from = 0 - width * length(WGCNAs)/2
      to = 0 + width * length(WGCNAs)/2
      x.coordinates = seq(from, to, length.out = length(WGCNAs))
      if (level == 1) 
        x.coordinates = 0
      for (nWGCNA in 1:length(WGCNAs)) {
        nModules = length(admittedModules[startsWith(admittedModules, 
                                                     WGCNAs[[nWGCNA]])])
        myCoords = append(myCoords, list(cbind(runif(nModules, 
                                                     x.coordinates[[nWGCNA]], x.coordinates[[nWGCNA]] + 
                                                       hjust), 3 - level + runif(nModules, -vjust, 
                                                                                 vjust))))
      }
    }
    layout = do.call(rbind, myCoords)
  }
  graph = graph_from_data_frame(d = filteredOverlapList, directed = FALSE)
  vcol = str_split_fixed(V(graph)$name, "_", 2)[, 1]
  conditions = unique(vcol)
  if (is.null(colors)) 
    palette = colors(length(conditions), random = TRUE)
  if (!is.null(colors)) 
    palette = colors
  for (condition in 1:length(conditions)) {
    vcol[vcol == conditions[[condition]]] = palette[[condition]]
  }
  V(graph)$color = vcol
  E(graph)$weight = -log10(E(graph)$p.adj)
  E(graph)$width = scales::rescale(E(graph)$weight, from = c(0, 320), 
                                   to = c(0, 5))
  ealpha = scales::rescale(-log10(E(graph)$p.adj), from = c(0, 320), 
                           to = c(0, 1))
  ecol = lapply(ealpha, function(x) rgb(1, 0, 0, x))
  E(graph)$color = unlist(ecol)
  modulesOfInterest = unique(c(filteredOverlapList$mod2[filteredOverlapList$mod1 == 
                                                          moduleOfInterest & filteredOverlapList$p.adj < alpha], 
                               filteredOverlapList$mod1[filteredOverlapList$mod2 == 
                                                          moduleOfInterest & filteredOverlapList$p.adj < alpha]))
  graph <- delete.edges(graph, which(!(filteredOverlapList$mod1 %in% 
                                         modulesOfInterest | filteredOverlapList$mod2 %in% modulesOfInterest)))
  par(mar = c(0, 0, 0, 0))
  plot = plot(graph,
              vertex.label.color = "black", 
              vertex.size = 5, 
              vertex.label = NA, 
              vertex.label.cex = 0.3, 
              layout = layout,
              rescale = F,asp = 0,xlim=range(layout[,1]),ylim = range(layout[,2]),vertex.label.dist= 1)
  legend("bottomleft",
         legend = conditions,
         pch    = 21,
         pt.bg  = palette,
         col    = palette,
         pt.cex = 1,
         cex    = 0.4,           
         bty    = "n")
  return(plot)
}
