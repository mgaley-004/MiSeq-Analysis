library(igraph)
library(phyloseq)
library(ggplot2)
library(ggnetwork)
library(phyloseq)

graph_perm_test = function(physeq, sampletype, grouping = 1:nsamples(physeq),
                           distance, type = c("mst", "knn", "threshold.value", "threshold.nedges"),
                           max.dist = .4, knn = 1, nedges = nsamples(physeq), keep.isolates = TRUE, nperm = 499) {
  type = match.arg(type)
  # make the network
  d = distance
  if(!validGrouping(sample_data(physeq), sampletype, grouping)) {
    stop("Not a valid grouping, all values of sampletype must
              be the same within each level of grouping")
  }
  switch(type,
         "threshold.value" = {
           neighbors = as.matrix(d) <= max.dist
           diag(neighbors) = 0
           net = graph.adjacency(neighbors, mode = "undirected", add.colnames = "name")   
         },
         "threshold.nedges" = {
           threshold = sort(as.vector(d))[nedges]
           neighbors = as.matrix(d) <= threshold
           diag(neighbors) = 0
           net = graph.adjacency(neighbors, mode = "undirected", add.colnames = "name")
         },
         "knn" = {
           neighbors = t(apply(as.matrix(d),1, function(x) {
             r = rank(x)
             nvec = ((r > 1) & (r < (knn + 2))) + 0
           }))
           neighbors = neighbors + t(neighbors)
           net = graph.adjacency(neighbors, mode = "undirected",
                                 add.colnames = "name", weighted = TRUE)
         },
         "mst" = {
           gr = graph.adjacency(as.matrix(d), mode = "undirected", weighted = TRUE,
                                add.colnames = "name")
           net = minimum.spanning.tree(gr, algorithm = "prim")
         }           
  )
  el = get.edgelist(net)
  sampledata = data.frame(sample_data(physeq))
  elTypes = el
  elTypes[,1] = sampledata[el[,1], sampletype]
  elTypes[,2] = sampledata[el[,2], sampletype]
  observedPureEdges = apply(elTypes, 1, function(x) x[1] == x[2])
  edgeType = sapply(observedPureEdges, function(x) if(x) "pure" else "mixed")
  # set these attributes for plotting later
  if(is.factor(sampledata[,sampletype]))
    V(net)$sampletype = as.character(sampledata[,sampletype])
  else
    V(net)$sampletype = sampledata[,sampletype]
  E(net)$edgetype = edgeType
  
  # find the number of pure edges for the non-permuted data
  nobserved = sum(observedPureEdges)
  origSampleData = sampledata[,sampletype]
  names(origSampleData) = rownames(sampledata)
  # find the permutation distribution of the number of pure edges
  permvec = numeric(nperm)
  for(i in 1:nperm) {
    sampledata[,sampletype] = permute(sampledata, grouping, sampletype)
    elTypes = el
    elTypes[,1] = sampledata[el[,1], sampletype]
    elTypes[,2] = sampledata[el[,2], sampletype]
    permPureEdges = apply(elTypes, 1, function(x) x[1] == x[2])
    permvec[i] = sum(permPureEdges)
  }
  pval = (sum(permvec >= nobserved) + 1) / (nperm + 1)
  if(!keep.isolates) {
    degrees = igraph::degree(net)
    net = igraph::induced_subgraph(net, which(degrees > 0))
  }
  out = list(observed = nobserved, perm = permvec, pval = pval,
             net = net, sampletype = origSampleData, type = type)
  class(out) = "psgraphtest"
  return(out)
}

print.psgraphtest <- function(x, ...) {
  cat("Output from graph_perm_test\n")
  cat("---------------------------\n")
  cat(paste("Observed test statistic: ", x$observed, " pure edges", "\n", sep = ""))
  cat(paste(nrow(get.edgelist(x$net)), " total edges in the graph", "\n", sep = ""))
  cat(paste("Permutation p-value: ", x$pval, "\n", sep = ""))
}


permute = function(sampledata, grouping, sampletype) {
  if(length(grouping) != nrow(sampledata)) {
    grouping = sampledata[,grouping]
  }
  x = as.character(sampledata[,sampletype])
  # gives the original mapping between grouping variables and sampletype
  labels = tapply(x, grouping, function(x) x[1])
  # permute the labels of the groupings
  names(labels) = sample(names(labels))
  return(labels[as.character(grouping)])
}


validGrouping = function(sd, sampletype, grouping) {
  if(!(sampletype %in% colnames(sd))) {
    stop("\'sampletype\' must be a column names of the sample data")
  }
  if(!(grouping %in% colnames(sd)) && (length(grouping) != nrow(sd))) {
    stop("\'grouping\' must be either a column name of the sample data
             or a vector with number of elements equal to the number of samples")
  }
  sd = data.frame(sd)
  if(length(grouping) != nrow(sd)) {
    grouping = sd[,grouping]
  }
  valid = all(tapply(sd[,sampletype], grouping, FUN = function(x)
    length(unique(x)) == 1))
  return(valid)
}


plot_test_network = function(graphtest, scaleshapes, scalecolors, layoutstring) {
  if(graphtest$type == "mst")
    layoutMethod = "kamadakawai"
  else
    layoutMethod = "fruchtermanreingold"
  #ggplot(graphtest$net,
   #      aes_string(x = "x", y = "y", xend = "xend", yend = "yend"), layout = layoutMethod) +
   # geom_edges(aes_string(linetype = "edgetype")) +
  #  geom_nodes(aes_string(color = "sampletype", size=0.75)) +
   # scale_linetype_manual(values = c(3,1)) + theme_blank()
  gnplot <- ggnetwork(graphtest$net)
  ggplot(gnplot,
         aes_string(x = "x", y = "y", xend = "xend", yend = "yend"), layout = layoutstring) +
    geom_edges(aes_string(linetype = "edgetype")) +
    geom_nodes(aes_string(color = "sampletype", size=0.6, shape="sampletype", stroke=1)) +
    scale_linetype_manual(values = c(3,1)) + theme_blank() + scale_shape_manual(values=scaleshapes) + scale_color_manual(values=scalecolors)
}

plot_permutations = function(graphtest, bins = 30) {
  p = qplot(graphtest$perm, geom = "histogram", bins = bins)
  if(packageVersion("ggplot2") >= "2.2.1.9000") {
    ymax = max(ggplot_build(p)$layout$panel_scales_y[[1]]$get_limits())
  } else {
    ymax = ggplot_build(p)$layout$panel_ranges[[1]][["y.range"]][2]
  }
  p + geom_segment(aes(x = graphtest$observed, y = 0,
                       xend = graphtest$observed, yend = ymax / 10), color = "red") +
    geom_point(aes(x = graphtest$observed, y = ymax / 10), color = "red") +
    xlab("Number of pure edges")                         
}
