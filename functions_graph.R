make_numeric_indices <- function(ind, n=NULL, unique=TRUE, sort=TRUE){
  if(is.logical(ind)){
    if(is.null(n)){
      n <- length(ind)
    }
    ind0 <- seq_len(n)
    ind <- ind0[ind]
  }
  ind <- as.integer(ind)
  if(unique){
    ind <- unique(ind)
  }
  if(sort){
    ind <- sort(ind)
  }
  return(ind)
}

getDanubeFlowGraph <- function(stationIndices=NULL, directed=FALSE){
  # If not specified, use all indices
  if(is.null(stationIndices)){
    stationIndices <- TRUE
  }
  # Convert to numerical indices
  stationIndices <- make_numeric_indices(stationIndices, nrow(danube$info))
  
  # Make sure danube data is available
  #danube <- graphicalExtremes::danube
  
  # Make (full) flow graph
  g <- igraph::graph_from_edgelist(danube$flow_edges, directed)
  
  # Keep only specified station indices
  igraph::V(g)$name <- as.character(seq_along(igraph::V(g)))
  g <- igraph::induced_subgraph(g, stationIndices)
  
  return(g)
}

plotDanubeIGraph <- function(
    stationIndices = NULL,
    graph = NULL,
    directed = NULL,
    labelStations = TRUE,
    vertexColors = NULL,
    vertexShapes = NULL,
    edgeColors = NULL,
    edge.width = NULL,
    ...
){
  #danube <- graphicalExtremes::danube
  
  # If not specified, use all indices
  if(is.null(stationIndices)){
    stationIndices <- TRUE
  }
  # Convert to numerical indices
  stationIndices <- make_numeric_indices(stationIndices, nrow(danube$info))
  
  # Handle vertex labelling
  labels <- NA # --> do not label
  if(labelStations){
    labels <- NULL # --> infer labels from graph
  }
  
  if(is.null(graph)){
    if(is.null(directed)){
      directed <- FALSE
    }
    graph <- getDanubeFlowGraph(stationIndices, directed)
  }
  # pos <- as.matrix(danube$info[,c('PlotCoordX', 'PlotCoordY')])
  pos <- as.matrix(danube$info[stationIndices,c('PlotCoordX', 'PlotCoordY')])
  pdf("/home/leej40/Documents/PTC/Code/Network_Graph.pdf", width = 6, height = 6)
  par(mar=c(0,0,0,0))
  igraph::plot.igraph(
    graph,
    layout = pos,
    vertex.label = labels,
    vertex.shape = vertexShapes,
    vertex.color = vertexColors,
    edge.color = edgeColors,
    vertex.frame.color = "black",
    vertex.label.color="black",
    edge.width=edge.width,
    vertex.label.cex = 1.5,
    margin = c(0, 0, 0, 0),
    asp=0,
    ...
  )
  dev.off()
}






