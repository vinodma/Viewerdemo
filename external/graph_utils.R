library(plyr)
library(igraph)

getcommunity_from_graphdb<-function(nodeval){
  #nodeval <- input$comm_id
  matchexp <- paste0("{\"", nodeval,"\":",sep="")
  con <- file("./www/data/database.json")
  writesubsr <- NULL
  open(con)
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    #print(grepl(matchexp,line,fixed=TRUE))
    if((grepl(matchexp,line,fixed=TRUE))==TRUE){
      f <-"./www/data/current_graph.json"
      pos = regexpr(matchexp, line,fixed=TRUE)
      writesubsr = substr(line,attr(pos,"match.length") + 1,nchar(line) -1)
      Sys.chmod(f, (file.info(f)$mode | "664"))
      file.remove(f)
      sink(file=f)
      cat(writesubsr)
      sink()
      
      
    }
  }
  close(con)
  
}


build_initial_graph <-function (file){
  # Uses igraph to parse the initial, full graph
  table <- read.csv(file, header = TRUE, sep = ",")
  edges <- table[c('entity1', 'entity2')]
  data1 <- table[c('entity1', 'type1')]
  data2 <- table[c('entity2', 'type2')]
  colnames(data1) <- c('entity', 'type')
  colnames(data2) <- c('entity', 'type')
  vertex_data <- unique(rbind(data1, data2))
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = vertex_data)
  return(g)
}


get_communities <- function(graph){
  # Runs louvain community detection
  return(multilevel.community(graph))
}

get_community_graph <- function(graph, communities){
  # Builds a graph of the communities
  V(graph)$comm <-communities$membership
  contracted <- contract.vertices(graph, communities$membership, "random")
  community_graph <- simplify(contracted, "random")
  V(community_graph)$name <- V(community_graph)$comm
  
  # Set the size of each node to be proportional to the community size
  counts <- count(communities$membership)
  V(community_graph)$size <- counts$freq
  
  
  return(community_graph)
}

subgraph_of_one_community <- function(graph, communities, community_id){
  # Builds a subgraph of one community from the original graph
  idx <- which(communities$membership == community_id)
  subgraph <- induced.subgraph(graph, idx)
  return(subgraph)
}