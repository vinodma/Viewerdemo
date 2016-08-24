library(DT)
library(shiny)
library(igraph)
library(plotly)
library(rstackdeque)

source("external/graph_utils.R", local = TRUE)
source("external/makenetjson.R", local = TRUE)
source("external/protein_label_dictionary.R",local = TRUE)

#initial_data <- "./www/data/ctd.csv"
#graph <- build_initial_graph(initial_data)
#communities <- get_communities(graph)
#htmlloaded = FALSE
#s1 <- rstack()
s2 <-rstack()
s3 <- rstack()
mp <<- NULL
sortedlabel<-NULL
protienDSpathway<<-data.frame()

function(input, output, session){ 
  global <- reactiveValues()
  global$is_comm_graph = TRUE
  global$currentCommId <- -1
  #global$viz_stack <- insert_top(s1, list(graph, communities))
  global$name <- insert_top(s2, "")
  global$commID <- insert_top(s3, -1)
  
  # reset button
  observeEvent(input$reset_button, {
    getcommunity_from_graphdb(-1)
    observe({
      session$sendCustomMessage(type = "updategraph",message="clear")
    })
    #global$viz_stack <- rstack()
    #global$viz_stack <- insert_top(global$viz_stack, list(graph, communities))
    #global$name <- insert_top(s2, "")
  })
  
  observeEvent(input$variable, {
    print(input$variable)
  })
  
  #Search button
  observeEvent(input$search_button,{
    searchelm <- input$searchentitiy
    lbllist <<- c()
    withProgress(message = "Searching ...",value = 0,{
      if(grepl(",",searchelm) == FALSE){
        getallparentforentity(searchelm)
      }
      else
      {
        res<-unlist(strsplit(searchelm,","))
        lapply(res,getallparentforentity)
        
      }
    })
    lbllist <- unique(lbllist)
    memcommunity<-paste(lbllist,collapse=",")
    memcommunity<-paste(searchelm,memcommunity,sep=",")
    #memcommunity=input$searchentitiy
    observe({
      session$sendCustomMessage(type = "commmemmsg" ,
                                message = list(id=memcommunity))
    })
  })
  
  # table click
  observe({
    row <- input$degree_table_rows_selected
    if (length(row)){
      print(row)
      session$sendCustomMessage(type = "commmemmsg" ,
                                message = list(id=tail(row, n=1)))
    }
  })
  
  # disease pathway table click
  observe({
    row <- input$plotgraph1_rows_selected
    last_selected_row = tail(row, n=1)
    print(last_selected_row)
    if( (!is.null(row)) && (length(row)>=1)){
      proteins<-protienDSpathway[protienDSpathway$Pathway==unlist(last_selected_row),]$Protein
      session$sendCustomMessage(type = "commmemmsg" ,
                                message = list(id=paste(proteins,collapse=",")))
    }
  })
  
  
  
  # back button
  observeEvent(input$back_button, {
    size <- length(global$viz_stack)
    if (size > 1){
      global$viz_stack <- without_top(global$viz_stack)
      global$name <- without_top(global$name)
    } 
  })
  
  # on-click from sigma.js
  observeEvent(input$comm_id, {
    print(input$comm_id)
    global$currentCommId<-input$comm_id 
    getcommunity_from_graphdb(input$comm_id)
    update_stats()
    observe({
      session$sendCustomMessage(type = "updategraph",message="xyz")
    })
  })
  
  
  
  # render with sigma the current graph (in json)
  output$graph_with_sigma <- renderUI({
    getcommunity_from_graphdb(-1)
    #makenetjson(data[[1]], "./www/data/current_graph.json", data[[2]]) 
    update_stats()
    
    observe({
      session$sendCustomMessage(type = "updategraph",message="")
    })
    
    return(includeHTML("./www/graph.html"))
  })
  
 # update the summary stats
  update_stats <- function(){
    con <- file("./www/data/current_graph.json")
    open(con)
    line <- readLines(con, n = 1, warn = FALSE)
    close(con)
    x<-fromJSON(line)
    edges<-x$edges[c('source','target')]
    vertex_data<-x$nodes[c('id','name','type')]
    if(nrow(vertex_data) > 1){
      graph <- graph_from_data_frame(edges, directed = FALSE, vertices =
vertex_data)
    
    nodes <- get.data.frame(graph, what="vertices")
    nodes$degree <- degree(graph)
    nodes$pagerank <- page_rank(graph)$vector
    colnames(nodes) <- c("Name", "Type", "Degree", "PageRank")

    global$nodes <- nodes
    }
    else
    {
      global$nodes <- NULL
    }
  }
  
  # Plot the degree distribution of the current graph
  output$degree_distribution <- renderPlotly({  
    if (!is.null(global$nodes)){
      plot_ly(global$nodes, x = Degree, type="histogram",  color="#FF8800")
    }
  })
  
  # Plot the pagerank distribution of the current graph
  output$pagerank_distribution <- renderPlotly({
    if (!is.null(global$nodes)){
      plot_ly(global$nodes, x = PageRank, type="histogram", color="#FF8800")
    }    
  })
  
  # Generate a table of node degrees
  output$degree_table <- DT::renderDataTable({
    if (!is.null(global$nodes)){
      table <- global$nodes[c("Name", "Degree", "PageRank")]
    }
  },
  options = list(order = list(list(1, 'desc'))),
  rownames = FALSE,
  selection = "single"
  )
  
  # Generate the current graph name (as a list of community labels)
  output$name <- renderText({
    name <- as.list(rev(global$name))
    name <- paste(name, collapse = "/", sep="/")
    #return(paste(c("Current Community", name)))
  })
  
  output$plotgraph1 <- DT::renderDataTable({
    protienDSpathway<<-data.frame()
    sortedlabel<-NULL
    lf<-NULL
    lbls<-NULL
    
    # This takes forever. If we can load a previously built object do it; otherwise don't hold your breath
    withProgress(message = "Loading ...",value = 0,{
      if(is.null(mp)){
        filename = 'mp.rds'
        if (file.exists(filename)){
          mp <<- NULL
          mp <<- readRDS(filename)
        } else {
          mp <<- getproteinlabeldict()
          saveRDS(mp, file=filename)
        }
      }
    })
    if(global$currentCommId==-1)
      return (NULL)
    finallist<-c()
    lbllist <<- c()
    
    withProgress(message = "Loading ...",value = 0,{
      getrawentititesfromComm(global$currentCommId)
    })
    
    table <- data.frame(Protein="No pathway data available")

    if (nrow(protienDSpathway)>1){
      labelfreq <- table(protienDSpathway)
      if (ncol(labelfreq)>1){
        z<-apply(labelfreq,1,sum)
        sortedlabel<-labelfreq[order(as.numeric(z), decreasing=TRUE),]
        table<-as.data.frame.matrix(sortedlabel)
      } else {
        table <- as.data.frame.matrix(labelfreq)
      }
      row.names(table) <- strtrim(row.names(table), 50)
    } 
    table
  },
  rownames = TRUE,
  selection = "single")
  
  
  
}
