  #Creates a flowchart. A stratified flowchart is presented unless otherwise specified
  create_flowchart <- function(model, type = 'stratified') {
    #Pick type of input into the function, depending on whether the type of flowchart is 
    if (type == 'stratified') {
      input_nodes <- names(model$compartment_values)
      type_of_flow <- model$flows
    } 
    else if (type == 'unstratified') {
      input_nodes <- model$compartment_types
      type_of_flow <- model$unstratified_flows
    }
    else {
      stop("Type needs to be either stratified or unstratified.")
    }
    #The inputs for the flowchart is ordered alphabetically
    sorted_nodes <- sort(input_nodes)
    #Inputs are collapsed into one string so it can be put into grViz
    nodes <- paste(sorted_nodes, collapse = ' ')
    #The pathways between nodes are set as empty
    connection_between_nodes <- ""
    #The pathway between nodes is populated from type_of_flow
    for (row in seq(nrow(type_of_flow))) {
      if (type_of_flow$implement[[row]]) {
        connection_between_nodes <- paste(connection_between_nodes, type_of_flow$from[[row]], "->", type_of_flow$to[[row]], ' ', sep = '')
      }}
    #The string necessary for grViz is created here
    input_for_grViz <- paste("digraph dot {
                   graph [layout = dot,
                   rankdir = LR] 
                   node [shape = box,
                   fontname = Helvetica, style = filled, color = red]
                  ",
                   nodes,
                   connection_between_nodes, '}')
    #~ are substituted for _ and input into the function
    grViz(str_replace_all(input_for_grViz, "~", "_"))
  }

  
  
  
  #Creates a flowchart. A stratified flowchart is presented unless otherwise specified
  create_flowchart <- function(model, type = 'stratified') {
    #Pick type of input into the function, depending on whether the type of flowchart is 
    if (type == 'stratified') {
      input_nodes <- names(sir_model$compartment_values)
      type_of_flow <- model$flows
    } 
    else if (type == 'unstratified') {
      input_nodes <- model$compartment_types
      type_of_flow <- model$unstratified_flows
    }
    else {
      stop("Type needs to be either stratified or unstratified.")
    }
    #The inputs for the flowchart is ordered alphabetically
    sorted_nodes <- sort(input_nodes)
    #Inputs are collapsed into one string so it can be put into grViz
    broken_down_nodes <- list()
    for (stem_value in 1:length(model$compartment_types)) {
      x_vector <- c()
      for (stem_type in 1:length(sorted_nodes)) {
        if (model$compartment_types[[stem_value]] == find_stem(sorted_nodes[[stem_type]])) {
          x_vector <- c(x_vector, sorted_nodes[[stem_type]])
        }
      }
      broken_down_nodes[[stem_value]] <- x_vector
    }
    settings <- ''
    for (list_different_nodes in 1:length(broken_down_nodes)) {
      nodes <- c()
      nodes <- paste(broken_down_nodes[[list_different_nodes]], collapse = ' ')
      settings <- paste(settings, 'node [shape = box,
      fontname = Helvetica, style = filled, color =', 
      c('BlanchedAlmond', 'Grey', 
        'RosyBrown', 'LavenderBlush',
        'Salmon', 'LightPink', 
        'PaleGreen', 'Thistle', 
        'Beige', 'PeachPuff', 
        'MintCream', 'AquaMarine', 
        'MistyRose', 'Tomato',
        'Honeydew', 'LightCyan')[[sample(1:16, 1, 
                                              replace = FALSE, 
                                              prob = NULL)]],
      ']', nodes)
    }
    #The pathways between nodes are set as empty
    connection_between_nodes <- ""
    #The pathway between nodes is populated from type_of_flow
    for (row in seq(nrow(type_of_flow))) {
      if (type_of_flow$implement[[row]]) {
        connection_between_nodes <- paste(connection_between_nodes, type_of_flow$from[[row]], "->", type_of_flow$to[[row]], ' ', sep = '')
      }}
    #The string necessary for grViz is created here
    input_for_grViz <- paste("digraph dot {
                             graph [layout = dot,
                             rankdir = LR]", 
                             settings,
                             connection_between_nodes, '}')
    #~ are substituted for _ and input into the function
    grViz(str_replace_all(input_for_grViz, "~", "_"))
  }
  
  
  
  
    
  
  create_flowchart(sir_model, type = 'unstratified')





Create_unstratified_flowchart <- function(model) {
  sorted_compartment_types_flowchart <- sort(model$compartment_types)
  nodes <- paste(sorted_compartment_types_flowchart, collapse = '; ')
  final <- ""
  for (row in seq(nrow(model$unstratified_flows))) {
    if (model$unstratified_flows$implement[[row]]) {
      final <- paste(final, model$unstratified_flows$from[[row]], "->", model$unstratified_flows$to[[row]], ' ', sep = '')
    }}
  hello <- paste("digraph dot {
                 
                 graph [layout = dot,
                 rankdir = LR] 
                 # several 'node' statements
                 node [shape = box,
                 fontname = Helvetica] \n",
                 nodes,
                 "\n",
                 final, '}')
  grViz(str_replace_all(hello, "~", "_"))
}


Create_stratified_flowchart <- function(model) {
  sorted_compartment_types_flowchart <- sort(names(model$compartment_values))
  nodes <- paste(sorted_compartment_types_flowchart, collapse = '; ')
  final <- ""
  for (row in seq(nrow(model$flows))) {
    if (model$flows$implement[[row]]) {
      final <- paste(final, model$flows$from[[row]], "->", model$flows$to[[row]], ' ', sep = '')
    }}
  hello <- paste("digraph dot {
                 graph [layout = dot,
                 rankdir = LR] 
                 # several 'node' statements
                 node [shape = box,
                 fontname = Helvetica] \n",
                 nodes,
                 "\n",
                 final, '}')
  grViz(str_replace_all(hello, "~", "_"))
}

M <- matrix(nrow = length(self$compartment_types), 
            ncol = length(self$compartment_types), 
            byrow = TRUE, data = 0)
for row in seq(nrow(self$flows)) {
  if (self$flows$implement) {
    M[self$flows$to, self$flows$from] <- self$flows$parameter
  }
}
plotmat(M, curve = 0, name = compartment_types, lwd = 1, box.lwd = 2, 
        cex.txt = 0.8, box.size = 0.1, box.type = "square", box.prop = 0.5)




M <- matrix(nrow = length(sorted_compartment_types_flowchart), 
            ncol = length(sorted_compartment_types_flowchart), 
            byrow = TRUE, data = 0)
dimnames(M) <- list(sorted_compartment_types_flowchart, 
                    sorted_compartment_types_flowchart)
for (row in seq(nrow(sir_model$flows))) {
  if (sir_model$flows$implement[[row]]) {
    arrow_output <- sir_model$flows$parameter[[row]]
    M[sir_model$flows$to[[row]], sir_model$flows$from[[row]]] <- arrow_output
  }
}

sorted_compartment_types_flowchart
paste(sorted_compartment_types_flowchart, collapse = '; ')

plotmat(M, curve = 0, name = sorted_compartment_types_flowchart, lwd = 1, box.lwd = 2, 
        cex.txt = 0.8, box.size = 0.15, box.type = "square", box.prop = 0.2)


                     
                      
                    


final <- ""
for (row in seq(nrow(sir_model$flows))) {
  if (sir_model$flows$implement[[row]]) {
    final <- paste(final, sir_model$flows$from[[row]], "->", sir_model$flows$to[[row]], ' ', sep = '')
  }}
hello <- paste("digraph boxes_and_circles {
                # a 'graph' statement
               graph [overlap = true, fontsize = 10]
               # several 'node' statements
               node [shape = box,
               fontname = Helvetica] \n",
               nodes,
               "\n",
               final, '}')
grViz(str_replace_all(hello, "~", "_"))




stem_names_for_flowchart <- c()
for (names in names(sir_model$compartment_values)) {
  stem_names_for_flowchart <- c(stem_names_for_flowchart, find_stem(names))
}  
sorted_compartment_types_flowchart <- sort(names(sir_model2$compartment_values))
nodes <- paste(sorted_compartment_types_flowchart, collapse = '; ')
final <- ""
for (row in seq(nrow(sir_model2$flows))) {
  if (sir_model2$flows$implement[[row]]) {
    final <- paste(final, sir_model2$flows$from[[row]], "->", sir_model2$flows$to[[row]], ' ', sep = '')
  }}
hello <- paste("digraph dot {

               graph [layout = dot,
               rankdir = LR] 
               # several 'node' statements
               node [shape = box,
               fontname = Helvetica] \n",
               nodes,
               "\n",
               final, '}')
grViz(str_replace_all(hello, "~", "_"))



Create_unstratified_flowchart <- function(model) {
  sorted_compartment_types_flowchart <- sort(model$compartment_types)
  nodes <- paste(sorted_compartment_types_flowchart, collapse = '; ')
  final <- ""
  for (row in seq(nrow(model$unstratified_flows))) {
    if (model$unstratified_flows$implement[[row]]) {
      final <- paste(final, model$unstratified_flows$from[[row]], "->", model$unstratified_flows$to[[row]], ' ', sep = '')
    }}
  hello <- paste("digraph dot {
                 
                 graph [layout = dot,
                 rankdir = LR] 
                 # several 'node' statements
                 node [shape = box,
                 fontname = Helvetica] \n",
                 nodes,
                 "\n",
                 final, '}')
  grViz(str_replace_all(hello, "~", "_"))
}


Create_stratified_flowchart <- function(model) {
  sorted_compartment_types_flowchart <- sort(names(model$compartment_values))
  nodes <- paste(sorted_compartment_types_flowchart, collapse = '; ')
  final <- ""
  for (row in seq(nrow(model$flows))) {
    if (model$flows$implement[[row]]) {
      final <- paste(final, model$flows$from[[row]], "->", model$flows$to[[row]], ' ', sep = '')
    }}
  hello <- paste("digraph dot {
                 
                 graph [layout = dot,
                 rankdir = LR] 
                 # several 'node' statements
                 node [shape = box,
                 fontname = Helvetica] \n",
                 nodes,
                 "\n",
                 final, '}')
  grViz(str_replace_all(hello, "~", "_"))
}

Create_stratified_flowchart(sir_model2)




Create_unstratified_flowchart(sir_model2)




stem_names_for_flowchart <- c()
for (names in names(sir_model2$compartment_values)) {
  stem_names_for_flowchart <- c(stem_names_for_flowchart, find_stem(names))
}  







  


    
stem_names_for_flowchart <- c()
for (names in names(sir_model$compartment_values)) {
  stem_names_for_flowchart <- c(stem_names_for_flowchart, find_stem(names))
}   
nodes <- paste(sorted_compartment_types_flowchart, collapse = '; ')  
for (row in seq(nrow(sir_model$flows))) {
  if (sir_model$flows$implement[[row]]) {
    sir_model$flows$from[[row]] <- sir_model$flows$to[[row]]
  }}
paste("
digraph boxes_and_circles {
      # a 'graph' statement
      graph [overlap = true, fontsize = 10]
      # several 'node' statements
      node [shape = box,
      fontname = Helvetica] \n",
      nodes,
      "\n # several 'edge' statements
      for (row in seq(nrow(sir_model$flows))) {
      if (sir_model$flows$implement[[row]]) {
      arrow_output <- sir_model$flows$parameter[[row]]
      M[sir_model$flows$to[[row]], sir_model$flows$from[[row]]] <- arrow_output
      })")



paste("digraph boxes_and_circles {
      
                # a 'graph' statement
                graph [overlap = true, fontsize = 10]
                
                # several 'node' statements
                node [shape = box,
                fontname = Helvetica] \n",
      nodes,
      "# several 'edge' statements
                for (row in seq(nrow(sir_model$flows))) {
                if (sir_model$flows$implement[[row]]) {
                arrow_output <- sir_model$flows$parameter[[row]]
                M[sir_model$flows$to[[row]], sir_model$flows$from[[row]]] <- arrow_output
                })")


cat(paste('h', '\nh'))
str_c('hello',nodes, sep = "/n")
'h'
'd
d'
count(stem_names_for_flowchart, vars = "id")
stem_names_for_flowchart

M[sir_model$flows$to[[1]], sir_model$flows$from[[1]]]


M <- matrix(nrow = length(sir_model$compartment_types), ncol = length(sir_model$compartment_types), byrow = TRUE, data = 0)


self$flows$from
sir_model$flows$from

c(sir_model$flows$from, sir_model$flows$to)




plotmat()

grViz("
digraph boxes_and_circles {
      
      # a 'graph' statement
      graph [overlap = true, fontsize = 10]
      
      # several 'node' statements
      node [shape = box,
      fontname = Helvetica]
      A; B; C; D; E; F
      
      node [shape = circle,
      fixedsize = true,
      width = 0.9] // sets as circles
      1; 2; 3; 4; 5; 6; 7; 8
      
      # several 'edge' statements
      A->1 B->2 B->3 B->4 C->A
      1->D E->A 2->4 1->5 1->F
      E->6 4->6 5->7 6->7 3->8
      }
      ")







M <- matrix(nrow = length(sir_model$compartment_types), 
            ncol = length(sir_model$compartment_types), 
            byrow = TRUE, data = 0)
dimnames(M) <- list(sir_model$compartment_types, sir_model$compartment_types)
for (row in seq(nrow(sir_model$flows))) {
  if (sir_model$flows$implement[[row]]) {
    arrow_output <- paste(sir_model$flows$type[[row]], sir_model$flows$parameter[[row]], sep = '_')
    M[sir_model$flows$to[[row]], sir_model$flows$from[[row]]] <- arrow_output
  }
}
plotmat(M, curve = 0, name = sir_model$compartment_types, lwd = 1, box.lwd = 2, 
        cex.txt = 0.8, box.size = 0.15, box.type = "square", box.prop = 0.2)






stem_names_for_flowchart

names(sir_model$compartment_values)

sort(names(sir_model$compartment_values))
