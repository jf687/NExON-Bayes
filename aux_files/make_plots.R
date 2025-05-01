# Functions to create network plots

create_layout <- function(network){
  layout_coords <- sna::gplot.layout.fruchtermanreingold(
    network::as.matrix.network.adjacency(
      network::network(network > 0.5, directed = FALSE)
    ),
    layout.par = list()
  )
  return(layout_coords)
}



create_network_plots <- function(data_list, layout_coords, color = "red",edge_color = "purple",
                                 title = "Network Plots") {
  N <- length(data_list)
  networks <- vector("list", N)
  plots <- vector("list", N)

  
  # Loop to generate networks and their corresponding plots
  for (i in seq_len(N)) {
    networks[[i]] <- network::network(data_list[[i]] > 0.5, directed = FALSE)
    plots[[i]] <- GGally::ggnet2(
      networks[[i]], edge.size = 0.3, alpha = 0.8, mode = layout_coords, 
      color = color, label = F, size = 1, size.min = 0, edge.color = edge_color, title = "test"
    ) +
      title(main = title)#+
      #ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 6))

  }
  
  # Print all plots using gridExtra
  print(gridExtra::grid.arrange(grobs = plots, ncol = N), title = "test")
}
