#' Make a bar chart from a named numeric vector or table
#'
#' @param vec A named numeric vector
#' @param vec_levels The order of the names in the plot
#' @param xlab The x-axis title
#' @param ylab The y-axis title
#' @param title The title of the plot
#' @param x_labels_angle  The angle to rotate the x-axis labels. 
#' @return A ggplot object
#' @export
ggplot_col_from_vector = function(vec, vec_levels = NULL, xlab = NULL, ylab = NULL, title = NULL, x_labels_angle = NULL){
  
  if(!is.null(vec_levels)){
    if(length(intersect(names(vec), vec_levels)) != length(union(names(vec), vec_levels))){
      stop("vec_levels should contain all the names of vec")
    } else {
      vec = vec[vec_levels]
    }
  } else {
    vec_levels = gtools::mixedsort(names(vec))
  }
  
  plot = ggplot(mapping = aes(x = factor(names(vec), levels = vec_levels), y = vec)) +
    geom_col() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), axis.title = element_text(size = 20), 
    axis.text = element_text(size = 18)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = xlab, y = ylab, title = title)
  
  if(!is.null(x_labels_angle)){
    plot = plot + theme(axis.text.x = element_text(angle = x_labels_angle, hjust = 1))
  }
  
  return(plot)
}

#' Customize the theme of a ggplot object
#'
#' @param plot A ggplot object
#' @param base_theme The base theme to use for the plot. Default is theme_bw().
#' @param title Title for the plot
#' @param plot_title_size Size of the plot title in pts. Default is 24.
#' @param xlab The x-axis title
#' @param ylab The y-axis title
#' @param axis_title_size Size of the axes titles in pts. Default is 20.
#' @param axis_text_size Size of the axes text in pts. Default is 18.
#' @param x_labels Optional labels to use for the x-axis
#' @param format_x_labels Logical value indicating whether to format the X axis labels by replacing "_" with a space and using title case. Default is FALSE.
#' @param x_labels_angle Number of degrees to rotate x-axis labels. Default is 0. 
#' @param show_legend A logical value indicating whether to show the legend or not. Default is TRUE. 
#' @param fill_title The title for fill in the legend
#' @param color_title The title for color in the legend
#' @param colors An optional vector of colors to use for the color mapping. 
#' @param fill_colors An optional vector of colors to use for the fill. 
#' @param fill_labels An optional vector of labels for the fill legend
#' @param legend_title_size Size of the legend titles in pts. Default is 20.
#' @param legend_text_size Size of the legend text in pts. Default is 18.
#' @param legend_key_size Size of the legend key in cm. Default is 1.
#' @param format_fill_labels Logical value indicating whether to format the legend labels by replacing "_" with a space and using title case. Default is FALSE.
#' @param scale_x A scale to use for the x-axis. If none provided, one is chosen by default.
#' @param scale_y A scale to use for the y-axis. If none provided, one is chosen by default. 
#' @param facet A variable to use for faceting. Default is not to facet. 
#' @param facet_nrow Number of rows to use for faceting. 
#' @param facet_ncol Number of rows to use for faceting. 
#' @param facet_scales scale for facet_wrap(). Default is "free". 
#' @param strip_text_size Size of the facet labels. Default is 20.
#' @param facet_labels Labels to use for panels when facetting
#' @return A ggplot object
#' @export
customize_ggplot_theme = function(plot, base_theme = theme_bw(), title = NULL, plot_title_size = 24,
  xlab = NULL, ylab = NULL, axis_title_size = 20, axis_text_size = 18, x_labels = NULL, format_x_labels = F, x_labels_angle = 0,  
  show_legend = T, fill_title = NULL, color_title = NULL, fill_colors = NULL, colors = NULL, fill_labels = waiver(), 
  legend_title_size = 20, legend_text_size = 18, legend_key_size = 1, format_fill_labels = F, 
  scale_x = NULL, scale_y = NULL, 
  facet = NULL, facet_nrow = NULL, facet_ncol = NULL, facet_scales = "free", strip_text_size = 20, facet_labels = NULL){
  
  plot = plot +
    base_theme +
    theme(plot.title = element_text(hjust = 0.5, size = plot_title_size), 
	    axis.title = element_text(size = axis_title_size), axis.text = element_text(size = axis_text_size), 
  	  legend.title = element_text(size = legend_title_size), legend.text = element_text(size = legend_text_size), 
      legend.key.size = unit(legend_key_size, "cm"), 
      strip.text = element_text(size = strip_text_size)) +
    labs(title = title, x = xlab, y = ylab, fill = fill_title,  color = color_title)
  
  # Remove legend if specified
  if(!show_legend){
    plot = plot + theme(legend.position = "None")
  }
  
  # Rotate x-axis lables if specified
  if(x_labels_angle != 0){
    plot = plot + theme(axis.text.x = element_text(angle = x_labels_angle, hjust = 1))
  }
  
  # # Format x axis labels if specified
  # if(format_fill_labels){
  #   xlabels = ggplot_build(plot)$layout$panel_params[[1]]$x$get_labels()
  #   new_xlabels = stringr::str_to_title(gsub("_", " ", xlabels))
  #   plot = plot + scale_x_discrete(labels = new_xlabels)
  # }
  
  # Use x_labels if provided
  if(!is.null(x_labels)){
    plot = plot + scale_x_discrete(labels = x_labels)
  }
  
  # Format x axis labels if specified
  if(format_x_labels){
    xlabels = ggplot_build(plot)$layout$panel_params[[1]]$x$get_labels()
    new_xlabels = stringr::str_to_title(gsub("_", " ", xlabels))
    plot = plot + scale_x_discrete(labels = new_xlabels)
  }
  
  # Find the type of geom the plot uses
  plot_geom = class(plot$layers[[1]]$geom)
  
  # Change fill colours and labels if specified
  if(!is.null(fill_colors)){
    plot = plot + scale_fill_manual(values = fill_colors, labels = fill_labels)
  } else {
    plot = plot + scale_fill_discrete(labels = fill_labels)
  }
  
  # Change colors if specified
  if(!is.null(colors)){
    plot = plot + scale_color_manual(values = colors)
  }
  
  # If scale_x and scale_y not provided. they are inferred
  if(is.null(scale_x)){
    if("GeomBar" %in% plot_geom & !"GeomCol" %in% plot_geom){
      scale_x = scale_x_continuous(expand = c(0,0), labels = scales::comma) 
    }
  }
  
  if(is.null(scale_y)){
    if("GeomBar" %in% plot_geom){
      scale_y = scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma)
    }
  }
  
  plot = plot + scale_x + scale_y
  
  if(!is.null(facet)){
    if(!is.null(facet_labels)){
      original_facet_labels = sort(unique(plot$data[[facet]]))
      labeller = as_labeller(setNames(facet_labels, original_facet_labels))
    } else {
      labeller = "label_value"
    }
    plot = plot + facet_wrap(facet, nrow = facet_nrow, ncol = facet_ncol, labeller = labeller, scales = facet_scales)
  }
  
  return(plot)
  
}


#' Save a list of plots to a PDF file
#'
#' @param plotlist A list of ggplot objects
#' @param nrows Number of rows in each sheet
#' @param ncols Number of columns in each sheet
#' @param width Width of each page. Default is 16.
#' @param height Heigth of each page. Default is 9.
#' @param filename Name of output file
#' @return Nothing
#' @export
pdf_save = function(plotlist, nrows = 1, ncols = 1, filename, width = 16, height = 9){
  on.exit(expr = dev.off())
  pdf(filename, width = width, height = height, title = basename(filename))
  for(x in seq(1, length(plotlist)/prod(nrows, ncols))) {
    print(ggpubr::ggarrange(plotlist = plotlist[x:(x+prod(nrows, ncols)-1)], nrow = nrows, ncol = ncols))
  }
  dev.off()
  on.exit()
}

#' Save a list of plots to separate files in a specified directory
#'
#' @param plotlist A list of ggplot objects
#' @param filenames A vector of filenames to use to save plots. Default is to use the variable names used for the plots. 
#' @param file_format Extension of file format to use if filenames not provided. Default is "jpeg".
#' @param directory Directory to save images. Default is current directory. 
#' @param width Width of the plots in inches
#' @param height Height of the plots in inches
#' @export
ggsave_list = function(plotlist, filenames = NULL, file_format = "jpeg", directory = ".", width = 16, height = 9){
  
  # Check that filenames has the same length as plotlist
  if(!is.null(filenames) & length(plotlist) != length(filenames)){
    stop("plotlist and filenames must have the same length")
  }
  
  # Check that filenames does not contain duplicates
  if(anyDuplicated(filenames)){stop("filenames cannot contain duplicates")}
  
  # If filenames are not provided, the names of the plot objects are used
  if(is.null(filenames)){
    env_addresses = sapply(ls(envir = .GlobalEnv), function(x) tryCatch(tracemem(get(x)), error = function(x) NA))
    env_addresses = setNames(names(env_addresses), env_addresses)
    plot_addresses = sapply(plotlist, tracemem)
    plot_names = env_addresses[plot_addresses]
    filenames = paste(directory, paste(plot_names, file_format, sep = "."), sep = "/")
    
    sapply(ls(envir = .GlobalEnv), function(x) tryCatch(untracemem(get(x)), error = function(x) NA))
    sapply(plotlist, untracemem)
  }
  
  sapply(seq_along(plotlist), function(x) ggsave(plot = plotlist[[x]], filename = filenames[x], width = width, height = height))
  return(NULL)
  
}

#' Save plots to a PDF file
#'
#' @param plotlist A list of ggplot objects
#' @param titles A vector of titles for the plots
#' @return A ggplot object
#' @export
title_ggplots = function(plotlist, titles){
  
  if(length(plotlist) != length(titles)){stop("plotlist and titles must be the same length")}
    
  lapply(seq_along(plotlist), function(x) plotlist[[x]] + ggtitle(titles[[x]]))
}

#' Encode p-values with symbols
#'
#' @param p_values A vector of p-values
#' @param cutpoints A vector of cutpoints for significance. Default is c(0, 0.001, 0.01, 0.05, 1)
#' @param symbol
#' @return A vector of symbols representing p-values
#' @export
sig_sym = function(p_values, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbol = "*"){
  
  symbols = sapply(3:0, function(x) paste(rep(symbol, x), collapse = ""))
  
  significance_symbols = as.character(symnum(p_values, corr = FALSE, na = FALSE, cutpoints = cutpoints, symbols = symbols))
  return(significance_symbols)
}

#' Display all available shapes and their associated number
#'
#' @export
display_shapes = function(){
  df_shapes = data.frame(shape = 0:24)
  shape_plot = ggplot(df_shapes, aes(0, 0, shape = shape)) +
    geom_point(aes(shape = shape), size = 5, fill = 'red') +
    scale_shape_identity() +
    facet_wrap(~shape) +
    theme_void()
  return(shape_plot)
}

#' Create an upset plot with a title and export
#'
#' @param input_list A named list to use as input for the plot
#' @param title Title for the plot
#' @param title_size Size of the title in points. Default is 30.
#' @param filename Name of file to save plot to. If not provided, file not saved. 
#' @return NULL
#' @export
make_upset_plot = function(input_list, title = NULL, title_size = 30, filename = NULL,
  intersection_size_title_size = 2, intersection_size_tick_labels_size = 2, bar_numbers = 2,
  set_size_title_size = 2, set_size_tick_labels_size = 2, set_names_size = 2){
  
  # add names to list if they are not set
  if(is.null(names(input_list))){
    names(input_list) = paste0("set", 1:length(input_list))
  }
  
  text.scale = c(intersection_size_title_size, intersection_size_tick_labels_size, bar_numbers,
    set_size_title_size, set_size_tick_labels_size, set_names_size)
  
  upset_plot = UpSetR::upset(UpSetR::fromList(input_list), order.by = c("freq", "degree"), 
    decreasing = c(T, T), text.scale = text.scale) 
  
  if(is.null(filename)){
    print(upset_plot)
    grid::grid.text(title, x = 0.65, y = 0.95, gp = grid::gpar(fontsize = title_size))
  } else {
    jpeg(filename = filename, width = 1920, height = 1080, units = "px", bg = "white")
    print(upset_plot)
    grid::grid.text(title, x = 0.65, y = 0.95, gp = grid::gpar(fontsize = title_size))
    dev.off()
  }
  
}

#' Arrange a list of plots
#'
#' @param plotlist A list of ggplot objects
#' @param title Title for arranged plots
#' @param title_size Font size for title. Default is 24
#' @param axis_title_size Size for the axes titles. Default is 18. 
#' @param remove_subtitles A logical value indicating whether or not to remove titles from individual plots. Default is FALSE. 
#' @param new_subtitles An optional vector with new titles to add to the plots. 
#' @param legend Character specifying where to put legend. Default is to the right. Set to "none" to remove.
#' @param subtitle_size Size of the subtitles. Default is 18
#' @param remove_individual_axes_titles A logical value indicating whether to remove individual axes titles. Default is TRUE.
#' @param xaxis_title Title to use for x-axis. 
#' @param yaxis_title Title to use for y-axis
#' @param nrow Number of rows for combined plot. 
#' @param ncol Number of columns for combined plot
#' @param align How to align plots. Same as for ggarrange(). Default is to align vertically and horizontally ("hv").
#' @return A ggarrange object
#' @export
arrange_ggplots = function(plotlist, title = NULL, remove_subtitles = F, new_subtitles = NULL, legend = "right",
  title_size = 24, subtitle_size = 18, remove_individual_axes_titles = T, xaxis_title = NULL, yaxis_title = NULL, axis_title_size = 18, 
  nrow = length(plotlist), ncol = length(plotlist)/nrow, align = "hv"){
  
  # Check that new_subtitles is the same length as plotlist if provided
  if(!is.null(new_subtitles)){
    if(length(new_subtitles) != length(plotlist)){
      stop("new_subtitles must have the same length as plotlist if provided")
    }
  }
  
  # Remove titles from individual plots if specified
  if(remove_subtitles){
    plotlist = lapply(plotlist, function(x) x + labs(title = NULL))
  }
  
  # Remove titles from individual plots if specified
  if(remove_individual_axes_titles){
    plotlist = lapply(plotlist, function(x) x + ggpubr::rremove("xylab"))
  }
  
  # Add new subtitles to plots if provided
  if(!is.null(new_subtitles)){
    plotlist = lapply(seq_along(plotlist), function(x) 
      plotlist[[x]] + ggtitle(new_subtitles[x]))
  }
  
  # Adjust subtitle sizes
  plotlist = lapply(seq_along(plotlist), function(x) 
      plotlist[[x]] + theme(plot.title = element_text(hjust = 0.5, size = subtitle_size)))
  
  # Arrange the plots
  arranged_plot = ggpubr::ggarrange(plotlist = plotlist, nrow = nrow, ncol = ncol, align = align, 
    common.legend = T, legend = legend, label.x = 0.5)
  
  # Add titles to the arranged plot
  arranged_plot = ggpubr::annotate_figure(arranged_plot, 
    top = ggpubr::text_grob(title, size = title_size),
    left = ggpubr::text_grob(yaxis_title, color = "black", rot = 90, size = axis_title_size),
    bottom = ggpubr::text_grob(xaxis_title, color = "black", size = axis_title_size))
  
  return(arranged_plot)
}

colour_list = list(
  ggplot_red_blue = c("#F8766D", "#00BFC4"),
  two_blues = c("#1F8AC0", "#104C91"),
  red_navy_pair = c("#D01C1FFF", "#4B878BFF"),
  gold_green_red = c("#DDCC77", "#117733", "#882255"),
  three_greens = c("#c4e6c3", "#4da284", "#1d4f60"),
  three_blues = c("#d1eeea", "#68abb8", "#2a5674"),
  three_reds = c("#FCBBA1", "#FB6A4A", "#A50F15"),
  three_purples = c("#BCBDDC", "#807DBA", "#54278F"),
  rainbox_six = c("#9E0142", "#FDAE61", "#FEE08B", "#ABDDA4", "#3288BD", "#5E4FA2"),
  nine_greens = RColorBrewer::brewer.pal(9, "YlGn"),
  purple_and_gold_dark = c("#7B5C90", "#bfab25"),
  purple_and_gold_light = c("#A28CB1", "#D2C465")
)

#' A list of colour vectors
#'
#' @format A named list of colour vectors
"colour_list"


# show_palettes = function(){
#   
#   par(mfrow = c(ceiling(sqrt(length(colour_palettes))), floor(sqrt(length(colour_palettes)))))
#   lapply(names(colour_palettes), function(x) {scales::show_col(colour_palettes[[x]]); title(main = x, outer = T)})
#   par(mfrow=c(1,1)) 
# }
