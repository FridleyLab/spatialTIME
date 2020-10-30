#' Generate plot of TMA point process
#'
#' @description This function generates plot of point process in rectangular or circular window.
#' @param mif MIF object created using create_MIF().
#' @param plot_title Character string or vector of character strings of variable name(s) to serve as plot title(s).
#' @param mnames Character vector containing marker names.
#' @param filename Character string of file name to store plots. Plots are generated as single .pdf file.
#' @param mlabels Character vector of label for marker names to display in plot legend. 
#' @param mcolors Character vector of color names to display markers in the plot.
#' @param cell_type Character vector of cell type ???
#' @param path Different path than file name or to use in conjunction with filename ???
#' @param dark_mode Plot using dark color scheme 
#' 
#' @return A list of ggplots (the length of dlist) - one for each TMA or ROI
#'    
#' @export
#'
plot_immunoflo <- function(
  mif,
  plot_title, 
  mnames, 
  mlabels = NULL, 
  # pretty_labels = TRUE,
  mcolors = NULL, 
  dark_mode = FALSE,
  cell_type = NULL, 
  filename = NULL,
  path = NULL
  ) {
  
  ### changes to make
  # 1. input will be new data type 
  # 2. plot to pdf unless given subset of samples 
  
  
  
  # convert to list of dataframes - throw error message if missing
  if (missing(mif)) stop("MIF is missing; please provide the appropriate data")
  if (class(mif) != "mif") stop("Please use a mif object")
  # if (is.data.frame()) dlist = list(dlist) - need to change to MIF object
  
  # if (missing(filename)) stop("filename is missing; filename must be a string")
  
  # check if marker names and labels are equal in length
  if (!is.null(mlabels) && length(mnames) != length(mlabels)) {
    stop("gene names and labels must be equal in length")
  }
  
  # set marker label to marker name if no marker label provided
  if (is.null(mlabels)) {
    mlabels = mnames
    # pretty up labels 
    mlabels = stringr::str_wrap(mlabels,20)
  }
  
  # progress bar for creating plots
  pb <- dplyr::progress_estimated(length(mif[["spatial"]]))
  
  # plots
  plot <- lapply(mif[["spatial"]], function(x){
    # update progress bar
    pb$tick()$print()
    # data to generate plot
    plot_data <- x %>% 
      janitor::clean_names() %>%
      dplyr::select(
        !!plot_title, x_min, x_max, y_min, y_max, !!mnames, !!cell_type) %>% 
      tidyr::pivot_longer(cols = !!mnames,
                          names_to = "marker", values_to = "indicator") %>% 
      dplyr::mutate(xloc = (x_min + x_max) / 2,
                    yloc = (y_min + y_max) / 2,
                    marker = factor(
                      marker, levels = mnames, labels = mlabels)) 
    
    # plot title
    plot_title <- if (length(plot_title) == 1) {
      paste0("ID: ", unique(x[[plot_title]]))
    } else {
      paste0("ID: ", paste(unique(x[, plot_title]), collapse = ", "))
    }
    
    # color palette
    if (is.null(mcolors)) {
      # viridis prints better (BW) and easily read (colorblindness)
      mcolors = viridisLite::viridis(length(mnames), option = "viridis")
      # mcolors = RColorBrewer::brewer.pal(length(mnames), "Paired")
    }
    
    basic_plot <- plot_data %>% 
      dplyr::filter(.data$indicator == 1) %>% 
      ggplot2::ggplot(ggplot2::aes(x = xloc, 
                                   y = yloc, 
                                   color = marker, 
                                   shape = !!as.name(cell_type))) +
      # ggplot2::geom_point(data = filter(plot_data, indicator == 0),
      #                     # aes(fill = "grey70"),
      #                     color = "gray70") +
      ggplot2::geom_point(data = plot_data[plot_data$indicator == 0,],
                          # aes(fill = "grey70"),
                          color = "gray70") +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(5)) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(5)) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::scale_color_manual(NULL, values = mcolors, drop = FALSE) +
      ggplot2::scale_shape_manual(NULL, values = c(3, 16), drop = FALSE) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(axis.title = ggplot2::element_blank(),
                    panel.grid = ggplot2::element_blank())
    
    if(dark_mode == TRUE){
      basic_plot <- basic_plot + ggdark::theme_dark_mode()
    }
    
    return(basic_plot)
    
  })
  
  # output to pdf if filename is specified 
  if(!is.null(filename)){
    grDevices::pdf(sprintf("%s.pdf",filename), height = 10, width = 10)
    on.exit(dev.off())
    invisible(
      lapply(seq_along(plot), function(x) {
        print(plot[[x]])
      })
    )
    dev.off()
  }
  
  return(plot)
  
}
