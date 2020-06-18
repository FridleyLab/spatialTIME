#' Generate plot of TMA point process
#'
#' @describtion This function generates plot of point process in rectangular or circular window.
#' @param dlist List of TMA/ROI data frames or single TMA/ROI data frame.
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
#' @examples
#' plot_immunoflo(tma.data[[1]],
#'  .id = "subid.x", mnames = marker_names, wshape = "r")
#'
plot_immunoflo <- function(
  dlist,
  plot_title, 
  mnames, 
  # filename, 
  mlabels = NULL, 
  # pretty_labels = TRUE,
  mcolors = NULL, 
  dark_mode = FALSE,
  cell_type = NULL#, 
  # path = NULL
  ) {
  # convert to list of dataframes - throw error message if missing
  if (missing(dlist)) stop("dlist is missing; must provide data for plotting")
  if (is.data.frame(dlist)) dlist = list(dlist)
  
  # if (missing(filename)) stop("filename is missing; filename must be a string")
  
  # check if marker names and labels are equal in length
  if (!is.null(mlabels) && length(mnames) != length(mlabels)) {
    stop("mnames and mlabels must be equal in length")
  }
  
  # set marker label to marker name if no marker label provided
  if (is.null(mlabels)) {
    mlabels = mnames
    # pretty up labels 
    mlabels = stringr::str_wrap(mlabels,20)
  }
  
  # progress bar for creating plots
  pb <- dplyr::progress_estimated(length(dlist))
  
  # plots
  plot <- lapply(dlist, function(x){
    # update progress bar
    pb$tick()$print()
    # data to generate plot
    plot_data <- x %>% 
      dplyr::select(
        !!plot_title, XMin, XMax, YMin, YMax, !!mnames, !!cell_type) %>% 
      tidyr::pivot_longer(cols = !!mnames,
                          names_to = "marker", values_to = "indicator") %>% 
      dplyr::mutate(xloc = (XMin + XMax) / 2,
                    yloc = (YMin + YMax) / 2,
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
      dplyr::filter(indicator == 1) %>% 
      ggplot2::ggplot(aes(x = xloc, y = yloc, color = marker,
                          shape = !!as.name(cell_type))) +
      ggplot2::geom_point(data = filter(plot_data, indicator == 0),
                          # aes(fill = "grey70"),
                          color = "gray70") +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(5)) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(5)) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::scale_color_manual(NULL, values = mcolors, drop = FALSE) +
      ggplot2::scale_shape_manual(NULL, values = c(3, 16), drop = FALSE) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank(),
            panel.grid = element_blank())
    
    if(dark_mode == TRUE){
      basic_plot <- basic_plot + theme_dark_mode()
    }
    
    return(basic_plot)
    
  })
  
  # # file name with full path
  # if (!is.null(path)) {
  #   filename <- file.path(path, sprintf("%s.pdf", filename))
  # } else {
  #   filename <- sprintf("%s.pdf", filename)
  # }
  # 
  # # progress bar for saving plots
  # pb <- dplyr::progress_estimated(length(plot))
  # 
  # # save plots with probress bar
  # pdf(filename, height = 10, width = 10)
  # on.exit(dev.off())
  # invisible(
  #   lapply(seq_along(plot), function(x) {
  #     # update progress bar
  #     pb$tick()$print()
  #     
  #     print(plot[[x]])
  #   })
  # )
}
