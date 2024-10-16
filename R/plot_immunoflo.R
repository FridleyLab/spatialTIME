#' Generate plot of TMA point process
#'
#' @description This function generates plot of point process in rectangular or circular window.
#' @param mif MIF object created using create_MIF().
#' @param plot_title Character string or vector of character strings of variable name(s) to serve as plot title(s).
#' @param mnames Character vector containing marker names.
#' @param filename Character string of file name to store plots. Plots are generated as single .pdf file.
#' @param mcolors Character vector of color names to display markers in the plot.
#' @param cell_type Character vector of cell type
#' @param path Different path than file name or to use in conjunction with filename ???
#' @param xloc,yloc columns in the spatial files containing the x and y locations of cells. Default is `NULL` which will result in `xloc` and `yloc` being calculated from `XMin`/`YMin` and `XMax`/`YMax`
#' 
#' @return mif object and the ggplot objects can be viewed form the derived slot of the mif object
#' 
#' @importFrom grDevices dev.off
#'    
#' @export
#'
#' @examples
#' #Create mif object
#' library(dplyr)
#' x <- create_mif(clinical_data = example_clinical %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' sample_data = example_summary %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' spatial_list = example_spatial,
#' patient_id = "deidentified_id", 
#' sample_id = "deidentified_sample")
#' 
#' mnames_good <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
#' "FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
#' "PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")
#' 
#' x <- plot_immunoflo(x, plot_title = "deidentified_sample", mnames = mnames_good, 
#' cell_type = "Classifier.Label")
#' 
#' x[["derived"]][["spatial_plots"]][[4]]


plot_immunoflo <- function(
  mif,
  plot_title, 
  mnames, 
  # pretty_labels = TRUE,
  mcolors = NULL, 
  cell_type = NULL, 
  filename = NULL,
  path = NULL,
  xloc = NULL,
  yloc = NULL) {
  
  ### changes to make
  # 1. input will be new data type 
  # 2. plot to pdf unless given subset of samples 
  
  # convert to list of dataframes - throw error message if missing
  if (missing(mif)) stop("MIF is missing; please provide the appropriate data")
  if (!is(mif, "mif")) stop("Please use a mif object")
  # if (is.data.frame()) dlist = list(dlist) - need to change to MIF object
  
  # if (missing(filename)) stop("filename is missing; filename must be a string")
  
  # plots
  plot <- pbmcapply::pbmclapply(mif[["spatial"]], function(x){
    #make the xloc and yloc columns
    if(is.null(xloc) | is.null(yloc)){
      x = x %>%
        dplyr::mutate(xloc = (XMax + XMin)/2,
                      yloc = (YMax + YMin)/2)
    } else {
      #rename columns to follow xloc and yloc names
      x = x %>%
        dplyr::rename("xloc" = !!xloc, 
                      "yloc" = !!yloc)
    }
    
    # data to generate plot
    plot_data <- x %>%
      dplyr::select(dplyr::any_of(c(!!plot_title, !!mnames, !!cell_type)),
                    xloc, yloc) %>% 
      tidyr::pivot_longer(cols = !!mnames,
                          names_to = "marker", values_to = "indicator") %>% 
      dplyr::mutate(marker = factor(marker, levels = mnames))
    
    # plot title
    plot_title <- if (length(plot_title) == 1) {
      paste0("ID: ", unique(x[[plot_title]]))
    } else {
      paste0("ID: ", paste(unique(x[, plot_title]), collapse = ", "))
    }
    
    # color palette
    if (is.null(mcolors)) {
      # viridis prints better (BW) and easily read (colorblindness)
      #mcolors = viridisLite::viridis(length(mnames), option = "viridis")
      
      #Set2 is not particularly pretty, but it is colorblind friendly
      #Ram had used 'Paired' which is also colorblind friendly.
      
      mcolors = RColorBrewer::brewer.pal(length(mnames), "Paired")
    }
    
    if(is.null(cell_type)){
      basic_plot <- plot_data %>% 
        dplyr::filter(indicator == 1) %>% 
        ggplot2::ggplot(ggplot2::aes(x = xloc, 
                                     y = yloc, 
                                     color = marker)) +
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
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(axis.title = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank()) 
    }else{
    basic_plot <- plot_data %>% 
      dplyr::filter(indicator == 1) %>% 
      ggplot2::ggplot(ggplot2::aes(x = xloc, 
                                   y = yloc, 
                                   color = marker, 
                                   shape = cell_type)) +
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
    }
    basic_plot = basic_plot + 
      ggplot2::scale_y_reverse()
    return(basic_plot)
    
  }, mc.cores = 1)
  
  # output to pdf if filename is specified 
  if(!is.null(filename)){
    grDevices::pdf(sprintf("%s.pdf",filename), height = 10, width = 10)
    on.exit(dev.off())
    invisible(
      lapply(seq_along(plot), function(x) {
        print(plot[[x]])
      })
    )
    grDevices::dev.off()
  }
  
  mif$derived$spatial_plots = plot
  
  return(mif)
  
}
