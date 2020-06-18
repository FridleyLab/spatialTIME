# Dark theme for plotting
theme_dark_mode <- function () { 
  ggdark::dark_mode() %+replace% 
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}