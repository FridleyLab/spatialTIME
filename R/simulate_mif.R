#' Simulate Multiplex Immunoflourescent data
#'
#' @description Simulate IF data 
#' @param k numeric value ...
#' @param xmin numeric value ...
#' @param xmax numeric value ...
#' @param ymin numeric value ...
#' @param ymax numeric value ...
#' @param sdmin numeric value ...
#' @param sdmax numeric value ...
#' @param k_cell_stroma numeric value ...
#' @param k_cell_tumor numeric value ...
#' @param upper_lim_tumor numeric value ...

simulate_mif <- function(k = 10,
                         xmin = 0,
                         xmax = 10,
                         ymin = 0,
                         ymax = 10,
                         sdmin = 0.5, 
                         sdmax = 2,
                         k_cell_stroma = 3,
                         k_cell_tumor = 5,
                         upper_lim_tumor = 0.1){
  
  center_peak <- modes_stroma(k = k, xmin = xmin, xmax = xmax, ymin = ymin, 
                              ymax = ymax, sdmin = sdmin, sdmax = sdmax)
  center_hole <- hole(xmin = xmin, xmax = xmax, ymin = ymin, 
                      ymax = ymax, sdmin = sdmin, sdmax = sdmax)
  center_pos_cell_tumor <- modes_stroma(k = k_cell_tumor, xmin = xmin, 
                                        xmax = xmax, ymin = ymin, 
                                        ymax = ymax, sdmin = 10, sdmax = 12)
  center_pos_cell_stroma <- modes_stroma(k = k_cell_stroma, xmin = xmin, 
                                         xmax = xmax, ymin = ymin, 
                                         ymax = ymax, sdmin = 2, sdmax = 3)
  
  pp <- spatstat::rpoispp(25, win = spatstat::owin(c(0,10),c(0,10))) %>% 
    data.frame() 
  
  # for(i in 1:nrow(pp)){
  #   
  #   dist <- sapply(1:nrow(center_peak), dist_function, 
  #                  other_data = center_peak)
  #   
  #   pp$prob_stroma[i] = max(dist)
  #   
  #   pp$prob_stroma_label[i] = sample(c('Stroma', 'Tumor'), 1,
  #                                    prob = c(pp$prob_stroma[i],
  #                                             (1-pp$prob_stroma[i])))
  #   
  #   dist <- sapply(1:nrow(center_pos_cell_tumor), dist_function,
  #                  other_data = center_pos_cell_tumor)
  #   
  #   pp$prob_positive_tumor[i] = max(max(dist) * (upper_lim_tumor - 
  #                                                  pp$prob_stroma[i]), 0)
  #   
  #   dist <- sapply(1:nrow(center_pos_cell_stroma), dist_function,
  #                  other_data = center_pos_cell_stroma)
  #   
  #   pp$prob_positive_stroma[i] = max(dist) * pp$prob_stroma[i]
  #   
  #   # If you want to have a constant proability for all cells to be positive, 
  #   # then set pp$prob_positive_stroma[i] to a constant value
  #   pp$prob_positive_label[i] = ifelse(pp$prob_stroma_label[i] == 'Stroma',
  #                                      sample(c('Positive', 'Negative'), 1,
  #                                             prob = c(pp$prob_positive_stroma[i], (1-pp$prob_positive_stroma[i]))),
  #                                      sample(c('Positive', 'Negative'), 1,
  #                                             prob = c(pp$prob_positive_tumor[i], (1-pp$prob_positive_tumor[i]))))
  #   
  # }
  
  return(pp)
  
}


# #########################################################################################################################
# #########################################################################################################################
# #Generates the plots for entire region
# 
# grid = expand.grid(x = seq(0,10,0.1), y = seq(0,10,0.1))
# for(i in 1:nrow(grid)){
#   dist = sapply(1:nrow(center_peak), function(a){
#     diff = c((grid[i,1] - center_peak[a,1]), (grid[i,2] - center_peak[a,2]))
#     sigma = matrix(c(center_peak$sd.x[a]^2, rep(center_peak$rho[a], 2),
#                      center_peak$sd.y[a]^2), nrow = 2, ncol = 2)
#     z = t(diff) %*% solve(sigma) %*% diff
#     closest = exp(-z)
#   })
#   grid$prob_stroma[i] = max(dist)
#   
#   dist = sapply(1:nrow(center_pos_cell_tumor), function(a){
#     diff = c((grid[i,1] - center_pos_cell_tumor[a,1]), (grid[i,2] - center_pos_cell_tumor[a,2]))
#     sigma = matrix(c(center_pos_cell_tumor$sd.x[a]^2, rep(center_pos_cell_tumor$rho[a], 2),
#                      center_pos_cell_tumor$sd.y[a]^2), nrow = 2, ncol = 2)
#     z = t(diff) %*% solve(sigma) %*% diff
#     closest = exp(-z)
#   })
#   grid$prob_positive_tumor[i] = max(max(dist) * (upper_lim_tumor - grid$prob_stroma[i]),0)
#   
#   dist = sapply(1:nrow(center_pos_cell_stroma), function(a){
#     diff = c((grid[i,1] - center_pos_cell_stroma[a,1]), (grid[i,2] - center_pos_cell_stroma[a,2]))
#     sigma = matrix(c(center_pos_cell_stroma$sd.x[a]^2, rep(center_pos_cell_stroma$rho[a], 2),
#                      center_pos_cell_stroma$sd.y[a]^2), nrow = 2, ncol = 2)
#     z = t(diff) %*% solve(sigma) %*% diff
#     closest = exp(-z)
#   })
#   grid$prob_positive_stroma[i] = max(dist) * grid$prob_stroma[i]
#   
# }
# 
# contours_stroma = ggplot(data = grid, aes(x = x, y = y, z = 100*prob_stroma)) + 
#   geom_contour_filled(breaks = seq(0,100,10)) + 
#   theme_bw() + labs(fill = 'Probability of \nStroma Cell') + 
#   theme(panel.grid = element_blank(), legend.position = 'right')  + 
#   theme(legend.title = element_text(size = 16),legend.text = element_text(size = 16))
# 
# contours_pos_stoma = ggplot(data = grid, aes(x = x, y = y, z = 100*prob_positive_stroma)) + 
#   geom_contour_filled(breaks = seq(0,100,10)) + 
#   theme_bw() + labs(fill = 'Probability of \nPositive Cell\nin Stroma ') + 
#   theme(panel.grid = element_blank(), legend.position = 'right')  + 
#   theme(legend.title = element_text(size = 16),legend.text = element_text(size = 16))
# 
# contours_pos_tumor = ggplot(data = grid, aes(x = x, y = y, z = 100*prob_positive_tumor)) + 
#   geom_contour_filled(breaks = seq(0,100,10)) + 
#   theme_bw() + labs(fill = 'Probability of \nPositive cell\nin Tumor') + 
#   theme(panel.grid = element_blank(), legend.position = 'right')  + 
#   theme(legend.title = element_text(size = 16),legend.text = element_text(size = 16)) 
# 
# pp_plot = ggplot(data = pp, aes(x = x, y = y, color = interaction(prob_stroma_label,prob_positive_label))) + 
#   geom_point()
# 
# ggarrange(plotlist = list(contours_stroma, contours_pos_stoma, contours_pos_tumor, pp_plot), nrow = 2, ncol = 2)
# table(pp$prob_stroma_label, pp$prob_positive_label)
# 
