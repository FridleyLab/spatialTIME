## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

set.seed(126)
library(spatialTIME)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyselect)
library(RColorBrewer)
library(pheatmap)

getColors = function(n){
  #https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  sample(col_vector, n)
}

## ----create-------------------------------------------------------------------

# Make sure the variable types are the same for deidentified_id and 
# deidentified_sample in their corresponding datasets
x <- create_mif(clinical_data = example_clinical %>% 
                  mutate(deidentified_id = as.character(deidentified_id)),
                sample_data = example_summary %>% 
                  mutate(deidentified_id = as.character(deidentified_id)),
                spatial_list = example_spatial,
                patient_id = "deidentified_id", 
                sample_id = "deidentified_sample")

x #prints a summary of how many patients, samples, and spatial files are present


## ----plot_bad, fig.width=14, fig.height=6,fig.align = 'left', fig.cap='Notice that using the "bad" (B) order that one would have the impression that there is a huge number of CD3+ and there are no cells that are CD3+ FOXP3+ or CD3+ CD8+. While the "good" (G) order can  better make this distiction as well as show cells that are only positive for CD3 and not postive for FOXP3 or CD8.'----
mnames_bad <- c("CD3..CD8.","CD3..FOXP3.","CD3..Opal.570..Positive",
                "CD8..Opal.520..Positive","FOXP3..Opal.620..Positive", 
                "PDL1..Opal.540..Positive", "PD1..Opal.650..Positive")

# Used to make the legends in both plots below be in same order and use the 
# same coloring scheme for the purpose making a common legend

values =  RColorBrewer::brewer.pal(length(mnames_bad), "Accent")
names(values) = mnames_bad

#add an element in the `derived` object position
x<- plot_immunoflo(x, plot_title = "deidentified_sample",  mnames = mnames_bad,
                   cell_type = "Classifier.Label")

bad_names <- x[["derived"]][["spatial_plots"]][[4]] + 
  theme(legend.position = 'bottom') + 
  scale_color_manual(breaks = mnames_bad,
                     values = values,
                     labels = mnames_bad %>%
                       gsub("..Opal.*", "+", .) %>% 
                       gsub("\\.\\.", "+", .) %>% 
                       gsub("\\.", "+", .)) 

mnames_good <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
                 "FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
                 "PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")

x <- plot_immunoflo(x, plot_title = "deidentified_sample", mnames = mnames_good, 
                    cell_type = "Classifier.Label")

good_names <- x[["derived"]][["spatial_plots"]][[4]] + 
  theme(legend.position = 'bottom') + 
  scale_color_manual(breaks = mnames_good, 
                     values = values[match(mnames_good, names(values))],
                     labels = mnames_good %>%
                       gsub("..Opal.*", "+", .) %>% 
                       gsub("\\.\\.", "+", .) %>% 
                       gsub("\\.", "+", .))

x$sample %>% filter(deidentified_sample == 'TMA3_[9,K].tif') %>% select(c(2, 4:15)) %>%
  pivot_longer(cols = 2:13, names_to = 'Marker', values_to = 'Count')

## -----------------------------------------------------------------------------
gridExtra::grid.arrange(bad_names, good_names, ncol=2)

## ----ripleys, warning = FALSE, fig.width=10, fig.height=6, fig.align = 'center'----
x <- ripleys_k(mif = x, mnames = mnames_good, 
               num_permutations = 10, method = "K",
               edge_correction = 'translation', r_range = 0:100,
               permute = TRUE,
               keep_permutation_distribution = FALSE, overwrite = TRUE, workers = 1)

# This will keeps the colors in evx$ery plot for the remainder of the vignette compatible 
values = RColorBrewer::brewer.pal(length(unique(x$derived$univariate_Count[[x$sample_id]])), "Accent")
names(values) = unique(x$derived$univariate_Count$deidentified_sample)

x$derived$univariate_Count  %>%
  filter(Marker != 'PDL1..Opal.540..Positive') %>%
  ggplot(aes(x = r, y = `Degree of Clustering Permutation`)) +
  geom_line(aes(color = get(x$sample_id)), show.legend = FALSE) +
  facet_wrap(Marker~., scales = 'free') + theme_bw() + 
  scale_color_manual(values = values)
  


## ----ripleys_exact, warning = FALSE, fig.width=10, fig.height=6, fig.align = 'center'----
x <- ripleys_k(mif = x, mnames = mnames_good, 
               num_permutations = 10, method = "K",
               edge_correction = 'translation', r_range = 0:100,
               permute = FALSE,
               keep_permutation_distribution = FALSE, overwrite = TRUE, workers = 1)

# This will keeps the colors in evx$ery plot for the remainder of the vignette compatible 
values = RColorBrewer::brewer.pal(length(unique(x$derived$univariate_Count[[x$sample_id]])), "Accent")
names(values) = unique(x$derived$univariate_Count$deidentified_sample)

x$derived$univariate_Count  %>%
  filter(Marker != 'PDL1..Opal.540..Positive') %>%
  ggplot(aes(x = r, y = `Degree of Clustering Exact`)) +
  geom_line(aes(color = get(x$sample_id)), show.legend = FALSE) +
  facet_wrap(Marker~., scales = 'free') + theme_bw() + 
  scale_color_manual(values = values)
  

## ----fig.width=10, fig.height=6, fig.align = 'center'-------------------------
x$derived$univariate_Count %>%
  filter(Marker == 'CD3..CD8.') %>%
  inner_join(x$clinical,.) %>%
  ggplot(aes(shape = status, y = `Degree of Clustering Exact`, x =r)) +
  geom_point(aes(color = get(x$sample_id))) +
  theme_bw() + scale_color_manual(values = values)

## ----fig.width=10, fig.height=6, fig.align = 'center'-------------------------
x <- bi_ripleys_k(mif = x, mnames = mnames_good, 
                  num_permutations = 10, 
                  permute=TRUE,
                  edge_correction = 'translation', r_range = 0:100,
                  keep_permutation_distribution = FALSE, workers = 1)

x$derived$bivariate_Count  %>%
  filter(Anchor == 'CD3..FOXP3.',
         Counted == 'CD3..CD8.') %>%
  ggplot(aes(x = r, y = `Degree of Clustering Permutation`)) +
  geom_line(aes(color = get(x$sample_id)), show.legend = TRUE) +
  theme_bw() + scale_color_manual(values = values)

## ----fig.width=10, fig.height=6, fig.align = 'center'-------------------------
x <- bi_ripleys_k(mif = x, mnames = mnames_good, 
                  num_permutations = 10, 
                  permute=FALSE,
                  edge_correction = 'translation', r_range = 0:100,
                  keep_permutation_distribution = FALSE, 
                  overwrite = TRUE, workers = 1)

x$derived$bivariate_Count  %>%
  filter(Anchor == 'CD3..FOXP3.',
         Counted == 'CD3..CD8.') %>%
  ggplot(aes(x = r, y = `Degree of Clustering Exact`)) +
  geom_line(aes(color = get(x$sample_id)), show.legend = TRUE) +
  theme_bw() + scale_color_manual(values = values)

## ----NN, warning = FALSE, fig.width=10, fig.height=6, fig.align = 'center'----

x <- NN_G(mif = x, mnames = mnames_good, num_permutations = 10,
                edge_correction = 'rs', r = 0:100, workers = 1)

x$derived$univariate_NN  %>%
  filter(Marker != 'PDL1..Opal.540..Positive') %>%
  ggplot(aes(x = r, y = `Degree of Clustering Permutation`)) +
  geom_line(aes(color = get(x$sample_id))) +
  facet_wrap(Marker~., scales = 'free') + theme_bw() + 
  scale_color_manual(values = values)
  


## ----fig.width=10, fig.height=6, fig.align = 'center'-------------------------
x <- bi_NN_G(mif = x, mnames = c("CD3..CD8.", "CD3..FOXP3."), num_permutations = 10,
               edge_correction = 'rs', r = 0:100, workers = 1, overwrite = TRUE)

x$derived$bivariate_NN  %>%
  filter(Anchor == 'CD3..FOXP3.') %>%
  ggplot(aes(x = r, y = `Degree of Clustering Permutation`)) +
  geom_line(aes(color = deidentified_sample), show.legend = TRUE) +
  theme_bw() +  scale_color_manual(values = values)

## ----fig.width=10, fig.height=6, fig.align = 'center', fig.cap='Cell locations for a particular TMA core where the blue cells are CD3+ cells and the grey cells are CD3- cells. Notice that there is a much larger quantity on the left than the right half indicating that there is some clustering occuring.'----
x <- plot_immunoflo(x, plot_title = "deidentified_sample", mnames = mnames_good[1], 
                    cell_type = "Classifier.Label")

x[["derived"]][["spatial_plots"]][[4]]

## ----fig.width=10, fig.height=6, fig.align = 'center', fig.cap = 'Histogram displaying how much the permuted estimates of K (100 permutations) can vary. This also illustrates how different the permutation distribution and its mean (green vertical line) can be from the theoretical value (red vertical line) or observed value (blue vertical line). Lastly, this plot confirms what we see visually above in that the CD3+ cells are not uniformly distributed throughout this slide, but clustered.'----
x <- ripleys_k(mif = x, mnames = mnames_good[1], num_permutations = 100,
               edge_correction = 'translation', r = c(0,50),
               permute = TRUE,
               keep_permutation_distribution = TRUE, workers = 1, overwrite = TRUE)

#Computed the mean of the permutation distribution for a particular image
perm_mean  = x$derived$univariate_Count %>%
  filter(deidentified_sample == 'TMA3_[9,K].tif',
         r == 50)  %>%
  summarize(`Permuted CSR` = mean(`Permuted CSR`)) %>%
  unlist()

x$derived$univariate_Count %>%
  filter(deidentified_sample == 'TMA3_[9,K].tif',
         r == 50) %>%
  ggplot(aes(x = `Permuted CSR`)) + 
  geom_histogram(color = 'white', bins = 50) + 
  geom_vline(aes(xintercept = `Theoretical CSR`), color = 'red', size = 1.5) + 
  geom_vline(aes(xintercept = `Observed K`), color = 'blue', size = 1.5) + 
  geom_vline(xintercept = perm_mean, color = 'green', size = 1.5) + 
  theme_bw()



## -----------------------------------------------------------------------------
x_tumor = subset_mif(x, classifier = 'Classifier.Label', level = 'Tumor', markers = mnames_good)

table(x$spatial[[1]]$Classifier.Label)
table(x_tumor$spatial[[1]]$Classifier.Label)

## ----fig.width=10, fig.height=6, fig.align = 'center', fig.cap='Scatter plots showing the difference in degree of spatial clustering for CD3+, and CD8+ cells. These markers were selected since all 5 examples spatial files had at least 2 to these cells.'----
x <- ripleys_k(mif = x, mnames = mnames_good, num_permutations = 10,
               edge_correction = 'translation', r = c(0,50), permute = T,
               keep_permutation_distribution = FALSE, workers = 1, overwrite = TRUE)

inner_join(x$clinical, x$derived$univariate_Count) %>%
  filter(grepl('CD3..O|CD8..O', Marker), r == 50) %>% 
  ggplot(aes(x = status, y = `Degree of Clustering Permutation`)) + 
  geom_point(aes(color = deidentified_sample), show.legend = FALSE) + 
  facet_wrap(Marker ~., scales = 'free') 


## ----fig.width=10, fig.height=6, fig.align = 'center', fig.cap='Heatmaps showing a potential descriptive visualization for identifying difference patterns within multiple markers across different populations.'----

Rip_K_df = x$derived$univariate_Count %>% 
  filter(!(Marker %in% c('PDL1..Opal.540..Positive', 
                         'CD3..FOXP3.',
                         'PD1..Opal.650..Positive')), r == 50) %>%
  select(deidentified_sample, Marker, `Degree of Clustering Permutation`) %>%
  pivot_wider(names_from = deidentified_sample, values_from = `Degree of Clustering Permutation`,
              id_cols = Marker)

Rip_K_matrix = Rip_K_df %>% 
  ungroup() %>%
  select(-Marker) %>%
  as.matrix()

rownames(Rip_K_matrix) = Rip_K_df$Marker

annotation = inner_join(x$clinical, x$derived$univariate_Count) %>%
  select(deidentified_sample, status) %>%
  filter(!duplicated(.)) %>%
  data.frame()
  

rownames(annotation) = annotation$deidentified_sample

annotation = annotation %>%
  select(status)
  
pheatmap(Rip_K_matrix, treeheight_row = 0, treeheight_col = 0, 
                   cluster_cols = FALSE, annotation = annotation)

## -----------------------------------------------------------------------------
x <- ripleys_k(mif = x, mnames = mnames_good[1], num_permutations = 100,
               edge_correction = 'translation', r = c(0,50),permute = TRUE, 
               keep_permutation_distribution = TRUE, workers = 1, overwrite = TRUE,
               xloc = 'XMin', yloc = 'YMin')

