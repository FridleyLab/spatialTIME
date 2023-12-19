## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tidyverse)
devtools::load_all()
mif = create_mif(clinical_data = example_clinical %>% 
                   mutate(deidentified_id = as.character(deidentified_id)),
                 sample_data = example_summary %>% 
                   mutate(deidentified_id = as.character(deidentified_id)),
                 spatial_list = example_spatial,
                 patient_id = "deidentified_id", 
                 sample_id = "deidentified_sample")

markers = colnames(mif$spatial[[1]]) %>% 
  grep("CD3|Pos", ., value = T) %>% 
  grep("Cyto|Nucle", ., value = T, invert = T)
markers = markers[c(1,2,4,5,8)]

## -----------------------------------------------------------------------------
mif = ripleys_k(mif = mif,
                mnames = markers[1:2], 
                r_range = 0:100, 
                num_permutations = 25,
                edge_correction = 'translation', 
                permute = TRUE, 
                keep_permutation_distribution = FALSE, 
                workers = 1, 
                overwrite = TRUE, 
                xloc = NULL, 
                yloc = NULL)
mif$derived$univariate_Count %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Permutation`, color = deidentified_sample)) +
  facet_grid(~Marker)

## -----------------------------------------------------------------------------
mif = bi_ripleys_k(mif = mif, 
                   mnames = markers[1:2],
                   r_range = 0:100,
                   num_permutations = 25,
                   edge_correction = "translation",
                   permute = TRUE,
                   keep_permutation_distribution = FALSE,
                   workers = 1, 
                   overwrite = TRUE,
                   xloc = NULL, 
                   yloc = NULL)

mif$derived$bivariate_Count %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Permutation`, color = deidentified_sample)) +
  facet_grid(~Anchor)

## -----------------------------------------------------------------------------
mif = NN_G(mif = mif, 
           mnames = markers[1:2], 
           r_range = 0:100, 
           num_permutations = 25, 
           edge_correction = "rs", 
           keep_perm_dis = FALSE, 
           workers = 1, 
           overwrite = TRUE, 
           xloc = NULL,
           yloc = NULL)

mif$derived$univariate_NN %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Permutation`, color = deidentified_sample)) +
  facet_grid(~Marker)

## -----------------------------------------------------------------------------
mif = bi_NN_G(mif = mif, 
           mnames = markers[1:2], 
           r_range = 0:100, 
           num_permutations = 25, 
           edge_correction = "rs", 
           keep_perm_dis = FALSE, 
           workers = 1, 
           overwrite = TRUE, 
           xloc = NULL,
           yloc = NULL)

mif$derived$bivariate_NN %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Permutation`, color = deidentified_sample)) +
  facet_grid(~Anchor)

## -----------------------------------------------------------------------------
mif = pair_correlation(mif = mif, 
                       mnames = markers[1:2],
                       r_range = 0:100, 
                       num_permutations = 25, 
                       edge_correction = "translation", 
                       keep_permutation_distribution = FALSE, 
                       workers = 1, 
                       overwrite = TRUE, 
                       xloc = NULL, 
                       yloc = NULL)

mif$derived$univariate_pair_correlation %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Correlation Permuted`, color = deidentified_sample)) +
  facet_grid(~Marker)

## -----------------------------------------------------------------------------
mif = bi_pair_correlation(mif = mif, 
                       mnames = markers[1:2],
                       r_range = 0:100, 
                       num_permutations = 25, 
                       edge_correction = "translation", 
                       keep_permutation_distribution = FALSE, 
                       workers = 1, 
                       overwrite = TRUE, 
                       xloc = NULL, 
                       yloc = NULL)

mif$derived$bivariate_pair_correlation %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Correlation Permuted`, color = deidentified_sample)) +
  facet_grid(~From)

## -----------------------------------------------------------------------------
mif = interaction_variable(mif = mif,
                           mnames = markers[1:2],
                           r_range = 0:100,
                           num_permutations = 25,
                           keep_permutation_distribution = FALSE,
                           workers = 1,
                           overwrite = TRUE,
                           xloc = NULL,
                           yloc = NULL)

mif$derived$interaction_variable %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Interaction Permuted`, color = deidentified_sample)) +
  facet_grid(~From)

## -----------------------------------------------------------------------------
mif = dixons_s(mif = mif, 
               mnames = markers[1:2], 
               num_permutations = 25, 
               type = "Z", 
               workers = 1, 
               overwrite = TRUE, 
               xloc = NULL, 
               yloc = NULL)

mif$derived$Dixon_Z %>%
    filter(From != To) %>%
    ggplot() +
    geom_point(aes(x = Z, y = S, color = deidentified_sample)) +
    facet_grid(~From)

