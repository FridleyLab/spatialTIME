test_that("spatial_exp_to_mif() function", {
   # Create mif object
  
   ovarian <- VectraPolarisData::HumanOvarianCancerVP()
  
   ova_mif <- spatial_exp_to_mif(spatial_exp = ovarian,
                                 patient_id = "sample_id",
                                 markers = c("phenotype_cd68", 
                                             "phenotype_cd3", 
                                             "phenotype_cd8"))
   #TODO: Put tests that fit the requirements of what a MIF object should have
  
   spe_lung <- VectraPolarisData::HumanLungCancerV3()
  
   spe_mif <- spatial_exp_to_mif(spatial_exp = spe_lung,
                                 patient_id = "slide_id",
                                 sample_id = "slide_id",
                                 markers = c(
                                   "phenotype_cd4",
                                   "phenotype_cd8",
                                   "phenotype_cd14",
                                   "phenotype_cd19",
                                   "phenotype_ck",
                                   "phenotype_other"))
   
   #TODO: Put tests that fit the requirements of what a MIF object should have
   
})
