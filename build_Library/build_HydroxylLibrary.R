# devtools::install_github("dongdongdong7/pubcmpR")
library(magrittr)
hmdbCmpTb <- pubcmpR::load_hmdbCmpTb()
reticulate::source_python("build_Library/whetherFunction.py")
reticulate::source_python("RT_prediction/process_240919/der_smiles.py")

logical_vec <- sapply(1:nrow(hmdbCmpTb), function(i) {
  print(paste0(i, " / ", "217920"))
  smiles <- hmdbCmpTb[i, ]$smiles
  if(is.na(smiles)) return(FALSE)
  else{
    res <- tryCatch(whether_COH_(smiles), 
                    error = function(e) {
                      print(e)
                      FALSE
                    })
    return(res)
  }
})

hmdbCmpTb_Hydroxyl <- hmdbCmpTb[which(logical_vec), ]
# remove lipid and mass > 1000
hmdbCmpTb_Hydroxyl <- hmdbCmpTb_Hydroxyl %>%
  dplyr::filter(monisotopic_molecular_weight < 1000)

# kingdom_vec <- sapply(hmdbCmpTb_Hydroxyl$taxonomy, function(x) {x$kingdom})
# super_class_vec <- sapply(hmdbCmpTb_Hydroxyl$taxonomy, function(x) {x$super_class})
# 
# hmdbCmpTb_Hydroxyl <- hmdbCmpTb_Hydroxyl[-c(which(kingdom_vec == "Inorganic compounds"),
#                                       which(super_class_vec == "Lipids and lipid-like molecules")), ]

HydroxylLibrary <- hmdbCmpTb_Hydroxyl %>% 
  dplyr::select(accession, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey)
smiles_der_vec <- sapply(HydroxylLibrary$smiles, function(x) {
  Der_Hy(x)
})
inchi_der_vec <- sapply(smiles_der_vec, function(x) {
  res <- tryCatch(rinchi::get.inchi(x), 
                  error = function(e) {
                    print(e)
                    NA
                  })
})
inchikey_der_vec <- sapply(smiles_der_vec, function(x) {
  res <- tryCatch(rinchi::get.inchi.key(x),
                  error = function(e) {
                    print(e)
                    NA
                  })
})
HydroxylLibrary$smiles_der <- smiles_der_vec
HydroxylLibrary$inchi_der <- inchi_der_vec
HydroxylLibrary$inchikey_der <- inchikey_der_vec
HydroxylLibrary <- HydroxylLibrary %>% 
  dplyr::filter(!is.na(inchi_der) & !is.na(inchikey_der))
#saveRDS(HydroxylLibrary, "./build_Library/HydroxylLibrary.rds")

# Predict RT
Retip::prep.wizard()
xgb <- readRDS("./xgb.rds")
training <- readRDS("./training.rds")
HydroxylLibrary <- readRDS("./HydroxylLibrary.rds")
rp_ext <- HydroxylLibrary %>% 
  dplyr::select(accession, inchikey_der, smiles_der)
colnames(rp_ext) <- c("NAME", "INCHKEY", "SMILES")
rp_ext_desc <- Retip::getCD(rp_ext)
#saveRDS(rp_ext_desc, "./rp_ext_desc_Hydroxyl.rds")
db_rt_ext <- Retip::proc.data(rp_ext_desc)
rp_ext_pred_xgb <- Retip::RT.spell(training, rp_ext_desc, model = xgb)

HydroxylLibrary <- HydroxylLibrary %>% 
  dplyr::left_join(rp_ext_pred_xgb, by = c("accession" = "NAME")) %>% 
  dplyr::select(accession, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey, smiles_der, inchi_der, inchikey_der, RTP)

openxlsx::write.xlsx(HydroxylLibrary, "./HydroxylLibrary.xlsx")
