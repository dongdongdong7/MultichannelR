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
    res <- tryCatch(whether_NH12_(smiles), 
                    error = function(e) {
                      print(e)
                      FALSE
                    })
    return(res)
  }
})
hmdbCmpTb_Amido <- hmdbCmpTb[which(logical_vec), ]
# remove lipid and mass > 1000
hmdbCmpTb_Amido <- hmdbCmpTb_Amido %>%
  dplyr::filter(monisotopic_molecular_weight < 1000) %>% 
  dplyr::filter(!stringr::str_detect(name, "PE[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PS[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PA[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PG[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PC[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PGP[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "DG[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "MG[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PEP[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "[Cc]er[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PI[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PIP[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "PIP2[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "SM[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "CL[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "Phosphatidylethanolamine[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "Phosphatidylserine[\\(\\s-]")) %>%
  dplyr::filter(!stringr::str_detect(name, "Phytosphingosine[\\(\\s-]"))

kingdom_vec <- sapply(hmdbCmpTb_Amido$taxonomy, function(x) {x$kingdom})
super_class_vec <- sapply(hmdbCmpTb_Amido$taxonomy, function(x) {x$super_class})

hmdbCmpTb_Amido <- hmdbCmpTb_Amido[-c(which(kingdom_vec == "Inorganic compounds"),
                                      which(super_class_vec == "Lipids and lipid-like molecules")), ]

# which(kingdom_vec == "Inorganic compounds")
# which(super_class_vec == "Lipids and lipid-like molecules")
# which(is.na(super_class_vec))

AmidoLibrary <- hmdbCmpTb_Amido %>% 
  dplyr::select(accession, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey)
smiles_der_vec <- sapply(AmidoLibrary$smiles, function(x) {
  Der_AP(x)
})
inchi_der_vec <- sapply(smiles_der_vec, rinchi::get.inchi)
inchikey_der_vec <- sapply(smiles_der_vec, rinchi::get.inchi.key)
AmidoLibrary$smiles_der <- smiles_der_vec
AmidoLibrary$inchi_der <- inchi_der_vec
AmidoLibrary$inchikey_der <- inchikey_der_vec
#saveRDS(AmidoLibrary, "./AmidoLibrary.rds")

# Predict RT
Retip::prep.wizard()
xgb <- readRDS("./xgb.rds")
training <- readRDS("./training.rds")
AmidoLibrary <- readRDS("./AmidoLibrary.rds")
rp_ext <- AmidoLibrary %>% 
  dplyr::select(accession, inchikey_der, smiles_der)
colnames(rp_ext) <- c("NAME", "INCHKEY", "SMILES")
rp_ext_desc <- Retip::getCD(rp_ext)
# rp_ext_desc <- readRDS("./rp_ext_desc_Amido.rds")
# saveRDS(rp_ext_desc, "./rp_ext_desc_Amido.rds")
db_rt_ext <- Retip::proc.data(rp_ext_desc)
rp_ext_pred_xgb <- Retip::RT.spell(training, rp_ext_desc, model = xgb)

AmidoLibrary <- AmidoLibrary %>% 
  dplyr::left_join(rp_ext_pred_xgb, by = c("accession" = "NAME")) %>% 
  dplyr::select(accession, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey, smiles_der, inchi_der, inchikey_der, RTP)

openxlsx::write.xlsx(AmidoLibrary, "./AmidoLibrary.xlsx")
