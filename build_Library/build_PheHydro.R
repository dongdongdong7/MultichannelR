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
    res <- tryCatch(whether_cOH_(smiles), 
                    error = function(e) {
                      print(e)
                      FALSE
                    })
    return(res)
  }
})

hmdbCmpTb_PheHydro <- hmdbCmpTb[which(logical_vec), ]
# remove lipid and mass > 1000
hmdbCmpTb_PheHydro <- hmdbCmpTb_PheHydro %>%
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

kingdom_vec <- sapply(hmdbCmpTb_PheHydro$taxonomy, function(x) {x$kingdom})
super_class_vec <- sapply(hmdbCmpTb_PheHydro$taxonomy, function(x) {x$super_class})

hmdbCmpTb_PheHydro <- hmdbCmpTb_PheHydro[-c(which(kingdom_vec == "Inorganic compounds"),
                                      which(super_class_vec == "Lipids and lipid-like molecules")), ]

PheHydroLibrary <- hmdbCmpTb_PheHydro %>% 
  dplyr::select(accession, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey)
smiles_der_vec <- sapply(PheHydroLibrary$smiles, function(x) {
  Der_AP(x)
})
inchi_der_vec <- sapply(smiles_der_vec, rinchi::get.inchi)
inchikey_der_vec <- sapply(smiles_der_vec, rinchi::get.inchi.key)
PheHydroLibrary$smiles_der <- smiles_der_vec
PheHydroLibrary$inchi_der <- inchi_der_vec
PheHydroLibrary$inchikey_der <- inchikey_der_vec
saveRDS(PheHydroLibrary, "./build_Library/PheHydroLibrary.rds")

# Predict RT
Retip::prep.wizard()
xgb <- readRDS("./xgb.rds")
training <- readRDS("./training.rds")
PheHydroLibrary <- readRDS("./PheHydroLibrary.rds")
rp_ext <- PheHydroLibrary %>% 
  dplyr::select(accession, inchikey_der, smiles_der)
colnames(rp_ext) <- c("NAME", "INCHKEY", "SMILES")
rp_ext_desc <- Retip::getCD(rp_ext)
# rp_ext_desc <- readRDS("./rp_ext_desc_PheHydro.rds")
#saveRDS(rp_ext_desc, "./rp_ext_desc_PheHydro.rds")
db_rt_ext <- Retip::proc.data(rp_ext_desc)
rp_ext_pred_xgb <- Retip::RT.spell(training, rp_ext_desc, model = xgb)

PheHydroLibrary <- PheHydroLibrary %>% 
  dplyr::left_join(rp_ext_pred_xgb, by = c("accession" = "NAME")) %>% 
  dplyr::select(accession, name, chemical_formula, monisotopic_molecular_weight, smiles, inchi, inchikey, smiles_der, inchi_der, inchikey_der, RTP)

openxlsx::write.xlsx(PheHydroLibrary, "./PheHydroLibrary.xlsx")
