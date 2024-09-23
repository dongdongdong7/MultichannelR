library(magrittr)
reticulate::source_python("library/hmdb/detect_alcohol.py")
data("metabolitesList", package = "MetaboRich")
data("setsList", package = "MetaboRich")
test_smiles <- metabolitesList$hmdb[4, ]$smiles
whether_COH(smiles = test_smiles)
pb <- utils::txtProgressBar(max = nrow(metabolitesList$hmdb), style = 3)
COH_logicalVec <- sapply(1:nrow(metabolitesList$hmdb), function(i) {
  utils::setTxtProgressBar(pb, i)
  smiles <- metabolitesList$hmdb[i, ]$smiles
  tryCatch({
    whether_COH(smiles = smiles)
  }, error = function(e) {
    return(FALSE)
  })
})

hmdb_alcohol <- metabolitesList$hmdb[COH_logicalVec, ] %>%
  dplyr::select(hmdb_id, name, formula, monoisotop_mass, smiles, inchi, inchikey) %>%
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
  dplyr::filter(!stringr::str_detect(name, "Phytosphingosine[\\(\\s-]")) %>%
  dplyr::filter(monoisotop_mass <= 1000)
nrow(hmdb_alcohol)
openxlsx::write.xlsx(hmdb_alcohol, file = "hmdb_alcohol.xlsx")
