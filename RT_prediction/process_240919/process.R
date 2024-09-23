library(magrittr)
# data_AP <- openxlsx::read.xlsx("RT_prediction/process_240919/20240817_dx-DEANS-AP+Hy_tR.xlsx", sheet = 1)
# data_Hy <- openxlsx::read.xlsx("RT_prediction/process_240919/20240817_dx-DEANS-AP+Hy_tR.xlsx", sheet = 5)
# data_Hy <- data_Hy[, -19]
# data <- rbind(data_AP, data_Hy)
# Library_rawAmido <- openxlsx::read.xlsx("Library/rawAmido.xlsx", sheet = 1)
# Library_rawPheHydro <- openxlsx::read.xlsx("Library/rawPheHydro.xlsx", sheet = 1)
# Library_rawHydroxyl <- openxlsx::read.xlsx("Library//rawHydroxyl.xlsx", sheet = 1) %>% 
#   dplyr::select(name, formula, smiles, monoisotop_mass, inchi, hmdb_id)
# Library_rawHydroxyl$inchi <- "HMDB"
# colnames(Library_rawHydroxyl) <- colnames(Library_rawAmido)
# Library <- rbind(Library_rawAmido, Library_rawPheHydro, Library_rawHydroxyl)
# data2 <- data %>%
#   dplyr::left_join(Library, by = c("HMDB" = "ID")) %>%
#   dplyr::select(Class, Name.x, `d0/min`, CAS, Classification, HMDB, KEGG, Pubchem, ntag, SMILES) %>% 
#   dplyr::distinct(HMDB, KEGG, Pubchem, .keep_all = TRUE)
# openxlsx::write.xlsx(data2, file = "RT_prediction/process_240919/data2.xlsx")

# data2 <- openxlsx::read.xlsx("RT_prediction/process_240919/data2.xlsx")
# reticulate::source_python("RT_prediction/process_240919/der_smiles.py")
# smiles <- data2$SMILES[211]
# smiles <- "C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N"
# smiles_der <- Der(smiles)
# 
# # load a smiles string
# atmcontainers <- purrr::map(
#   c(smiles, smiles_der),
#   depict::parse_smiles
# )
# many_containers <- depict::atomcontainer_list_to_jarray(atmcontainers)
# pic <- depict::depiction() |>
#   depict::set_zoom(2) |>
#   depict::depict(many_containers) |>
#   depict::get_svg_string()
# 
# svgtools::display_svg(pic, height = 500, width = 500)

data2 <- openxlsx::read.xlsx("RT_prediction/process_240919/data2.xlsx")
reticulate::source_python("RT_prediction/process_240919/der_smiles.py")
smilesDerVec <- as.character(sapply(1:nrow(data2), function(i) {
  class_tmp <- data2[i, ]$Class
  smiles_tmp <- data2[i, ]$SMILES
  if(class_tmp == "Hy"){
    smiles_der <- Der_Hy(smiles_tmp)
  }else{
    smiles_der <- Der_AP(smiles_tmp)
  }
  return(smiles_der)
}))
data2$SMILES_Der <- smilesDerVec
inchikeyDerVec <- as.character(sapply(data2$SMILES_Der, rinchi::get.inchi.key))
data2$INCHIKEY_Der <- inchikeyDerVec
openxlsx::write.xlsx(data2, file = "RT_prediction/process_240919/rawRT_info.xlsx")

atmcontainers <- purrr::map(
  data2$SMILES_Der,
  depict::parse_smiles
)
many_containers <- depict::atomcontainer_list_to_jarray(atmcontainers)
pic <- depict::depiction() |>
  depict::set_zoom(2) |>
  depict::depict(many_containers) |>
  depict::get_svg_string()
write(pic, "RT_prediction/process_240919/pic.svg")
svgtools::display_svg(pic)
