library(magrittr)
# 处理AmidoLibrary
{
  rawAmido <- readr::read_tsv("./raw_data/AmidoLibrary.txt") # 15414
  rawAmido_HMDB <- rawAmido %>% dplyr::filter(Source == "HMDB")
  rawAmido_MoNA <- rawAmido %>% dplyr::filter(Source == "MoNA")
  rawAmido_T3DB <- rawAmido %>% dplyr::filter(Source == "T3DB")

  rawAmido_HMDB <- rawAmido_HMDB %>%
    dplyr::distinct(Formula, SMILES, .keep_all = TRUE)
  rawAmido_MoNA <- rawAmido_MoNA %>%
    dplyr::distinct(Formula, SMILES, .keep_all = TRUE)
  rawAmido_T3DB <- rawAmido_T3DB %>%
    dplyr::distinct(Formula, SMILES, .keep_all = TRUE)

  delete_MoNA_idx <- c()
  delete_T3DB_idx <- c()
  for(i in 1:nrow(rawAmido_HMDB)){
    print(i)
    tmp1 <- rawAmido_HMDB[i, ]
    tmp2 <- rawAmido_MoNA[which(rawAmido_MoNA$Formula == tmp1$Formula & rawAmido_MoNA$SMILES == tmp1$SMILES), ]
    tmp3 <- rawAmido_T3DB[which(rawAmido_T3DB$Formula == tmp1$Formula & rawAmido_T3DB$SMILES == tmp1$SMILES), ]
    if(length(tmp2) == 0 & length(tmp3) == 0){
      next
    }else{
      if(length(tmp2) != 0){
        idx <- which(rawAmido_MoNA$Formula == tmp1$Formula & rawAmido_MoNA$SMILES == tmp1$SMILES)
        delete_MoNA_idx <- c(delete_MoNA_idx, idx)
      }else if(length(tmp3) != 0){
        idx <- which(rawAmido_T3DB$Formula == tmp1$Formula & rawAmido_T3DB$SMILES == tmp1$SMILES)
        delete_T3DB_idx <- c(delete_T3DB_idx, idx)
      }
    }
  }
  delete_MoNA_idx <- unique(delete_MoNA_idx)
  delete_T3DB_idx # NULL
  rawAmido_MoNA <- rawAmido_MoNA[-delete_MoNA_idx, ]
  rawAmido <- rbind(rawAmido_HMDB, rawAmido_MoNA, rawAmido_T3DB) %>% # 14548
    dplyr::arrange(Mass)

  rawAmido <- rawAmido %>%
    dplyr::filter(!stringr::str_detect(Name, "PE[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PS[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PA[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PG[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PC[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PGP[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "DG[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PEP[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "[Cc]er[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PI[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "Phosphatidylethanolamine[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "Phosphatidylserine[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "Phytosphingosine[\\(\\s-]")) %>%
    dplyr::filter(Mass <= 1000)
  openxlsx::write.xlsx(rawAmido, "./raw_data/rawAmido.xlsx")
}

# 处理PheHydroLibrary
{
  rawPheHydro <- readr::read_tsv("./raw_data/PheHydroLibrary.txt")
  rawPheHydro_HMDB <- rawPheHydro %>% dplyr::filter(Source == "HMDB")
  rawPheHydro_MoNA <- rawPheHydro %>% dplyr::filter(Source == "MoNA")
  rawPheHydro_T3DB <- rawPheHydro %>% dplyr::filter(Source == "T3DB")

  rawPheHydro_HMDB <- rawPheHydro_HMDB %>%
    dplyr::distinct(Formula, SMILES, .keep_all = TRUE)
  rawPheHydro_MoNA <- rawPheHydro_MoNA %>%
    dplyr::distinct(Formula, SMILES, .keep_all = TRUE)
  rawPheHydro_T3DB <- rawPheHydro_T3DB %>%
    dplyr::distinct(Formula, SMILES, .keep_all = TRUE)

  delete_MoNA_idx <- c()
  delete_T3DB_idx <- c()
  for(i in 1:nrow(rawPheHydro_HMDB)){
    print(i)
    tmp1 <- rawPheHydro_HMDB[i, ]
    tmp2 <- rawPheHydro_MoNA[which(rawPheHydro_MoNA$Formula == tmp1$Formula & rawPheHydro_MoNA$SMILES == tmp1$SMILES), ]
    tmp3 <- rawPheHydro_T3DB[which(rawPheHydro_T3DB$Formula == tmp1$Formula & rawPheHydro_T3DB$SMILES == tmp1$SMILES), ]
    if(length(tmp2) == 0 & length(tmp3) == 0){
      next
    }else{
      if(length(tmp2) != 0){
        idx <- which(rawPheHydro_MoNA$Formula == tmp1$Formula & rawPheHydro_MoNA$SMILES == tmp1$SMILES)
        delete_MoNA_idx <- c(delete_MoNA_idx, idx)
      }else if(length(tmp3) != 0){
        idx <- which(rawPheHydro_T3DB$Formula == tmp1$Formula & rawPheHydro_T3DB$SMILES == tmp1$SMILES)
        delete_T3DB_idx <- c(delete_T3DB_idx, idx)
      }
    }
  }
  delete_MoNA_idx <- unique(delete_MoNA_idx)
  delete_T3DB_idx # NULL
  rawPheHydro_MoNA <- rawPheHydro_MoNA[-delete_MoNA_idx, ]
  rawPheHydro <- rbind(rawPheHydro_HMDB, rawPheHydro_MoNA, rawPheHydro_T3DB) %>% # 12537
    dplyr::arrange(Mass)

  rawPheHydro <- rawPheHydro %>%
    dplyr::filter(!stringr::str_detect(Name, "PE[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PS[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PA[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PG[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PC[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PGP[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "DG[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PEP[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "[Cc]er[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "PI[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "Phosphatidylethanolamine[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "Phosphatidylserine[\\(\\s-]")) %>%
    dplyr::filter(!stringr::str_detect(Name, "Phytosphingosine[\\(\\s-]")) %>%
    dplyr::filter(Mass <= 1000)
  openxlsx::write.xlsx(rawPheHydro, file = "./raw_data/rawPheHydro.xlsx")
}
