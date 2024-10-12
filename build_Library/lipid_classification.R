hmdbCmpTb <- pubcmpR::load_hmdbCmpTb()

kingdom_vec <- sapply(hmdbCmpTb$taxonomy, function(x) {x$kingdom})
super_class_vec <- sapply(hmdbCmpTb$taxonomy, function(x) {x$super_class})
class_vec <- sapply(hmdbCmpTb$taxonomy, function(x) {x$class})
sub_class_vec <- sapply(hmdbCmpTb$taxonomy, function(x) {x$sub_class})

classification <- list()
classification$`Steroids and steroid derivatives` <- unique(sub_class_vec[which(class_vec == "Steroids and steroid derivatives")])
classification$`Steroids and steroid derivatives` <- classification$`Steroids and steroid derivatives`[!is.na(classification$`Steroids and steroid derivatives`)]
classification$`Fatty Acyls` <- unique(sub_class_vec[which(class_vec == "Fatty Acyls")])
classification$Glycerophospholipids <- unique(sub_class_vec[which(class_vec == "Glycerophospholipids")])
classification$Glycerophospholipids <- classification$Glycerophospholipids[!is.na(classification$Glycerophospholipids)]
classification$`Prenol lipids` <- unique(sub_class_vec[which(class_vec == "Prenol lipids")])
classification$Glycerolipids <- unique(sub_class_vec[which(class_vec == "Glycerolipids")])[!is.na(unique(sub_class_vec[which(class_vec == "Glycerolipids")]))]
classificationTb <- tibble::tibble(class = c(purrr::list_c(lapply(names(classification), function(x) {c(rep(x, length(classification[[x]])))})), "Endocannabinoids", "Saccharolipids"), 
                                   sub_class = c(purrr::list_c(classification), NA, NA))
#openxlsx::write.xlsx(classificationTb, file = "./build_Library/lipids_classification.xlsx")
