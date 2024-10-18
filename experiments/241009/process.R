library(magrittr)
# report_Hy_mix1 <- openxlsx::read.xlsx("./experiments/241009/report_Hy_mix1.xlsx", sheet = 1)
# report_Hy_mix2 <- openxlsx::read.xlsx("./experiments/241009/report_Hy_mix2.xlsx", sheet = 1)
# report_Hy_mix3 <- openxlsx::read.xlsx("./experiments/241009/report_Hy_mix3.xlsx", sheet = 1)
# report_Hy_mix4 <- openxlsx::read.xlsx("./experiments/241009/report_Hy_mix4.xlsx", sheet = 1)
# report_Hy_mix5 <- openxlsx::read.xlsx("./experiments/241009/report_Hy_mix5.xlsx", sheet = 1)
# report_Hy <- rbind(report_Hy_mix1, report_Hy_mix2, report_Hy_mix3, report_Hy_mix4, report_Hy_mix5) %>% 
#   dplyr::arrange(Mass, Rt)
# nrow(report_Hy) # 89
# openxlsx::write.xlsx(report_Hy, file = "./experiments/241009/report_Hy.xlsx")

report_Hy <- openxlsx::read.xlsx("./experiments/241009/report_Hy.xlsx", sheet = 1) # 82
# mass = 124.0901289
res1 <- pubcmpR::search_pubchem_online(inputs = pubcmpR::massSearch_CID(mass = 124.0901289, tol_mass = 0.001))
res1$cid
# mass = 256.2832452
res2 <- pubcmpR::search_pubchem_online(inputs = pubcmpR::massSearch_CID(mass = 256.2832452, tol_mass = 0.001))
res2$cid
# mass = 284.3145817
res3 <- pubcmpR::search_pubchem_online(inputs = pubcmpR::massSearch_CID(mass = 284.3145817, tol_mass = 0.001))
res3$cid

# false positive
raw_table <- openxlsx::read.xlsx("./RT_prediction/process_240919/20240922_rawRT_info.xlsx", sheet = 1)
length(which(sapply(report_Hy$ID, function(x) {
  any(strsplit(x, split = "\\|")[[1]] %in% raw_table$HMDB)
}))) # 64
61 / 82 # 74% true metabolites 看准确率不看假阳性 准确率加PubChem
61 / 73 # 84% 计算鉴定率时可以把PubChem的物质删去
report_Hy[which(sapply(report_Hy$ID, function(x) {
  any(strsplit(x, split = "\\|")[[1]] %in% raw_table$HMDB)
})), ]$ID
