standard <- openxlsx::read.xlsx("./20240808_AP+Hy_HMDB信息.xlsx", sheet = 2)
standard_id <- standard$X8[!is.na(standard$X8)][-1]
report <- openxlsx::read.xlsx("./report_all.xlsx")

report_ex <- report[sapply(report$ID, function(x) {
  if(is.na(x)) return(FALSE)
  else return(any(stringr::str_detect(standard_id, x)))
}), ]

maeVec <- sapply(1:nrow(report_ex), function(i) {
  values <- report_ex[i, 8:25][!is.na(report_ex[i, 8:25])]
  ratios <- outer(values, values, "/")
  ratios_upper <- ratios[upper.tri(ratios)]
  differences <- ratios_upper - 1
  mae <- mean(abs(differences))
  return(mae)
})
report_ex$mae <- maeVec

rmseVec <- sapply(1:nrow(report_ex), function(i) {
  values <- report_ex[i, 8:25][!is.na(report_ex[i, 8:25])]
  ratios <- outer(values, values, "/")
  ratios_upper <- ratios[upper.tri(ratios)]
  differences <- ratios_upper - 1
  rmse <- sqrt(mean(differences^2))
  return(rmse)
})
report_ex$rmse <- rmseVec
openxlsx::write.xlsx(report_ex, file = "report_ex.xlsx")
