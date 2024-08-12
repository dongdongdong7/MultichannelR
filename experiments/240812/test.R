standard <- openxlsx::read.xlsx("D:/fudan/Projects/2024/MultichannelR/Progress/build_package/MultichannelR/experiments/240808/20240808_AP+Hy_HMDB信息.xlsx", sheet = 2)
standard_id <- standard$X8[!is.na(standard$X8)][-1]
report <- openxlsx::read.xlsx("D:/fudan/Projects/2024/MultichannelR/Progress/build_package/MultichannelR/experiments/240812/report.xlsx")

report_ex <- report[sapply(report$ID, function(x) {
  if(is.na(x)) return(FALSE)
  else return(any(stringr::str_detect(standard_id, x)))
}), ]

openxlsx::write.xlsx(report_ex, file = "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/MultichannelR/experiments/240812/report_ex.xlsx")
