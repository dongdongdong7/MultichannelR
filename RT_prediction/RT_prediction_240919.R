library(Retip)
library(magrittr)
prep.wizard()
rawRT_info <- openxlsx::read.xlsx("RT_prediction/rawRT_info.xlsx", sheet = 1)
data <- rawRT_info %>% 
  dplyr::select(Name.x, INCHIKEY_Der, SMILES_Der, `d0/min`)
colnames(data) <- c("NAME", "INCHKEY", "SMILES", "RT")
descs <- getCD(data)
db_rt <- proc.data(descs)
db <- db_rt[, -1]
descrCor <-  cor(db)
highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = .98)
filteredDescr <- db[,-highlyCorDescr]
db_rt <- db_rt[, c("RT", colnames(filteredDescr))]

set.seed(101)
in_training <- caret::createDataPartition(db_rt$XLogP, p = .8, list = FALSE)
training <- db_rt[in_training, ]
testing <- db_rt[-in_training, ]
rf  <- fit.rf(training)
lightgbm <- fit.lightgbm(training, testing)
xgb <- fit.xgboost(training)
stat <- get.score(testing,rf, lightgbm, xgb)
stat
p.model(testing, m = rf, title = "Random Forest")

t1 <- as.matrix(testing)
ncolt1 <- ncol(t1)
prd <- stats::predict(xgb, t1[, 2:ncolt1])
prd <- data.frame(prd)
names(prd) <- c("RTP")
rt_obs <- testing$RT
rt_pred <- prd$RTP
df <- data.frame(rt_obs, rt_pred, rt_obs - rt_pred)
names(df) <- c("RT_obs", "RT_pred", "RT_ERR")
res_rf_fh <- data.frame(round((caret::postResample(testing$RT, 
                                                   prd$RTP)), 2))
names(res_rf_fh) <- c("Stats")
stats_table <- cbind(Stats = c("RMSE", "R2", "MAE"), res_rf_fh$Stats)
stats_table <- gridExtra::tableGrob(stats_table, rows = NULL, 
                                    cols = NULL, theme = gridExtra::ttheme_minimal(core = list(fg_params = list(hjust = 0, 
                                                                                                                x = 0.1)), base_colour = "#384049"))
stats_table <- gtable::gtable_add_grob(stats_table, grobs = grid::rectGrob(gp = grid::gpar(fill = NA, 
                                                                                           lwd = 2)), t = 1, b = nrow(stats_table), l = 1, r = ncol(stats_table))
ggplot2::ggplot(df, aes(RT_obs, RT_pred)) + 
  ggplot2::geom_point(aes(RT_obs, RT_pred), colour = "#5E676F") + 
  ggplot2::labs(title = paste0("Predicted RT vs Real - xgb"), x = "Observed RT", y = "Predicted RT") + 
  ggplot2::xlim(0,  16) + ggplot2::ylim(0, 16) + 
  ggplot2::annotation_custom(stats_table, xmin = 1, xmax = 2, ymin = 16 - 4, ymax = 16 - 1) + 
  ggplot2::theme_classic() + 
  ggplot2::theme(plot.title = element_text(color = "#384049", face = "bold", hjust = 0.5), axis.line = element_line(colour = "#384049"), axis.text = element_text(colour = "#384049"), axis.title = element_text(colour = "#384049")) + 
  ggplot2::geom_abline(intercept = 0, 
                       slope = 1, color = "#D02937")
p.model.features(xgb, mdl_type = "xgb")

saveRDS(training, file = "./RT_prediction/training.rds")
saveRDS(xgb, file = "./RT_prediction/xgb.rds")
