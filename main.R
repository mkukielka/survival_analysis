source('models.R')
source('utils.R')
load.packages()

# register parralel computing
cl <- makeCluster(8, type = "FORK")
doParallel::registerDoParallel(cl)

# generate Cox smodels
results <- foreach(cancer_type = c('luad', 'lusc'), .combine=bind_rows) %:%
  foreach(stats_name = c('prevalence', 'microenvironment') , .combine=bind_rows) %:%
  foreach(add_pathological_stage = c(TRUE, FALSE), .combine=bind_rows) %:%
  foreach(add_smoking_data = c(TRUE, FALSE), .combine=bind_rows) %:%
  foreach(mutations_variant = c(0.05, 0.1), .combine=bind_rows) %dopar% {
    generate_cox_models(cancer_type, stats_name,
                        windowed_data = F, unified=T,
                        add_pathological_stage=add_pathological_stage,
                        add_smoking_data=add_smoking_data,
                        mutations_variant=mutations_variant)
  }
  
# save results as csv
results <- results[order(results$mean_cindex, decreasing = T),]
write.table(results, 'stats.tsv', sep='\t', row.names=FALSE)

# Save results to spreadsheet
openxlsx::write.xlsx(list("master table"= results), file = "stats.xlsx")
