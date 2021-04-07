source('constants.R')
source('datasets.R')
source('datasets_unified.R')
source('utils.R')
source('plots.R')


generate_cox_models <- function(
  cancer_type, stats_name, windowed_data = FALSE, add_smoking_data = F,
  mutations_variant=0.05, add_pathological_stage = F, all_mutations = F, 
  high_low_classes = F, k = 10, verbose = F, unified = T, custom_seed=34) {
  # set seed to obtain reproducible results
  set.seed(custom_seed)

  # prepare results df
  results <- data.frame(matrix(0, nrow = 0, ncol = 5))
  colnames(results) <- c('cancer_type', 'model', 'mean_cindex',
                         'sd_cindex', 'best_cindex')

  # assign stats type based on stats name
  stats_type <- get.stats.type(stats_name, windowed_data)

  # Load desired dataset
  data <- load.dataset(cancer_type, stats_type, windowed_data,
                       unified = unified,
                       add_smoking_data = add_smoking_data,
                       add_pathological_stage = add_pathological_stage,
                       mutations_variant=mutations_variant,
                       high_low_classes=high_low_classes)
  patients <- data$patients
  patients.mutations <- data$patients.mutations
  patients.mutations.stats <- data$patients.mutations.stats

  # based on tissue statistics and cancer type determine appropriate tissue subset
  tissues <- TISSUES
  if(tolower(stats_name) == 'prevalence') {
    tissues <- ALL.TISSUES
  }
  tissue.subsets <- determine.tissue.subsets(cancer_type, stats_name) 
  
  # If windowed data is used, add prefix '-windowed' to cancer type
  if (windowed_data) {
    cancer_type <- paste0(cancer_type, '-windowed')
  }

  if (!windowed_data & stats_name == 'prevalence') {
    # Metadata/Clinical data
    if (verbose) {
      print(paste('Generating model:', cancer_type, 'for clinical'))
    }
    metadata.formula <- prepare.formula(colnames(patients))
    metadata.cox <- train.cox.model.counts(cox_data = patients,
                                           cox_formula = metadata.formula,
                                           k = k,
                                           verbose = verbose)
    # append results
    results <- append.cox.results(results = results,
                                  cancer_type = cancer_type,
                                  model = 'clinical',
                                  cox_model = metadata.cox)

    # Metadata + mutations
    if (verbose) {
      print(paste('Generating model:', cancer_type, 'for clinical + mutations'))
    }
    mutations.formula <- prepare.formula(colnames(patients.mutations))
    mutations.cox <- train.cox.model.counts(cox_data = patients.mutations,
                                            cox_formula = mutations.formula,
                                            k = k,
                                            verbose = verbose)
    # append results
    results <- append.cox.results(results = results,
                                  cancer_type = cancer_type,
                                  model = 'clinical + mutations',
                                  cox_model = mutations.cox)

  }

  # Counts/Stats + survival
  if (verbose) {
    print(paste('Generating model:', cancer_type, 'for', stats_name))
  }
  stats.formula <- prepare.formula(tissues)
  stats.cox <- train.cox.model.counts(cox_data = patients.mutations.stats,
                                      cox_formula = stats.formula,
                                      k = k,
                                      verbose = verbose)

  # append results
  results <- append.cox.results(results = results,
                                cancer_type = cancer_type,
                                model = stats_name,
                                cox_model = stats.cox)

  # Counts/Stats + metadata
  if (verbose) {
    print(paste('Generating model:', cancer_type, 'for clinical +', stats_name))
  }
  metadata.stats.formula <- prepare.formula(c(colnames(patients), tissues))
  metadata.stats.cox <- train.cox.model.counts(cox_data = patients.mutations.stats,
                                               cox_formula = metadata.stats.formula,
                                               k = k,
                                               verbose = verbose)

  # append results
  results <- append.cox.results(results = results,
                                cancer_type = cancer_type,
                                model = paste0('clinical + ', stats_name),
                                cox_model = metadata.stats.cox)

  for (tissue.subset in tissue.subsets) {
    metadata.stats.formula <- prepare.formula(c(colnames(patients), tissue.subset))
    metadata.stats.cox <- train.cox.model.counts(cox_data = patients.mutations.stats,
                                                 cox_formula = metadata.stats.formula,
                                                 k = k,
                                                 verbose = verbose)
    # append results
    results <- append.cox.results(results = results,
                                  cancer_type = cancer_type,
                                  model = paste0('clinical + ', stats_name, ' (', paste(tissue.subset, collapse = ', '), ')'),
                                  cox_model = metadata.stats.cox)
  }

  # Metadata + mutations + tissue statistics
  if (verbose) {
    print(paste0('Generating model: ', cancer_type, ' for clinical + mutations + ', stats_name))
  }
  metadata.mutations.stats.formula <- prepare.formula(c(colnames(patients.mutations), tissues))
  all.cox <- train.cox.model.counts(cox_data = patients.mutations.stats,
                                    cox_formula = metadata.mutations.stats.formula,
                                    k = k,
                                    verbose = verbose)

  # append results
  results <- append.cox.results(results = results,
                                cancer_type = cancer_type,
                                model = paste0('clinical + mutations + ', stats_name),
                                cox_model = all.cox)
  
  for (tissue.subset in tissue.subsets) {
    metadata.mutations.stats.formula <- prepare.formula(c(colnames(patients.mutations), tissues))
    all.cox <- train.cox.model.counts(cox_data = patients.mutations.stats,
                                      cox_formula = metadata.mutations.stats.formula,
                                      k = k,
                                      verbose = verbose)
    
    # append results
    results <- append.cox.results(results = results,
                                  cancer_type = cancer_type,
                                  model = paste0('clinical + mutations + ', stats_name, ' (', paste(tissue.subset, collapse = ', '), ')'),
                                  cox_model = all.cox)
  }

  # Spatial metrics are only calculated for tissue prevalence
  if (stats_type == 'stats') {
    # Metadata + spatial metrics
    for (metric in METRICS) {
      if (verbose) {
        print(paste('Generating model:', cancer_type, 'for clinical +', metric))
      }

      patients.mutations.stats.tmp <- patients.mutations.stats
      if (windowed_data & metric == 'Morisita_Horn') {
        # since Morisita Horn aggregates data from slides, remove duplicated patients
        patients.mutations.stats.tmp <- patients.mutations.stats.tmp[
          !duplicated(patients.mutations.stats.tmp$patient_id),]
      }
      patients.metrics.formula <- prepare.formula(c(colnames(patients), metric))
      metrics.cox <- train.cox.model.counts(cox_data = patients.mutations.stats.tmp,
                                            cox_formula = patients.metrics.formula,
                                            k = k,
                                            verbose = verbose)

      # append results
      results <- append.cox.results(results = results,
                                    cancer_type = cancer_type,
                                    model = paste('clinical +', metric),
                                    cox_model = metrics.cox)

      # Metadata + mutations + spatial metrics
      if (verbose) {
        print(paste('Generating model', cancer_type, 'for clinical + mutations +', metric))
      }

      patients.mutations.metrics.formula <- prepare.formula(c(colnames(patients.mutations), metric))
      metrics.all.cox <- train.cox.model.counts(cox_data = patients.mutations.stats.tmp,
                                                cox_formula = patients.mutations.metrics.formula,
                                                k = k,
                                                verbose = verbose)

      # append results
      results <- append.cox.results(results = results,
                                    cancer_type = cancer_type,
                                    model = paste0('clinical + mutations + ', metric),
                                    cox_model = metrics.all.cox)
    }
  }
  results['mutations'] = paste('p-value <', mutations_variant)
  results['dataset'] = paste0(
    ifelse(!add_pathological_stage & !add_smoking_data, 'base', ''),
    ifelse(add_pathological_stage, 'p', ''),
    ifelse(add_smoking_data, 's', ''))
  results
}


train.cox.model.counts <- function(cox_data, cox_formula, k = 10, verbose = TRUE) {
  all_c_index <- c()
  for(i in 1:k) {
    # Cross validation, split data into k separated subsets (train + test)
    ids <- levels(cox_data$patient_id)
    data_folds <- createFolds(ids, k = k, list = TRUE, returnTrain = FALSE)
    
    for (current_fold in data_folds) {
      # subset data to train and test
      test_fold <- ids[current_fold]
      train <- cox_data[!cox_data$patient_id %in% test_fold, , drop = F]
      test <- cox_data[cox_data$patient_id %in% test_fold, , drop = F]
      
      # fit cox model
      generated_cox_model <- coxph(formula = cox_formula, data = train)
      
      # C-index of predictions
      current_c_index <- concordance(generated_cox_model, newdata = test)$concordance
      all_c_index <- c(all_c_index, round(current_c_index, 3))
    }  
  } 
  mean_c_index <- round(mean(all_c_index), 2)
  sd_c_index <- round(sd(all_c_index), 2)
  best_c_index <- round(max(all_c_index), 2)

  if (verbose) {
    print(paste('Mean C-index:', mean_c_index, '±', sd_c_index))
  }

  return(list(mean_cindex = mean_c_index,
              sd_cindex = sd_c_index,
              best_cindex = best_c_index,
              cox_data = cox_data,
              surv_formula = cox_formula)
  )
}

