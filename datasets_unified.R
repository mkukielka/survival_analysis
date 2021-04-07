load.unified.dataset <- function(cancer_type, stats_type, windowed_data = F,
                                 add_smoking_data = F, add_pathological_stage = F,
                                 mutations_variant = 0.05, high_low_classes=F) {
  # Clinical data
  cancer <- read.csv(paste0('TCGA/training_dataset_', tolower(cancer_type), '_stage.csv'), row.names = 1)
  cancer$sex <- ifelse(cancer$sex == 1, 'male', 'female')
  cancer$age <- ifelse(cancer$age == 1, 'age >= 65', 'age < 65')

  # Add time and vital status
  patients <- read.csv(paste0('TCGA/tcga_patient_data_', cancer_type, '.csv'),
                       row.names = 1, stringsAsFactors = FALSE)
  colnames(patients)[1] <- 'patient_id'
  patients$patient_id <- factor(patients$patient_id)
  time <- patients$days_to_death
  time.na <- which(is.na(time))
  time[time.na] <- patients$days_to_last_followup[time.na]
  patients$time <- time

  cancer <- merge(cancer, patients[, c('patient_id', 'vital_status', 'time')], by = 'patient_id')

  # Metadata
  metadata.cols <- c('patient_id', 'time', 'vital_status', 'age', 'sex')
  if (add_smoking_data) {
    metadata.cols <- c(metadata.cols, 'smoking_status')
  }
  if (add_pathological_stage) {
    metadata.cols <- c(metadata.cols, 'disease_stage')
  }
  metadata <- cancer[, metadata.cols]
  metadata <- metadata[complete.cases(metadata),]

  common_patients <- intersect(metadata$patient_id, patients$patient_id)

  if (mutations_variant == 0.05) {
    # mutations variant < 0.05
    if (tolower(cancer_type) == 'luad') {
      mutations.cols <- c('EGFR', 'STK11', 'TP53')
    } else {
      # lusc
      mutations.cols <- c('CDKN2A', 'NFE2L2', 'TP53')
    }
  } else if (mutations_variant == 0.1) {
    # mutations variant < 0.1
    if (tolower(cancer_type) == 'luad') {
      mutations.cols <- c('CDKN2A', 'EGFR', 'STK11', 'TP53', 'RET')
    } else {
      # lusc
      mutations.cols <- c('CDKN2A', 'KRAS', 'NFE2L2', 'RB1', 'TP53')
    }
  } else {
    # Mutations
    mutations.cols <- c('AKT1', 'ALK', 'BRAF', 'CD274', 'CDKN2A', 'DDR2', 'EGFR', 'ERBB2',
                        'FGFR1', 'FGFR3', 'KEAP1', 'KRAS', 'MAP2K1', 'MET', 'NF1', 'NFE2L2',
                        'NRAS', 'NTRK1', 'PIK3CA', 'PTEN', 'RB1', 'RET', 'ROS1', 'SMARCA4',
                        'STK11', 'TP53', 'U2AF1L4')
  }

  metadata.mutations.cols <- c(metadata.cols, mutations.cols)

  metadata.mutations <- cancer[cancer$patient_id %in% common_patients, metadata.mutations.cols]
  for (mutation in mutations.cols) {
    metadata.mutations[mutation] <- factor(metadata.mutations[, mutation, drop = T], levels = c(0, 1),
                                           labels = c('absent', 'present'))
  }
  # keep only desired columns
  metadata.mutations <- metadata.mutations[, metadata.mutations.cols]
  metadata.mutations <- metadata.mutations[complete.cases(metadata.mutations),]
  
  tissues <- TISSUES
  if(stats_type == 'stats') {
    tissues <- ALL.TISSUES
  }
  
  # Prevalence/Microenvironment counts
  if (windowed_data) {
    metadata.mutations.stats <- prepare.tissue.statistics.windowed(metadata.mutations, stats_type)
    final_patients <- unique(metadata.mutations.stats$patient_id)
    metadata <- metadata[metadata$patient_id %in% (final_patients),]
    metadata.mutations <- metadata.mutations[metadata.mutations$patient_id %in% (final_patients),]
  } else {
    stats.cols <- colnames(cancer)[grepl('percentage', colnames(cancer)) | grepl('count', colnames(cancer))]
    metadata.mutations.stats <- metadata.mutations
    metadata.mutations.stats[metadata.mutations.stats$patient_id %in% common_patients, stats.cols] <- cancer[
      cancer$patient_id %in% common_patients, stats.cols]

    already.processed.cols <- c('ITLR', 'Morisita_Horn')
    metadata.mutations.stats[metadata.mutations.stats$patient_id %in% common_patients, already.processed.cols] <-
      cancer[cancer$patient_id %in% common_patients, already.processed.cols]

    metadata.mutations.stats <- prepare.tissue.metrics(metadata.mutations.stats, stats_type)

    # keep only desired columns
    metadata.mutations.stats <- metadata.mutations.stats[, c(metadata.cols, mutations.cols, tissues, METRICS)]
    metadata.mutations.stats <- metadata.mutations.stats[complete.cases(metadata.mutations.stats),]
  }
  
  if(tolower(cancer_type) != 'combined' & high_low_classes) {
    for (tissue in c(tissues, METRICS)) {
      tryCatch({
        y <- surv_cutpoint(
          metadata.mutations.stats,
          time = "time",
          event = "vital_status",
          variables = c(tissue)
        )
        new.levels <- ifelse(metadata.mutations.stats[, c(tissue), drop=T] > y$cutpoint[[1]], 'WYSOKI', 'NISKI')
        metadata.mutations.stats[tissue] <- as.factor(new.levels)
      }, error = function(error_condition) {
        avg.tissue <- mean(metadata.mutations.stats[, c(tissue), drop = T])
        new.levels <- ifelse(metadata.mutations.stats[, c(tissue), drop=T] > avg.tissue, 'WYSOKI', 'NISKI')
        metadata.mutations.stats[tissue] <- as.factor(new.levels)
      })
    }
  }
  
  return(list(patients = metadata,
              patients.mutations = metadata.mutations,
              patients.mutations.stats = metadata.mutations.stats))
}


prepare.tissue.metrics <- function(data.stats, stats_type) {
  # add new empty columns for tissue statistics/metrics
  data.stats[, c('Lymphocyte_Ratio', 'Shannon', 'Simpson')] <- 0
  
  tissues <- TISSUES
  if(stats_type == 'stats') {
    tissues <- ALL.TISSUES
  }
  
  # Tissue prevalence/microenvironment
  stats.cols <- unlist(lapply(tissues, function(x) {
    paste0(x, ifelse(stats_type == 'stats', '_percentage', '_nb_percentage'))
  }))

  data.stats[, tissues] <- data.stats[, stats.cols]

  if (stats_type == 'stats') {
    for (i in 1:nrow(data.stats)) {

      # Tissue prevalence
      tissue.proportions <- data.stats[i, stats.cols]

      tumor <- data.stats[i, 'TUMOR_count']
      immune <- data.stats[i, 'IMMUNE_count']
      stroma <- data.stats[i, 'STROMA_count']

      # ITLR, DTLR
      #data.stats[i, 'ITLR'] <- data.stats[i, 'IMMUNE_nb_count'] / tumor
      #data.stats[i, 'DTLR'] <- data.stats[i, 'IMMUNE_nb_count'] / tumor

      # Lymphocyte Ratio
      data.stats[i, 'Lymphocyte_Ratio'] <- immune / (tumor + immune + stroma)

      # Morisita Horn
      #data.stats[i, 'Morisita_Horn'] <- (2 * tumor * immune) / ((tumor^2) + (immune^2))

      # Shannon, Simpson
      simpson <- 0
      shannon <- 0
      for (tissue in tissue.proportions) {
        if (tissue > 0) {
          # filter NA values
          shannon <- shannon + (tissue * log(tissue))
        }
        simpson <- simpson + (tissue^2)
      }
      data.stats[i, 'Shannon'] <- -shannon
      data.stats[i, 'Simpson'] <- simpson
    }
  }
  data.stats
}

load.unified.combined <- function(stats_type, windowed_data, 
                                  add_smoking_data = F,
                                  add_pathological_stage = F,
                                  mutations_variant = 0.05,
                                  high_low_classes = F) {
  # LUAD
  data.luad <- load.unified.dataset('luad', stats_type, windowed_data,
                                    add_smoking_data = add_smoking_data,
                                    add_pathological_stage = add_pathological_stage,
                                    mutations_variant = mutations_variant)
  data.luad$patients$cancer_type <- 'LUAD'
  data.luad$patients.mutations$cancer_type <- 'LUAD'
  data.luad$patients.mutations.stats$cancer_type <- 'LUAD'

  # LUSC
  data.lusc <- load.unified.dataset('lusc', stats_type, windowed_data,
                                    add_smoking_data = add_smoking_data,
                                    add_pathological_stage = add_pathological_stage,
                                    mutations_variant = mutations_variant)
  data.lusc$patients$cancer_type <- 'LUSC'
  data.lusc$patients.mutations$cancer_type <- 'LUSC'
  data.lusc$patients.mutations.stats$cancer_type <- 'LUSC'

  # clinical data
  patients <- rbind(data.luad$patients, data.lusc$patients)
  patients$cancer_type <- factor(patients$cancer_type)

  # clinical data + mutations
  common.cols <- intersect(colnames(data.lusc$patients.mutations),
                           colnames(data.luad$patients.mutations))
  patients.mutations <- rbind(data.lusc$patients.mutations[, common.cols],
                              data.luad$patients.mutations[, common.cols])
  patients.mutations$cancer_type <- factor(patients.mutations$cancer_type)

  # clinical + mutations + tissue statistics
  common.cols <- intersect(colnames(data.lusc$patients.mutations.stats),
                           colnames(data.luad$patients.mutations.stats))
  patients.mutations.stats <- rbind(data.lusc$patients.mutations.stats[, common.cols],
                                    data.luad$patients.mutations.stats[, common.cols])
  patients.mutations.stats$cancer_type <- factor(patients.mutations.stats$cancer_type)
  
  if(high_low_classes) {
    for (tissue in c(TISSUES, METRICS)) {
      tryCatch({
        y <- surv_cutpoint(
          metadata.mutations.stats,
          time = "time",
          event = "vital_status",
          variables = c(tissue)
        )
        new.levels <- ifelse(as.numeric(patients.mutations.stats[, c(tissue), drop=T]) > y$cutpoint[[1]], 'WYSOKI', 'NISKI')
        patients.mutations.stats[tissue] <- as.factor(new.levels)
      }, error = function(error_condition) {
        avg.tissue <- mean(as.numeric(patients.mutations.stats[, c(tissue), drop = T]))
        new.levels <- ifelse(as.numeric(patients.mutations.stats[, c(tissue), drop=T]) > avg.tissue, 'WYSOKI', 'NISKI')
        patients.mutations.stats[tissue] <- as.factor(new.levels)
      })
    }
  }
  return(list(patients = patients,
              patients.mutations = patients.mutations,
              patients.mutations.stats = patients.mutations.stats))
}
