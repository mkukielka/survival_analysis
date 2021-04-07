source('constants.R')

load.packages <- function() {
  suppressPackageStartupMessages({
    library("tidyverse")
    library("Xmisc")
    library("survival")
    library("dplyr")
    library("survminer")
    library("mltools")
    library("tidyr")
    library("dplyr")
    library("stringr")
    library("caret")
    library("rlist")
    library("mefa")
    library("doParallel")
    library("openxlsx")
  })
}


load.dataset <- function(cancer_type, stats_type, add_pathological_stage = F,
                         windowed_data = F, add_smoking_data = F, unified = T,
                         mutations_variant = 0.05, high_low_classes = F) {
  cancer_type <- tolower(cancer_type)
  if (cancer_type %in% c('luad', 'lusc')) {
    if (unified) {
      load.unified.dataset(cancer_type, stats_type, windowed_data,
                           add_smoking_data = add_smoking_data,
                           add_pathological_stage = add_pathological_stage,
                           mutations_variant = mutations_variant,
                           high_low_classes = high_low_classes)
    } else {
      load.cancer.dataset(cancer_type, stats_type, windowed_data,
                          add_smoking_data = add_smoking_data,
                          add_pathological_stage = add_pathological_stage)
    }
  } else if (cancer_type == 'combined') {
    if (unified) {
      load.unified.combined(stats_type, windowed_data,
                            add_smoking_data = add_smoking_data,
                            add_pathological_stage = add_pathological_stage,
                            mutations_variant = mutations_variant,
                            high_low_classes = high_low_classes)
    } else {
      load.combined.dataset(stats_type, windowed_data,
                            add_smoking_data = add_smoking_data,
                            add_pathological_stage = add_pathological_stage)
    }
  } else {
    stop('Unsupported dataset! Use luad, luad or combined.')
  }
}


prepare.formula <- function(desired_variables) {
  variables.to.drop <- c('patient_id', 'time', 'vital_status')
  filtered.variables <- desired_variables[!desired_variables %in% variables.to.drop]
  variables.str <- paste0('`', str_c(filtered.variables, collapse = '` + `'), '`')
  as.formula(paste('Surv(time, vital_status)', variables.str, sep = " ~ "))
}


assign_age_group <- function(age) {
  if (age < 65) {
    return('lower than 65')
  } else {
    return('greater than 65')
  }
}


assign.smoking.status <- function(packs) {
  if (is.na(packs)) {
    return(NA)
  } else if (packs == 0) {
    return('non-smoker')
  } else if (packs <= 10) {
    return('light smoker (<10 packs annually)')
  } else {
    return('heavy smoker (>10 packs annually)')
  }
}

assign.pathologic.stage <- function(stage) {
  if (is.na(stage)) {
    return(NA)
  }
  stage.stripped <- rstrip(stage, 'ab')
  if (stage.stripped == 'stage i') {
    return('early (stage I)')
  } else if (stage.stripped == 'stage ii' | stage == 'stage iiia') {
    return('itermediate (stage II + IIIa)')
  } else if (stage.stripped == 'stage iv' | stage == 'stage iiib') {
    return('advanced (stage IIIB + IV)')
  } else if (stage.stripped == 'stage iii') {
    return('itermediate (stage II + IIIa)')
  }
}

get.stats.type <- function(stats_name, windowed_data) {
  stats_name <- tolower(stats_name)
  if (stats_name == 'prevalence') {
    'stats'
  } else if (stats_name == 'microenvironment') {
    'counts'
  } else {
    stop("Unknown statistics! Use prevalence or microenvironment.")
  }
}

append.cox.results <- function(results, cancer_type, model, cox_model, c_index_cutoff=0.65) {
  rbind(results, list(cancer_type = cancer_type,
                      model = model,
                      mean_cindex = cox_model$mean_cindex,
                      sd_cindex = cox_model$sd_cindex,
                      best_cindex = cox_model$best_cindex),
        stringsAsFactors = F)
}

determine.tissue.subsets <- function(cancer_type, stats_type) {
  if(tolower(cancer_type) == 'luad') {
    if(tolower(stats_type) == 'prevalence') {
      # luad prevalence
      return(list(c('TUMOR', 'VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'MIXED'),
                  c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA')))
    } else {
      # luad microenvironment
      return(list(c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'BRONCHI'),
                  c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA')))
    }
  } else if (tolower(cancer_type) == 'lusc') {
    if(tolower(stats_type) == 'prevalence') {
      # lusc prevalence
      return(list(c('TUMOR', 'VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA'), 
                  c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA')))
    } else {
      # lusc microenvironment
      return(list(c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'LUNG'),
                  c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA')))
    }
  } else {
    return(c(PREVALENCE.TISSUE.SUBSETS, MICROENVIRONMENT.TISSUE.SUBSETS))
  }
}
