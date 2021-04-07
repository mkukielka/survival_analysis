ALL.TISSUES <- c('TUMOR', 'STROMA', 'MIXED', 'IMMUNE', 'VESSEL', 'BRONCHI', 'NECROSIS', 'LUNG')
TISSUES <- c('STROMA', 'MIXED', 'IMMUNE', 'VESSEL', 'BRONCHI', 'NECROSIS', 'LUNG')
METRICS <- c('ITLR', 'Morisita_Horn', 'Shannon', 'Simpson', 'Lymphocyte_Ratio')

MICROENVIRONMENT.TISSUE.SUBSETS <- list(
  c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA'),
  c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'BRONCHI'), # luad micro
  c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'LUNG') # lusc micro
)


PREVALENCE.TISSUE.SUBSETS <- list(
  c('TUMOR', 'VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'MIXED'), # luad prevalence
  c('TUMOR', 'VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA') # lusc prevalence
)