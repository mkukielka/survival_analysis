# Survival Analysis

Survival analysis for lung cancer patients based on clinical data, 
histopathological images and genomic mutations based on data 
originating from TCGA project. Two cancer types are taken into consideration: 
lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC).

### Project structure:
* main.R - perform whole analysis with params variations
* datasets_unified.R - load and prepare clinical, mutations and tissue data
* models.R - generate cox models
* utils.R - utilities used across analysis
* constants.R - stores names of used tissue types and spatial metrics
* master_plots.R - script, which generate plots used in master thesis
