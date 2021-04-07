source('models.R')
source('constants.R')
source('utils.R')
load.packages()

custom_font_size <- 28

stats_type <- 'prevalence'
lusc <- load.dataset(
  'lusc',
  stats_type,
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats
lusc$typ_raka = 'LUSC'
lusc$disease_stage <- factor(lusc$disease_stage, levels = c('early', 'middle', 'late'))
lusc$smoking_status <- factor(lusc$smoking_status, levels = c('non_smoker', 'light_smoker', 'heavy_smoker'))
luad <- load.dataset(
  'luad',
  stats_type,
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats
luad$typ_raka = 'LUAD'
luad$disease_stage <- factor(luad$disease_stage, levels = c('early', 'middle', 'late'))
luad$smoking_status <- factor(luad$smoking_status, levels = c('non_smoker', 'light_smoker', 'heavy_smoker'))
data <- rbind(luad, lusc)

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ typ_raka, data = data)
ggsurvplot(
  survival.fit,
  xlab = "Lata",
  data = data,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  legend.labs = c("LUAD", "LUSC"),
  conf.int = T,
  risk.table = T,
  legend.title = "Typ raka",
  break.time.by = 365.25 * 5,
  pval.size = 8,
  pval = T,
)

mutations <-
  c('TP53', 'KRAS', 'EGFR', 'NFE2L2', 'CDKN2A', 'RET', 'STK11', 'RB1')



mutations.dataset <- luad
pdf(
  'mutations.pdf',
  width = 16,
  height = 18,
  encoding = 'ISOLatin2'
)
splots <- list()
survival.fit <-
  survfit(Surv(time = time, vital_status) ~ TP53, data = mutations.dataset)
splots[[1]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ STK11, data = mutations.dataset)
splots[[2]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ EGFR, data = mutations.dataset)
splots[[3]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )


survival.fit <-
  survfit(Surv(time = time, vital_status) ~ RET, data = mutations.dataset)
splots[[4]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )


survival.fit <-
  survfit(Surv(time = time, vital_status) ~ CDKN2A, data = mutations.dataset)
splots[[5]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )

arrange_ggsurvplots(
  splots,
  print = TRUE,
  ncol = 2,
  nrow = 3,
  risk.table.height = 0.2,
  fontsize = 1.5
)
# dev.off()
#

mutations.dataset <- lusc
# pdf(
#   'lusc_mutations.pdf',
#   width = 16,
#   height = 18,
#   encoding = 'ISOLatin2'
# )
splots <- list()
survival.fit <-
  survfit(Surv(time = time, vital_status) ~ TP53, data = mutations.dataset)
splots[[1]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ NFE2L2, data = mutations.dataset)
splots[[2]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ CDKN2A, data = mutations.dataset)
splots[[3]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )


survival.fit <-
  survfit(Surv(time = time, vital_status) ~ RB1, data = mutations.dataset)
splots[[4]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )


survival.fit <-
  survfit(Surv(time = time, vital_status) ~ KRAS, data = mutations.dataset)
splots[[5]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = mutations.dataset,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Mutacja",
    break.time.by = 365.25 * 5,
    font.custom_font_size = c(18, "plain", "black"),
    font.y = c(18, "plain", "black"),
    font.tickslab = c(18, "plain", "black"),
    font.legend = c(18, "plain", "black"),
    pval.size = 8
  )


arrange_ggsurvplots(
  splots,
  print = TRUE,
  ncol = 2,
  nrow = 3,
  risk.table.height = 0.2,
  fontsize = 1.5
)
dev.off()


####################################################################
####################################################################
####################################################################
####################################################################

font.custom_font_size.meta <- c(custom_font_size, "plain", "black")
font.y.meta = c(custom_font_size, "plain", "black")
font.tickslab.meta = c(custom_font_size, "plain", "black")
font.legend.meta = c(custom_font_size, "plain", "black")
  
  
pdf(
  'metadata_2_per_page.pdf',
  width = 18,
  height = 20,
  encoding = 'ISOLatin2'
)

splots <- list()
survival.fit <-
  survfit(Surv(time = time, vital_status) ~ age, data = luad)
splots[[1]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = luad,
    title = 'LUAD',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Grupa wiekowa",
    legend.labs = c('wiek < 65 lat', 'wiek >= 65 lat'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ age, data = lusc)
splots[[2]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = lusc,
    title = 'LUSC',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Grupa wiekowa",
    legend.labs = c('wiek < 65 lat', 'wiek >= 65 lat'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

for(i in c(1, 2)) {
  splots[[i]]$table$theme$axis.text.x$size = custom_font_size
  splots[[i]]$table$theme$axis.text.y$size = custom_font_size
  splots[[i]]$table$theme$axis.title.x$size = custom_font_size
  splots[[i]]$table$theme$axis.title.y$size = custom_font_size
}

arrange_ggsurvplots(
  splots,
  ncol = 1,
  nrow = 2,
)

splots <- list()
survival.fit <-
  survfit(Surv(time = time, vital_status) ~ sex, data = luad)
splots[[1]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = luad,
    title = 'LUAD',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Płeć",
    legend.labs = c('kobieta', 'mężczyzna'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ sex, data = lusc)
splots[[2]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = lusc,
    title = 'LUSC',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Płeć",
    legend.labs = c('kobieta', 'mężczyzna'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

for(i in c(1, 2)) {
  splots[[i]]$table$theme$axis.text.x$size = custom_font_size
  splots[[i]]$table$theme$axis.text.y$size = custom_font_size
  splots[[i]]$table$theme$axis.title.x$size = custom_font_size
  splots[[i]]$table$theme$axis.title.y$size = custom_font_size
}

arrange_ggsurvplots(
  splots,
  ncol = 1,
  nrow = 2,
)

splots <- list()
survival.fit <-
  survfit(Surv(time = time, vital_status) ~ disease_stage, data = luad)
splots[[1]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = luad,
    title = 'LUAD',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Stadium raka",
    legend.labs = c('wczesne', 'lokalnie-zaawans.', 'zaawansowne'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ disease_stage, data = lusc)
splots[[2]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = lusc,
    title = 'LUSC',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Stadium raka",
    legend.labs = c('wczesne', 'lokalnie-zaawans.', 'zaawansowane'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

for(i in c(1, 2)) {
  splots[[i]]$table$theme$axis.text.x$size = custom_font_size
  splots[[i]]$table$theme$axis.text.y$size = custom_font_size
  splots[[i]]$table$theme$axis.title.x$size = custom_font_size
  splots[[i]]$table$theme$axis.title.y$size = custom_font_size
}

arrange_ggsurvplots(
  splots,
  ncol = 1,
  nrow = 2,
)

splots <- list()
survival.fit <-
  survfit(Surv(time = time, vital_status) ~ smoking_status, data = luad)
splots[[1]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = luad,
    title = 'LUAD',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Profile palących",
    legend.labs = c('niepalący', 'sporadyczny', 'nałogowy'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ smoking_status, data = lusc)
splots[[2]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUSC',
    data = lusc,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Profile palących",
    legend.labs = c('niepalący', 'sporadyczny', 'nałogowy'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

for(i in c(1, 2)) {
  splots[[i]]$table$theme$axis.text.x$size = custom_font_size
  splots[[i]]$table$theme$axis.text.y$size = custom_font_size
  splots[[i]]$table$theme$axis.title.x$size = custom_font_size
  splots[[i]]$table$theme$axis.title.y$size = custom_font_size
}

arrange_ggsurvplots(
  splots,
  ncol = 1,
  nrow = 2,
)
dev.off()

####################################################################
####################################################################
####################################################################
####################################################################

pdf(
  'metadata_4_per_page.pdf',
  width = 22,
  height = 22,
  encoding = 'ISOLatin2'
)

splots <- list()
survival.fit <-
  survfit(Surv(time = time, vital_status) ~ age, data = luad)
splots[[1]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = luad,
    title = 'LUAD',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Grupa wiekowa",
    legend.labs = c('wiek < 65 lat', 'wiek >= 65 lat'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ age, data = lusc)
splots[[3]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = lusc,
    title = 'LUSC',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Grupa wiekowa",
    legend.labs = c('wiek < 65 lat', 'wiek >= 65 lat'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )


survival.fit <-
  survfit(Surv(time = time, vital_status) ~ sex, data = luad)
splots[[2]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = luad,
    title = 'LUAD',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Płeć",
    legend.labs = c('kobieta', 'mężczyzna'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ sex, data = lusc)
splots[[4]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = lusc,
    title = 'LUSC',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Płeć",
    legend.labs = c('kobieta', 'mężczyzna'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

for(i in 1:length(splots)) {
  splots[[i]]$table$theme$axis.text.x$size = custom_font_size
  splots[[i]]$table$theme$axis.text.y$size = custom_font_size
  splots[[i]]$table$theme$axis.title.x$size = custom_font_size
  splots[[i]]$table$theme$axis.title.y$size = custom_font_size
  splots[[i]]$table = splots[[i]]$table + theme(plot.title = element_text(size=custom_font_size))
  splots[[i]]$table$labels$y = ""
}
arrange_ggsurvplots(
  splots,
  ncol = 2,
  nrow = 2,
  risk.table.height = 0.22
)

splots <- list()
survival.fit <-
  survfit(Surv(time = time, vital_status) ~ disease_stage, data = luad)
splots[[1]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = luad,
    title = 'LUAD',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Stadium raka",
    legend.labs = c('wczesne', 'lokalnie-zaawans.', 'zaawans.'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ disease_stage, data = lusc)
splots[[3]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = lusc,
    title = 'LUSC',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Stadium raka",
    legend.labs = c('wczesne', 'lokalnie-zaawans.', 'zaawans.'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ smoking_status, data = luad)
splots[[2]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    data = luad,
    title = 'LUAD',
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Profile palących",
    legend.labs = c('niepalący', 'sporadyczny', 'nałogowy'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )

survival.fit <-
  survfit(Surv(time = time, vital_status) ~ smoking_status, data = lusc)
splots[[4]] <-
  ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUSC',
    data = lusc,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Profile palących",
    legend.labs = c('niepalący', 'sporadyczny', 'nałogowy'),
    risk.table = T,
    break.time.by = 365.25 * 5,
    font.custom_font_size = font.custom_font_size.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
  )


for(i in 1:length(splots)) {
  splots[[i]]$table$theme$axis.text.x$size = custom_font_size
  splots[[i]]$table$theme$axis.text.y$size = custom_font_size
  splots[[i]]$table$theme$axis.title.x$size = custom_font_size
  splots[[i]]$table$theme$axis.title.y$size = custom_font_size
  splots[[i]]$table = splots[[i]]$table + theme(plot.title = element_text(size=custom_font_size))
  splots[[i]]$table$labels$y = ""
}
arrange_ggsurvplots(
  splots,
  ncol = 2,
  nrow = 2,
  risk.table.height = 0.22,
  
)
dev.off()


####################################################################
####################################################################
####################################################################
####################################################################

font.custom_font_size.meta <- c(20, "plain", "black")
font.y.meta = c(20, "plain", "black")
font.tickslab.meta = c(20, "plain", "black")
font.legend.meta = c(20, "plain", "black")

pdf(
  'tissues_luad.pdf',
  width = 22,
  height = 22,
  encoding = 'ISOLatin2'
)

data <- load.dataset(
  'luad',
  'stats',
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats

for (tissue in TISSUES) {
  y <- surv_cutpoint(
    patients.mutations.stats,
    time = "time",
    event = "vital_status",
    variables = c(tissue)
  )
  # new.levels <- ifelse(data[, c(tissue)] > mean(data[, c(tissue), drop = T]), 'WYSOKI', 'NISKI')
  new.levels <- ifelse(data[, c(tissue)] > y$cutpoint[[1]], 'WYSOKI', 'NISKI')
  data[paste0(tissue, '_LEVEL')] <- as.factor(new.levels)
}

c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'MIXED'), # luad prevalence
c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA'), # lusc prevalence

c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'BRONCHI'), # luad micro
c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'LUNG') # lusc micro

splots <- list()

survival.fit <- survfit(Surv(time = time, vital_status) ~ VESSEL_LEVEL, data = data)
splots[[3]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'LUAD',
  data = data,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom tkanki naczyniowej",
  legend.labs = c('wysoki', 'niski'),
  risk.table = T,
  break.time.by = 365.25 * 5,
  font.custom_font_size = font.custom_font_size.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
)


survival.fit <- survfit(Surv(time = time, vital_status) ~ NECROSIS_LEVEL, data = data)
splots[[2]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'LUAD',
  data = data,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom martwicy",
  legend.labs = c('wysoki', 'niski'),
  risk.table = T,
  break.time.by = 365.25 * 5,
  font.custom_font_size = font.custom_font_size.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
)

survival.fit <- survfit(Surv(time = time, vital_status) ~ IMMUNE_LEVEL, data = data)
splots[[1]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'LUAD',
  data = data,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom tkanki immunologicznej",
  legend.labs = c('wysoki', 'niski'),
  risk.table = T,
  break.time.by = 365.25 * 5,
  font.custom_font_size = font.custom_font_size.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
)

data <- load.dataset(
  'lusc',
  'stats',
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats

new.levels <- ifelse(data[, c('IMMUNE')] > mean(data[, c('IMMUNE'), drop = T]), 'WYSOKI', 'NISKI')
data['IMMUNE_LEVEL'] <- as.factor(new.levels)

survival.fit <- survfit(Surv(time = time, vital_status) ~ IMMUNE_LEVEL, data = data)
splots[[4]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'LUSC',
  data = data,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom tkanki immunologicznej",
  legend.labs = c('wysoki', 'niski'),
  risk.table = T,
  break.time.by = 365.25 * 5,
  font.custom_font_size = font.custom_font_size.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
)

for(i in 1:length(splots)) {
  splots[[i]]$table$theme$axis.text.x$size = custom_font_size
  splots[[i]]$table$theme$axis.text.y$size = custom_font_size
  splots[[i]]$table$theme$axis.title.x$size = custom_font_size
  splots[[i]]$table$theme$axis.title.y$size = custom_font_size
  splots[[i]]$table = splots[[i]]$table + theme(plot.title = element_text(size=custom_font_size))
  splots[[i]]$table$labels$y = ""
}
arrange_ggsurvplots(
  splots,
  ncol = 2,
  nrow = 2,
  risk.table.height = 0.22
)
# dev.off()

splots <- list()

data <- load.dataset(
  'luad',
  'counts',
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats

for (tissue in TISSUES) {
  new.levels <- ifelse(data[, c(tissue)] > mean(data[, c(tissue), drop = T]), 'WYSOKI', 'NISKI')
  data[paste0(tissue, '_LEVEL')] <- as.factor(new.levels)
}

survival.fit <- survfit(Surv(time = time, vital_status) ~ VESSEL_LEVEL, data = data)
splots[[2]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'LUAD',
  data = data,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom komórek naczyniowych",
  legend.labs = c('wysoki', 'niski'),
  risk.table = T,
  break.time.by = 365.25 * 5,
  font.custom_font_size = font.custom_font_size.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
)

survival.fit <- survfit(Surv(time = time, vital_status) ~ IMMUNE_LEVEL, data = data)
splots[[1]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'LUAD',
  data = data,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom komórek immunologicznych",
  legend.labs = c('wysoki', 'niski'),
  risk.table = T,
  break.time.by = 365.25 * 5,
  font.custom_font_size = font.custom_font_size.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
)

data <- load.dataset(
  'lusc',
  'counts',
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats

new.levels <- ifelse(data[, c('NECROSIS')] > mean(data[, c('NECROSIS'), drop = T]), 'WYSOKI', 'NISKI')
data['NECROSIS_LEVEL'] <- as.factor(new.levels)

survival.fit <- survfit(Surv(time = time, vital_status) ~ NECROSIS_LEVEL, data = data)
splots[[3]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'LUSC',
  data = data,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom martwicy",
  legend.labs = c('wysoki', 'niski'),
  risk.table = T,
  break.time.by = 365.25 * 5,
  font.custom_font_size = font.custom_font_size.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(28, 'plain', 'bold')
)


for(i in 1:length(splots)) {
  splots[[i]]$table$theme$axis.text.x$size = custom_font_size
  splots[[i]]$table$theme$axis.text.y$size = custom_font_size
  splots[[i]]$table$theme$axis.title.x$size = custom_font_size
  splots[[i]]$table$theme$axis.title.y$size = custom_font_size
  splots[[i]]$table = splots[[i]]$table + theme(plot.title = element_text(size=custom_font_size))
  splots[[i]]$table$labels$y = ""
}
arrange_ggsurvplots(
  splots,
  ncol = 2,
  nrow = 2,
  risk.table.height = 0.22
)
dev.off()

################################################################################
################################################################################
################################################################################
pdf('hazard_plots.pdf',
  width = 13,
  height= 17,
  encoding = 'ISOLatin2'
)


source('models.R')
source('utils.R')
load.packages()

cancer_type <- 'luad'
stats_name <- 'microenvironment'
add_pathological_stage <- T
add_smoking_data <- F
windowed_data <- F
unified <- T
k <- 10
verbose <- T

stats_type <- get.stats.type(stats_name, windowed_data)

# Load desired dataset
data <- load.dataset(cancer_type, stats_type, windowed_data,
                     unified = unified,
                     add_smoking_data = add_smoking_data,
                     add_pathological_stage = add_pathological_stage,
                     mutations_variant = 0.05)
patients <- data$patients
patients.mutations <- data$patients.mutations
patients.mutations.stats <- data$patients.mutations.stats

patients.mutations.stats$disease_stage <- as.factor(patients.mutations.stats$disease_stage)
levels(patients.mutations.stats$disease_stage) <- c('wczesne', 'zaawans.', 'lokalnie-zaawans.')
patients.mutations.stats$disease_stage <- factor(patients.mutations.stats$disease_stage,
                                                  levels = c('wczesne', 'lokalnie-zaawans.', 'zaawans.'))


colnames(patients.mutations.stats) <- c('patient_id', 'time', 'vital_status', 'grupa_wiekowa', 'płeć', 'stadium_raka',
              'EGFR', 'STK11', 'TP53', 'tkanka zrębu', 'tkanka mieszana', 'tkanka immunologiczna',
              'naczynia krwionośne', 'oskrzela', 'martwica', 'tkanka płucna',
              'ITLR', 'Morisita-Horn', 'Shannon', 'Simpson')

colnames(patients.mutations) <- c('patient_id', 'time', 'vital_status', 'grupa_wiekowa', 'płeć', 'stadium raka',
                                  'EGFR', 'STK11', 'TP53')

set.seed(18)

metadata.mutations.stats.formula <- prepare.formula(c(colnames(patients.mutations),
                                                      c('tkanka zrębu', 'tkanka mieszana', 'tkanka immunologiczna',
                                                        'naczynia krwionośne', 'oskrzela', 'martwica', 'tkanka płucna')))
all.cox <- train.cox.model.counts(cox_data = patients.mutations.stats,
                                  cox_formula = metadata.mutations.stats.formula,
                                  k = k,
                                  verbose = verbose)

generated_cox_model <- coxph(formula = metadata.mutations.stats.formula, 
                             data = patients.mutations.stats)
custom_font_size <- ggforest(generated_cox_model, data = patients.mutations.stats,
              fontsize = 1.35, main='', refLabel = 'odniesienie')
custom_font_size


##################################################################################
source('models.R')
source('utils.R')
load.packages()

cancer_type <- 'luad'
stats_name <- 'prevalence'
add_pathological_stage <- T
add_smoking_data <- F
windowed_data <- F
unified <- T
k <- 10
verbose <- T

stats_type <- get.stats.type(stats_name, windowed_data)

# Load desired dataset
data <- load.dataset(cancer_type, stats_type, windowed_data,
                     unified = unified,
                     add_smoking_data = add_smoking_data,
                     add_pathological_stage = add_pathological_stage,
                     mutations_variant = 0.1)
patients <- data$patients
patients.mutations <- data$patients.mutations
patients.mutations.stats <- data$patients.mutations.stats

patients.mutations.stats$disease_stage <- as.factor(patients.mutations.stats$disease_stage)
levels(patients.mutations.stats$disease_stage) <- c('wczesne', 'zaawans.', 'lokalnie-zaawans.')
patients.mutations.stats$disease_stage <- factor(patients.mutations.stats$disease_stage,
                                                 levels = c('wczesne', 'lokalnie-zaawans.', 'zaawans.'))

# patients.mutations.stats$smoking_status <- as.factor(patients.mutations.stats$smoking_status)
# levels(patients.mutations.stats$smoking_status) <- c('nałogowy', 'sporadyczny', 'niepalący')
# patients.mutations.stats$smoking_status <- factor(patients.mutations.stats$smoking_status,
#                                                  levels = c('niepalący', 'sporadyczny', 'nałogowy'))
# 
# colnames(patients.mutations.stats) <- c('patient_id', 'time', 'vital_status', 'grupa_wiekowa', 'płeć', 
#                                         'profile palaczy', 'stadium_raka',
#                                         'CDKN2A', 'EGFR', 'STK11', 'TP53', 'RET', 'tkanka nowotworowa', 
#                                         'tkanka zrębu', 'tkanka mieszana', 'tkanka immunologiczna',
#                                         'naczynia krwionośne', 'oskrzela', 'martwica', 'tkanka płucna',
#                                         'ITLR', 'Morisita-Horn', 'Shannon', 'Simpson', "współczynnik limfocytów (LR)")
# 
# colnames(patients.mutations) <- c('patient_id', 'time', 'vital_status', 'grupa_wiekowa', 'płeć', 'profile palaczy',
#                                   'stadium_raka', 'CDKN2A', 'EGFR', 'STK11', 'TP53', 'RET')


colnames(patients.mutations.stats) <- c('patient_id', 'time', 'vital_status', 'grupa_wiekowa', 'płeć',
                                        'stadium raka', 'CDKN2A', 'EGFR', 'STK11', 'TP53', 'RET', 
                                        'tkanka nowotworowa',  'tkanka zrębu', 'tkanka mieszana', 
                                        'tkanka immunologiczna', 'naczynia krwionośne', 'oskrzela', 
                                        'martwica', 'tkanka płucna', 'ITLR', 'Morisita-Horn', 
                                        'Shannon', 'Simpson', "współczynnik limfocytów (LR)")

colnames(patients.mutations) <- c('patient_id', 'time', 'vital_status', 'grupa_wiekowa', 'płeć',
                                  'stadium_raka', 'CDKN2A', 'EGFR', 'STK11', 'TP53', 'RET')

set.seed(18)

metadata.mutations.stats.formula <- prepare.formula(
  c(colnames(patients.mutations),
    c('tkanka nowotworowa', 'tkanka zrębu', 'tkanka mieszana', 'tkanka immunologiczna',
      'naczynia krwionośne', 'oskrzela', 'martwica', 'tkanka płucna')))
all.cox <- train.cox.model.counts(cox_data = patients.mutations.stats,
                                  cox_formula = metadata.mutations.stats.formula,
                                  k = k,
                                  verbose = verbose)

generated_cox_model <- coxph(formula = metadata.mutations.stats.formula, 
                             data = patients.mutations.stats)
custom_font_size <- ggforest(generated_cox_model, data = patients.mutations.stats,
              fontsize = 1.35, main='', refLabel='odniesienie')
custom_font_size

dev.off()


source('models.R')
source('constants.R')
source('utils.R')
load.packages()

x <- 28

c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'MIXED') # luad prevalence
c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA') # lusc prevalence

c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'BRONCHI') # luad micro
c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA', 'LUNG') # lusc micro

font.x.meta <- c(28, "plain", "black")
font.y.meta = c(28, "plain", "black")
font.tickslab.meta = c(28, "plain", "black")
font.legend.meta = c(28, "plain", "black")

pdf(
  'tissues.pdf',
  width = 22,
  height = 22,
  encoding = 'ISOLatin2'
)

for(i in c('stats', 'counts')) {
  tissues <- c('VESSEL', 'IMMUNE', 'NECROSIS', 'STROMA')
  z <- 4
  if(i == 'stats') {
    tissues <- c('TUMOR', tissues)
    z <- 5
  }

  luad <- load.dataset(
    'luad',
    i,
    windowed_data = F,
    unified = T,
    add_smoking_data = T,
    add_pathological_stage = T
  )$patients.mutations.stats


  for (tissue in tissues) {
    y <- surv_cutpoint(
      luad,
      time = "time",
      event = "vital_status",
      variables = c(tissue)
    )
    new.levels <- ifelse(luad[, c(tissue)] > y$cutpoint[[1]], 'wysoki', 'niski')
    luad[paste0(tissue, '_LEVEL')] <- as.factor(new.levels)
  }

  lusc <- load.dataset(
    'lusc',
    i,
    windowed_data = F,
    unified = T,
    add_smoking_data = T,
    add_pathological_stage = T
  )$patients.mutations.stats


  for (tissue in tissues) {
    y <- surv_cutpoint(
      lusc,
      time = "time",
      event = "vital_status",
      variables = c(tissue)
    )
    new.levels <- ifelse(lusc[, c(tissue)] > y$cutpoint[[1]], 'wysoki', 'niski')
    lusc[paste0(tissue, '_LEVEL')] <- as.factor(new.levels)
  }


  splots <- list()


  survival.fit <- survfit(Surv(time = time, vital_status) ~ VESSEL_LEVEL, data = luad)
  splots[[1]] <- ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUAD',
    data = luad,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Poziom naczyń krwionośnych",
    legend.labs = c('niski', 'wysoki'),
    # risk.table = T,
    break.time.by = 365.25 * 5,
    font.x = font.x.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
  )

  survival.fit <- survfit(Surv(time = time, vital_status) ~ VESSEL_LEVEL, data = lusc)
  splots[[z+1]] <- ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUSC',
    data = lusc,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Poziom naczyń krwionośnych",
    legend.labs = c('niski', 'wysoki'),
    # risk.table = T,
    break.time.by = 365.25 * 5,
    font.x = font.x.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
  )

  survival.fit <- survfit(Surv(time = time, vital_status) ~ IMMUNE_LEVEL, data = luad)
  splots[[2]] <- ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUAD',
    data = luad,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Poziom tkanki immunologicznej",
    legend.labs = c('niski', 'wysoki'),
    # risk.table = T,
    break.time.by = 365.25 * 5,
    font.x = font.x.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
  )

  survival.fit <- survfit(Surv(time = time, vital_status) ~ IMMUNE_LEVEL, data = lusc)
  splots[[z+2]] <- ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUSC',
    data = lusc,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Poziom tkanki immunologicznej",
    legend.labs = c('niski', 'wysoki'),
    # risk.table = T,
    break.time.by = 365.25 * 5,
    font.x = font.x.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
  )

  survival.fit <- survfit(Surv(time = time, vital_status) ~ NECROSIS_LEVEL, data = luad)
  splots[[3]] <- ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUAD',
    data = luad,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Poziom martwicy",
    legend.labs = c('niski', 'wysoki'),
    # risk.table = T,
    break.time.by = 365.25 * 5,
    font.x = font.x.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
  )

  survival.fit <- survfit(Surv(time = time, vital_status) ~ NECROSIS_LEVEL, data = lusc)
  splots[[z+3]] <- ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUSC',
    data = lusc,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Poziom martwicy",
    legend.labs = c('niski', 'wysoki'),
    # risk.table = T,
    break.time.by = 365.25 * 5,
    font.x = font.x.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
  )

  survival.fit <- survfit(Surv(time = time, vital_status) ~ STROMA_LEVEL, data = luad)
  splots[[4]] <- ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUAD',
    data = luad,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Poziom tkanki zrębu",
    legend.labs = c('niski', 'wysoki'),
    # risk.table = T,
    break.time.by = 365.25 * 5,
    font.x = font.x.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
  )


  survival.fit <- survfit(Surv(time = time, vital_status) ~ STROMA_LEVEL, data = lusc)
  splots[[z+4]] <- ggsurvplot(
    survival.fit,
    xlab = "Lata",
    title = 'LUSC',
    data = lusc,
    xscale = "d_y",
    risk.table.title = "Liczba zagrożonych pacjentów",
    ylab = "Szanse przeżycia",
    mark.time = TRUE,
    pval = T,
    conf.int = T,
    legend.title = "Poziom tkanki zrębu",
    legend.labs = c('niski', 'wysoki'),
    # risk.table = T,
    break.time.by = 365.25 * 5,
    font.x = font.x.meta,
    font.y = font.y.meta,
    font.tickslab = font.tickslab.meta,
    font.legend = font.legend.meta,
    pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
  )

  if(i == 'stats') {
    survival.fit <- survfit(Surv(time = time, vital_status) ~ TUMOR_LEVEL, data = luad)
    splots[[5]] <- ggsurvplot(
      survival.fit,
      xlab = "Lata",
      title = 'LUAD',
      data = luad,
      xscale = "d_y",
      risk.table.title = "Liczba zagrożonych pacjentów",
      ylab = "Szanse przeżycia",
      mark.time = TRUE,
      pval = T,
      conf.int = T,
      legend.title = "Poziom tkanki nowotorowrej",
      legend.labs = c('niski', 'wysoki'),
      # risk.table = T,
      break.time.by = 365.25 * 5,
      font.x = font.x.meta,
      font.y = font.y.meta,
      font.tickslab = font.tickslab.meta,
      font.legend = font.legend.meta,
      pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
    )

    survival.fit <- survfit(Surv(time = time, vital_status) ~ TUMOR_LEVEL, data = lusc)
    splots[[z+5]] <- ggsurvplot(
      survival.fit,
      xlab = "Lata",
      title = 'LUSC',
      data = lusc,
      xscale = "d_y",
      risk.table.title = "Liczba zagrożonych pacjentów",
      ylab = "Szanse przeżycia",
      mark.time = TRUE,
      pval = T,
      conf.int = T,
      legend.title = "Poziom tkanki nowotworowej",
      legend.labs = c('niski', 'wysoki'),
      # risk.table = T,
      break.time.by = 365.25 * 5,
      font.x = font.x.meta,
      font.y = font.y.meta,
      font.tickslab = font.tickslab.meta,
      font.legend = font.legend.meta,
      pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
    )
  }


  for(i in 1:length(splots)) {
    splots[[i]]$table$theme$axis.text.x$size = x
    splots[[i]]$table$theme$axis.text.y$size = x
    splots[[i]]$table$theme$axis.title.x$size = x
    splots[[i]]$table$theme$axis.title.y$size = x
    splots[[i]]$table = splots[[i]]$table + theme(plot.title = element_text(size=x))
    splots[[i]]$table$labels$y = ""
  }

  arrange_ggsurvplots(
    splots,
    ncol = 2,
    nrow = z
  )
}

luad <- load.dataset(
  'luad',
  'counts',
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats


for (tissue in c('BRONCHI')) {
  y <- surv_cutpoint(
    luad,
    time = "time",
    event = "vital_status",
    variables = c(tissue)
  )
  new.levels <- ifelse(luad[, c(tissue)] > y$cutpoint[[1]], 'WYSOKI', 'NISKI')
  luad[paste0(tissue, '_LEVEL')] <- as.factor(new.levels)
}

lusc <- load.dataset(
  'lusc',
  'counts',
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats


for (tissue in c('LUNG')) {
  y <- surv_cutpoint(
    lusc,
    time = "time",
    event = "vital_status",
    variables = c(tissue)
  )
  new.levels <- ifelse(lusc[, c(tissue)] > y$cutpoint[[1]], 'WYSOKI', 'NISKI')
  lusc[paste0(tissue, '_LEVEL')] <- as.factor(new.levels)
}


splots <- list()

survival.fit <- survfit(Surv(time = time, vital_status) ~ BRONCHI_LEVEL, data = luad)
splots[[1]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'Tkanka oskrzelowa w mikrośrodowisku nowotoworu dla gruczolakoraka',
  data = luad,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom tkanki oskrzelowej",
  legend.labs = c('niski', 'wysoki'),
  # risk.table = T,
  break.time.by = 365.25 * 5,
  font.x = font.x.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
)

survival.fit <- survfit(Surv(time = time, vital_status) ~ LUNG_LEVEL, data = lusc)
splots[[2]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'Tkanka płucna w mikrośrodowisku nowotoworu dla raka płaskonabłonkowego (LUSC)',
  data = lusc,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom tkanki płucnej",
  legend.labs = c('niski', 'wysoki'),
  # risk.table = T,
  break.time.by = 365.25 * 5,
  font.x = font.x.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
)

luad <- load.dataset(
  'luad',
  'stats',
  windowed_data = F,
  unified = T,
  add_smoking_data = T,
  add_pathological_stage = T
)$patients.mutations.stats


for (tissue in c('MIXED')) {
  y <- surv_cutpoint(
    luad,
    time = "time",
    event = "vital_status",
    variables = c(tissue)
  )
  new.levels <- ifelse(luad[, c(tissue)] > y$cutpoint[[1]], 'WYSOKI', 'NISKI')
  luad[paste0(tissue, '_LEVEL')] <- as.factor(new.levels)
}

survival.fit <- survfit(Surv(time = time, vital_status) ~ MIXED_LEVEL, data = luad)
splots[[3]] <- ggsurvplot(
  survival.fit,
  xlab = "Lata",
  title = 'Zliczenia tkanki mieszanej dla gruczokaraka (LUAD)',
  data = luad,
  xscale = "d_y",
  risk.table.title = "Liczba zagrożonych pacjentów",
  ylab = "Szanse przeżycia",
  mark.time = TRUE,
  pval = T,
  conf.int = T,
  legend.title = "Poziom tkanki mieszanej",
  legend.labs = c('niski', 'wysoki'),
  # risk.table = T,
  break.time.by = 365.25 * 5,
  font.x = font.x.meta,
  font.y = font.y.meta,
  font.tickslab = font.tickslab.meta,
  font.legend = font.legend.meta,
  pval.size = 10, fontsize = 10, font.title=c(30, 'plain', 'bold')
)


for(i in 1:length(splots)) {
  splots[[i]]$table$theme$axis.text.x$size = x
  splots[[i]]$table$theme$axis.text.y$size = x
  splots[[i]]$table$theme$axis.title.x$size = x
  splots[[i]]$table$theme$axis.title.y$size = x
  splots[[i]]$table = splots[[i]]$table + theme(plot.title = element_text(size=x))
  splots[[i]]$table$labels$y = ""
}
arrange_ggsurvplots(
  splots,
  ncol = 1,
  nrow = 3,
)
dev.off()
