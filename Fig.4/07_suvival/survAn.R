rm(list = ls())

library(survival)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(survminer) 
library(survival)


#~~~ BT + G +++++++++
data <- data.frame(
  PatientID = c("G523", "G549", "G564", "G566", "G583", "G620", "G851", "G876", "BT50", "BT67", "BT69", "BT89", "BT94", "BT238", "BT248", "BT301"),
  time = c(87.5, 122, 320.5, 317, 198, 128.5, 131, 154, 376, 154, 335, 130, 375, 317,  46, 171),
  status = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  group = c("ecDNA+", "ecDNA+", "ecDNA-", "ecDNA-", "ecDNA+", "ecDNA+", "ecDNA-", "ecDNA-", "ecDNA-", "ecDNA-", "ecDNA-", "ecDNA-", "ecDNA-", "ecDNA-", "ecDNA+", "ecDNA+")
)

xenograft <- read_excel("./inputs/xenograft.xlsx") %>% 
  as.data.frame() %>% 
  filter(!is.na(suppressWarnings(as.numeric(as.character(OrthotopicXenoSurvival_Days))))) %>%
  filter(PatientID %in% data$PatientID) %>%
  left_join(data, by = "PatientID") %>%
  select(PatientID, OrthotopicXenoSurvival_Days, Sphere_Formation_Percent, time, status, group)



#~~~ fit the Kaplan-Meier survival model
fit <- survfit(Surv(time, status) ~ group, data = data)

#~~~ plot Kaplan-Meier curve
ggsurvplot(
  fit,
  data = data,
  risk.table = F,            # Show risk table below the plot
  pval = TRUE,                  # Show p-value on the plot
  conf.int = F,              # Show confidence intervals for survival curves
  xlab = "Time to event (days)",
  ylab = "Overall survival",
  legend.labs = c(paste0("ecDNA- ", "(n=", length(which(data$group == "ecDNA-")) ,")"), paste0("ecDNA+ ", "(n=", length(which(data$group == "ecDNA+")) ,")")), # Names for the groups
  size = .6,
  palette = c("#377EB8", "red")   # Colors for the lines
)



#--- G only +++++++++
df <- data.frame(celllines = c("G523", "G549", "G564", "G566", "G583", "G620", "G637", "G851", "G876"),
                 xengfdays = c(87.5, 122, 320, 317, 198, 128, 359, 154, 131),
                 status = c(1,1,1,1,1,1,1,1,1),
                 celltypes = c("ecDNA(+)", "ecDNA(+)", "ecDNA(-)", "ecDNA(-)", "ecDNA(+)", "ecDNA(+)", "ecDNA(-)", "ecDNA(-)", "ecDNA(-)"))

survObj <- Surv(time = df$xengfdays,
                event = df$status)

fit <- survfit(survObj ~ df$celltypes, data = df)

ggsurvplot(fit,
           data = df,
           pval = TRUE,
           legend.labs = c("ecDNA(-)", "ecDNA(+)"),
           palette = c("black", "red"),
           risk.table = T,
           xlab = "Time (days)"
)
