# Analyzing multiple single CD8 staining images of CTCL 
# to conduct survival analysis
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

cat("Survival analysis of CD8 single-staining data, pre-analyzed with QuPath \n")

# Loading packages --------------------------------------------------------
libraries = c("argparse",
              "ggplot2",
              "plyr",
              "dplyr",
              "tidyr",
              "Seurat",
              "SpatialExperiment")
suppressPackageStartupMessages(invisible(lapply(libraries, require, character.only = TRUE)))
setwd("/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/")

source("scripts/MultiplexImaging_utils.R")
source("scripts/GlobalVars.R")

library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(knitr)
library(survival)
library(tibble)

output_dir = file.path("output/CTCL_response/Single_Staining")
dir.create(output_dir)

# Load Data
Survival = readxl::read_xlsx(file.path("../../ProClipi/Patients_HE_Flow.xlsx"))
Survival = Survival[which(Survival$`well definied infiltrate [1/0]` == 1),] 
Survival$`% cT` = as.numeric(Survival$`% cT`)
Survival = Survival %>% mutate(
  status = dplyr::recode(`Dead [0]/Alive [1]`, `1` = 0, `0` = 1, `Lost to follow up` = 0, `Check Follow up` = 0)
)
Survival$`Survival after biopsy [m]` = as.numeric(Survival$`Survival after biopsy [m]`)

Survival = Survival %>% mutate(
  cytotoxic_status = ifelse(`% cT` > median(`% cT`),  "High", "Low")
)

pdf(file.path(output_dir, "Distribution_CD8_in_cohort.pdf"), height = 1.5, width = 5)
Survival %>% ggplot(aes(y = `% cT`, fill = cytotoxic_status)) + 
  geom_histogram(bins = 45) + ggplot2::coord_flip() +
  theme_minimal() + 
  scale_fill_manual(values = alpha(c("High" = "#74C276", "Low" = "#B54A48"), 0.8)) +
  geom_hline(yintercept = median(Survival$`% cT`), lty = 3) +
  ylab("T cytotoxic cells (% of CD3+ cells)") + xlab("")
dev.off()
 

Survival = Survival %>% mutate(
  condition = gsub("[0-9]*","",Survival$Patient_code)
)

s1 <- survfit(Surv(`Survival after biopsy [m]`, status) ~ condition, data = Survival)
str(s1)


pdf(file.path(output_dir, "Survival_Curve_CD8.pdf"), height = 4.5, width = 5)
survfit2(Surv(`Survival after biopsy [m]`, status) ~ cytotoxic_status, data = Survival ) %>% 
  ggsurvfit(theme = theme_minimal(), size = 1) +
  scale_color_manual(values = c("High" = "#74C276", "Low" = "#B54A48")) + 
  scale_fill_manual(values = c("High" = "#74C276", "Low" = "#B54A48")) + 
  labs(
    x = "Months",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable() + ggsurvfit::add_pvalue()
dev.off()

s = survdiff(Surv(`Survival after biopsy [m]`, status) ~ cytotoxic_status, data = Survival %>%
           filter(`well definied infiltrate [1/0]` == 1))
