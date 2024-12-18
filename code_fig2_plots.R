### General information ----
# Title: Establishing a random forest model to predict the maturity of adult, human plasma cells using various cell surface marker
# Author: Tobit D. Steinmetz
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: d.t.steinmetz@umcg.nl
# Collaboration: please ask permission from the author before using this script
# Date created: 19-10-2023
# Date last adjustment: 03-12-2024
# RStudio version: 2023.06.1 Build 524
# References: 

### introduction ----
# this scripts analysis imputed and unimputed data used for training the ASC-ME model
# plots for individual maturity marker are generated for immune stage groups and age correlations
# 

### This R-script requires the following packages:
library(dplyr)
library(ggplot2)
library(readxl)

## contine with R environment from training the RF model
load(".../ASC_ME_prediction_model/env_train_model_final.RData")
setwd(".../ASC_ME_prediction_model")

###dpi grouping for imputed data
imp_markers<-data.imputed_mean %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 4.9, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "excl","early", "intermediate", "late", "very late"))) # Define the labels for the age groups
imp_markers<-cbind(imp_markers,StudyID=dataforimp$StudyID)
#
for (a in 4:11) {
ANOVA_plot<- lm(imp_markers[[a]] ~ imp_markers$dpi_group)  |> aov() |> TukeyHSD()
plot_imp_marker<-ggplot(imp_markers, aes(x = dpi_group, y = imp_markers[,a])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(imp_markers$StudyID)), width = 0.25, alpha = 1, size = 2,shape=1) +
  theme_bw() +
  labs(title = "one-way ANOVA with Tukey correction",
       x = "immune response stage",
       y = paste(colnames(imp_markers[a]), "(normalized mean expression)"),
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  xlim("baseline", "early", "intermediate", "late", "very late") +
  annotate("text", x="early", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][2,4], format = "e", digits=2)), size=5) +
  annotate("text", x="intermediate", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][3,4], format = "e", digits=2)), size=5) +
  annotate("text", x="late", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][4,4], format = "e", digits=2)), size=5) +
  annotate("text", x="very late", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][5,4], format = "e", digits=2)), size=5)
print(plot_imp_marker) 
ggsave(filename = paste("plot_imp_",colnames(imp_markers[a]),".png"), plot = plot_imp_marker)
}
## output for Figure 2

###age correlation for imputed data
imp_markers$Age<-as.numeric(imp_markers$Age)
imp_markers_BL<-imp_markers %>% filter(dpi_group=="baseline")
for (a in 4:11) {
  r_plot<-cor.test(imp_markers_BL$Age, imp_markers_BL[[a]])
  print(ggplot(imp_markers_BL, aes(x=Age, y=imp_markers_BL[,a])) +
    geom_point(shape=1) + #none-filled circles
    geom_smooth(method = "lm", color = "red", se = TRUE ) + #add regression line /w confidence interval
    labs(title = paste("pearson correlation, r =",format(r_plot$estimate, digits=3)," p =",format(r_plot$p.value, digits=3)),
         x = "sample age",
         y = paste(colnames(imp_markers_BL[a]), "(normalized mean expression)")) +
    theme_bw() +
    ylim(-1,5)+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.y = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.x = element_text(size = 20,face = "bold"),
          plot.title = element_text(size = 15,face = "bold"))
  )
  ggsave(filename = paste("plot_BL_",colnames(imp_markers_BL[a]),"_age.png", sep = ""))
}
## output for Extended Data Figure 2

###
###dpi grouping for unimputed data
unimp_markers <- read_excel("all_extended_data.xlsx") %>% 
  select(-c(No)) %>% 
  filter(!(sample=="NA")) %>%
  filter(((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT"))) %>%
  filter((StudyID=="AID")|(StudyID=="CER")|(StudyID=="SDY144")|(StudyID=="SDY180")|(StudyID=="SDY224")
         |(StudyID=="SDY272")|(StudyID=="SDY296")|(StudyID=="SDY301")|(StudyID=="SDY364")|(StudyID=="SDY387")
         |(StudyID=="SDY522")|(StudyID=="SDY80")|(StudyID=="SDY819")|(StudyID=="SDY984")|(StudyID=="SDY1086")
         |(StudyID=="SDY1288")|(StudyID=="SDY1669")|(StudyID=="SDY1734")|(StudyID=="ZY77")
  )
unimp_markers<-unimp_markers %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "early", "intermediate", "late", "very late"))) # Define the labels for the age groups
#
for (b in c(3:10,13:16,18,19)) {
  ANOVA_plot<- lm(unimp_markers[[b]] ~ unimp_markers[[30]])  |> aov() |> TukeyHSD()
  annotations <- list()
  groups <- c("early", "intermediate", "late", "very late")
  for (group in groups) {
    if (sum(!is.na(unimp_markers[unimp_markers$dpi_group == group, b])) > 0) {
      annotations <- append(annotations,annotate("text", x = group, y = max(unimp_markers[, b], na.rm = TRUE) * 1.1,
                                                 label = paste("p=", formatC(ANOVA_plot[["unimp_markers[[30]]"]][paste0(group, "-baseline"), 4],format = "e",digits = 2)
                                                 ),size = 5))}}
  plot_unimp_marker<-ggplot(unimp_markers, aes(x = dpi_group, y = unimp_markers[[b]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = as.factor(unimp_markers[[20]])), width = 0.25, alpha = 1, size = 2,shape=1) +
    theme_bw() +
    labs(title = "one-way ANOVA with Tukey correction",
         x = "immune response stage",
         y = paste(colnames(unimp_markers[b]), "(normalized mean expression)"),
         color = "StudyID") +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.y = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.x = element_text(size = 20,face = "bold"),
          plot.title = element_text(size = 15,face = "bold")) +
    ylim(min(unimp_markers[,b],na.rm = TRUE),max(unimp_markers[,b],na.rm = TRUE)*1.1) +
    xlim("baseline", "early", "intermediate", "late", "very late") +
    annotations
  print(plot_unimp_marker) 
ggsave(filename = paste("plot_unimp_",colnames(unimp_markers[b]),".png", sep = ""), plot = plot_unimp_marker)   
}
## output for Supplemental Figure 1A

#age correlation for unimputed data
unimp_markers$Age<-as.numeric(unimp_markers$Age)
unimp_markers_BL<-unimp_markers %>% filter(dpi_group=="baseline")
for (b in c(3:10,13:16,18,19)) {
  r_plot<-cor.test(unimp_markers_BL$Age, unimp_markers_BL[[b]])
  print(ggplot(unimp_markers_BL, aes(x=Age, y=unimp_markers_BL[[b]])) +
          geom_point(shape=1) + #none-filled circles
          geom_smooth(method = "lm", color = "red", se = TRUE ) + #add regression line /w confidence interval
          labs(title = paste("pearson correlation, r =",format(r_plot$estimate, digits=3)," p =",format(r_plot$p.value, digits=3)),
               x = "sample age",
               y = paste(colnames(unimp_markers_BL[b]), "(normalized mean expression)")) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.y = element_text(size = 20,face = "bold"),
                axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.x = element_text(size = 20,face = "bold"),
                plot.title = element_text(size = 15,face = "bold"))
  )
  ggsave(filename = paste("plot_BL.unimp_",colnames(unimp_markers_BL[b]),"_age.png", sep = ""))
}
## output for Supplemental Figure 1B

## End of project