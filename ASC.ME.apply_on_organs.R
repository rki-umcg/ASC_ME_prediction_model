### General information ----
# Title: Establishing a random forest model to predict the maturity of adult, human plasma cells using various cell surface marker
# Author: Tobit D. Steinmetz
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: d.t.steinmetz@umcg.nl
# Collaboration: please ask permission from the author before using this script
# Date created: 19-10-2023
# Date last adjustment: 02-12-2024
# RStudio version: 2023.06.1 Build 524
# References: 

### introduction ----
# this scripts contains the application of the ASC-ME model to bone marrow and gut tissue samples
# 

### This R-script requires the following packages:
library(randomForest)
library(mice)
library(readxl)
library(dplyr)
library(tidyverse)
library(dunn.test)

## contine with R environment from training the RF model
load(".../ASC_ME_prediction_model/env_train_model_final.RData")
setwd(".../ASC_ME_prediction_model")

#set seed for imputation reproducibility and import organ-specific data
set.seed(555)
data_new <- read_excel("full_organ_data.xlsx") %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,No,CD86)) %>% 
  filter(!(sample=="NA"))

#use imputation to complete missing validation data
#the whole imputed training data is added for this purpose to cope with new data that miss markers entirely
new_merging <- data_new %>% 
  select(sample,check,SubjectID,File,condition)
new_imp <- data_new %>% 
  select(-c(check,SubjectID,File,condition))
percent_miss<-(sum(is.na(new_imp)) / (nrow(new_imp) * ncol(new_imp))) *100
#
iter=100
data.norm<-cbind(data.imputed_mean[,c(1,2,4:11)],dataforimp[,c(11:16)])
data.total<-rbind(data.norm,new_imp)
imp.new <- mice(data.total, method = "rf", m = iter, seed = 555)

#get all imputed data (as a list)
new.imputed <- list()
for(m in 1:length(imp.new)){
  new.imputed[[m]] <- complete(imp.new, action = m)
} 
for(m in 1:length(imp.new)){  
  new.imputed[[m]] <- cbind(new.imputed[[m]][c((nrow(data.norm)+1):nrow(new.imputed[[m]])),], new_merging)
}

#get average imputed values for prediction markers
new.imputed_mean <- new.imputed[[1]][,c(1,2,14)]
for (o in 1:length(prediction_marker)) {
  new_average_marker <- matrix(0, nrow(new.imputed[[1]]))
  for (m in 1:length(new.imputed)) {
    new_average_marker <- new_average_marker + new.imputed[[m]][,o+2]
  }
  new_average_marker <- new_average_marker / length(new.imputed)
  new.imputed_mean <- cbind(new.imputed_mean,new_average_marker)
  colnames(new.imputed_mean)[o+3] <- paste(prediction_marker[o])
}

###
###use model on all imputed organ-specific data and average out predictions
prediction_new <- list()
for (m in 1:length(new.imputed)) {
  prediction_new[[m]] <- average_predictions(new.imputed[[m]], rf_ascmi)
}
new_predictions <- matrix(0, nrow(prediction_new[[1]]), ncol(prediction_new[[1]]))
for (m in 1:length(prediction_new)) {
  new_predictions <- new_predictions + prediction_new[[m]][,3]
}
new_predictions <- new_predictions / length(prediction_new)
predictions_new <- cbind(new_merging,new_imp[,c(11:12,15)], prediction_new[[1]][, 1:2], new_predictions)

# trim output table
new_output<-predictions_new %>%
  select(sample, StudyID,SubjectID, organ, group, dpi, p) %>%
  mutate(Age=new_imp$Age) %>%
  mutate(Gender=new_imp$Gender)
new_output$Age<-as.numeric(new_output$Age)

### generate plots to explore the ASC-ME model application on bone marrow and gut data

### maturity prediction by organ
data_organs<-new_output %>%
  filter((organ == "blood")|(organ == "BM")|(organ == "gut"))
#merge with imported blood samples from the training data
setup_BL<-read_excel("train_data_prediction.xlsx") %>% 
  filter(dpi_group=="BL") 
setup_BL2<-setup_BL[,c(1,13,16,15,17,2,12,3)]
  colnames(setup_BL2)[7]<-"p"
data_organ_plot<-rbind(data_organs[,-9], setup_BL2) %>%
  filter(p<470)
ANOVA_imp<- lm(data_organ_plot[[7]] ~ data_organ_plot[[4]])  |> aov() |> TukeyHSD()
plot_HC_organ<-ggplot(data_organ_plot, aes(x = organ, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(data_organ_plot[[2]])), width = 0.25, alpha = 1, size = 2) +
  theme_bw() +
  labs(title = "one-way ANOVA with Tukey correction",
       x = "",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  ylim(4,96)+
  annotate("text", x="BM", y=95, label=paste("p=",formatC(ANOVA_imp[["data_organ_plot[[4]]"]][1,4], format = "e", digits=2)), size=5)+
  annotate("text", x="gut", y=95, label=paste("p=",formatC(ANOVA_imp[["data_organ_plot[[4]]"]][2,4], format = "e", digits=2)), size=5)
print(plot_HC_organ) 
ggsave(filename = "plot_all_organ.png", plot = plot_HC_organ)
## output for Figure 6A

### maturity prediction along different gut locations
Festen_data<-data_organs %>%
  filter(StudyID == "Festen") %>%
  mutate(organ = SubjectID) %>%
  mutate(across(SubjectID, ~ gsub("Festen_01", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_02", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_03", "rectum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_04", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_05", "rectum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_06", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_07", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_08", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_09", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_10", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_11", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_12", "des.col", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_13", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_14", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_15", "rectum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_16", "rectum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_17", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_18", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_19", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_20", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_21", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_22", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_23", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_24", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_25", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_26", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_27", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_28", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_29", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_30", "rectum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_31", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_32", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_33", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_34", "rectum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_35", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_36", "rectum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_37", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_38", "des.col", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_39", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_40", "ileum", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_41", "sigmoid", .))) %>%
  mutate(across(SubjectID, ~ gsub("Festen_42", "sigmoid", .)))
kt<- kruskal.test(Festen_data[[7]] ~ Festen_data[[3]])
dt <- dunn.test(Festen_data[[7]], Festen_data[[3]], method = "bonferroni", list = TRUE)
plot_Festen_organ<-ggplot(Festen_data, aes(x = SubjectID, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 1, size = 2) +
  theme_bw() +
  labs(title = "Kruskal-Wallis test with Bonferroni correction",
       x = "",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  xlim("ileum", "des.col", "sigmoid", "rectum")+
  annotate("text", x="des.col", y=80, label=paste("p=",format(dt$P.adjusted[1], digits=3)), size=5)+
  annotate("text", x="sigmoid", y=80, label=paste("p=",format(dt$P.adjusted[5], digits=3)), size=5)+
  annotate("text", x="rectum", y=80, label=paste("p=",format(dt$P.adjusted[3], digits=3)), size=5)
print(plot_Festen_organ) 
ggsave(filename = "plot_Festen_organ.png", plot = plot_Festen_organ)
## output for Figure 6C

### maturity prediction in gut samples before (T0) and after (T4) vedolizumab treatment
Festen_data2<-data_organs %>%
       filter(StudyID == "Festen") %>%
       mutate(organ = SubjectID) %>%
       mutate(across(SubjectID, ~ gsub("Festen_01", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_02", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_03", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_04", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_05", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_06", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_07", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_08", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_09", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_10", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_11", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_12", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_13", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_14", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_15", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_16", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_17", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_18", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_19", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_20", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_21", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_22", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_23", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_24", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_25", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_26", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_27", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_28", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_29", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_30", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_31", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_32", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_33", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_34", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_35", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_36", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_37", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_38", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_39", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_40", "T0", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_41", "T4", .))) %>%
       mutate(across(SubjectID, ~ gsub("Festen_42", "T4", .)))
kt2<- kruskal.test(Festen_data2[[7]] ~ Festen_data2[[3]])
dt2<- dunn.test(Festen_data2[[7]], Festen_data2[[3]], method = "bonferroni", list = TRUE)
plot_Festen_organ<-ggplot(Festen_data2, aes(x = SubjectID, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 1, size = 2) +
  theme_bw() +
  labs(title = "Mann Whitney U-test",
       x = "",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  xlim("T0", "T4")+
  annotate("text", x=1.5, y=80, label=paste("p=0.0244"), size=5)
print(plot_Festen_organ) 
ggsave(filename = "plot_Festen2_organ.png", plot = plot_Festen_organ)
## output for Figure 6C

### maturity prediction for antigen-specific ASC data from study Z7A5
g.spec <- read_excel("Ag.spec.xlsx", range = "A1:K41") %>% filter(!(SubjectID=="KM557"))
p.spec <- Ag.spec$p
pval<-wilcox.test(p.spec[c(1:19)], p.spec[c(20:38)], paired = TRUE, exact = TRUE)
set.seed(123)  # to fix the jitter pattern
Ag.spec <- Ag.spec %>%
  mutate(jittered_x = jitter(as.numeric(factor(Ag)), amount = 0.25))
print(
  ggplot(Ag.spec, aes(x = Ag, y = p)) +  # Use Ag for box plot grouping
    geom_boxplot(outlier.shape = NA) +  # Box plot uses original Ag levels
    geom_line(aes(x = jittered_x, group = SubjectID), color = "gray", alpha = 0.6, size = 1) +  # Lines use jittered x
    geom_jitter(aes(x = jittered_x), width = 0, alpha = 1, size = 2) +  # Jittered points
    theme_bw() +
    labs(
      title = "Wilcoxon paired rank test",
      x = "",
      y = "predicted DPI",
      color = "SubjectID"
    ) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.text.x = element_text(size = 20, face = "bold", colour = "black"),
      axis.title.x = element_text(size = 20, face = "bold"),
      plot.title = element_text(size = 15, face = "bold")
    ) +
    ylim(10, 90) +
    annotate("text", x = 1.5, y = 85, label = paste("p=", format(pval$p.value, digits = 3)), size = 5)
)
ggsave(filename = "plot_Ag.spec_Z7A5.png")
## output for Figure 6B

## End of project