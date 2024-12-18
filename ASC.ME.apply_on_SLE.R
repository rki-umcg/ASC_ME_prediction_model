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
# this scripts contains the application of the ASC-ME model to samples from a longitudinal observation study
# in systemic lupus erythematosus (SLE) patients
# study ID = Z6YW (flowrepository)

### This R-script requires the following packages:
library(randomForest)
library(mice)
library(readxl)
library(dplyr)
library(tidyverse)

## contine with R environment from training the RF model
load(".../ASC_ME_prediction_model/env_train_model_final.RData")
setwd(".../ASC_ME_prediction_model")

#set seed for imputation reproducibility and import SLE data
set.seed(555)
data_new <- read_excel("all_extended_data.xlsx") %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,No,CD86)) %>% 
  filter(StudyID=="Z6YW") %>%
  filter(!(condition=="pSLE_T6")) %>%
  filter(!(condition=="pSLE_T6R"))

#use imputation to complete missing validation data
#the whole imputed training data is added for this purpose to cope with new data that miss markers entirely
new_merging <- data_new %>% 
  select(sample,check,SubjectID,File,condition)
new_imp <- data_new %>% 
  select(-c(check,SubjectID,File,condition))
percent_miss<-(sum(is.na(new_imp)) / (nrow(new_imp) * ncol(new_imp))) *100
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
###use model on all imputed SLE data and average out predictions
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

### generate plots to explore the ASC-ME model application on SLE patient data

### maturity prediction by SLE disease group
AID_data<-new_output %>%
  filter(!(group=="anchor")) %>%
  filter((StudyID=="Z6YW")) %>%
  filter(str_detect(sample, "_A_"))
ANOVA_AID<- lm(AID_data[[7]] ~ AID_data[[5]])  |> aov() |> TukeyHSD()
plot_AID_own<-ggplot(AID_data, aes(x = group, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(AID_data[[2]])), width = 0.25, alpha = 1, size = 2, shape=16) +
  theme_bw() +
  labs(#title = "one-way ANOVA with Tukey correction",
       x = "",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  annotate("text", x="HC", y=55, label=paste("p=",format(ANOVA_AID[["AID_data[[5]]"]][1,4], digits=3)), size=5)+
  annotate("text", x="SLE_LN", y=55, label=paste("p=",format(ANOVA_AID[["AID_data[[5]]"]][3,4], digits=3)), size=5)+
  xlim("HC","SLE","SLE_LN")
print(plot_AID_own) 
ggsave(filename = "plot_SLE_baseline.png", plot = plot_AID_own)
## output for Figure 6D

### maturity prediction in SLE patients during study follow-up
# import SLE patient characteristics
SLE_scores <- read_excel(".../ASC_ME_prediction_model/SLE_Z6YW_stats_sorted.xlsx")
# merge SLE data and generate plot & statistics
SLE_merged<-merge(new_output, SLE_scores, by="sample", all.x = TRUE)# %>% filter(MMF_time>0)
r_AID_plot<-cor.test(SLE_merged$p, SLE_merged$Timepoint, method = "pearson", exact = FALSE)
plot_AID_example<-ggplot(SLE_merged, aes(x = Timepoint, y = p)) +
  geom_jitter(width = 0.1, alpha = 1, size = 3) +
  theme_bw() +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = paste("pearson r =",format(r_AID_plot$estimate, digits=3)," p =",format(r_AID_plot$p.value, digits=3)),
       x = "Timepoint",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  scale_x_continuous(breaks = seq(1,7, by=1))
print(plot_AID_example) 
ggsave(filename = "plot_SLE_timeline.png", plot = plot_AID_example)
## output for Figure 6D

### maturity prediction in SLE patients during study follow-up
r_AID_plot<-cor.test(SLE_merged[SLE_merged$Steroid>0,7], SLE_merged[SLE_merged$Steroid>0,25], method = "pearson", exact = FALSE)
plot_AID_example<-ggplot(SLE_merged, aes(x = Steroid, y = p)) +
  geom_jitter(width = 0.2, alpha = 1, size = 3) +
  theme_bw() +
  geom_smooth(method = "lm", color = "red", se = TRUE, 
              data = SLE_merged %>% filter(Steroid > 0), 
              aes(x = Steroid, y = p)) +
  labs(title = paste("pearson r =",format(r_AID_plot$estimate, digits=3)," p =",format(r_AID_plot$p.value, digits=3)),
       x = "Steroid usage intensity",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  xlim(-0.2,3.2)
print(plot_AID_example) 
ggsave(filename = "plot_SLE_steroids.png", plot = plot_AID_example)
## output for Figure 6D

## End of project