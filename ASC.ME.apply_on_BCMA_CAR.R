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

#set seed for imputation reproducibility and import BCMA CAR-T cell data
set.seed(555)
data_new <- read_excel("all_extended_data.xlsx") %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,No,CD86)) %>% 
  filter(StudyID=="Z4KB") 

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
###use model on all imputed BCMA CAR-T data and average out predictions
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

### generate plots to explore the ASC-ME model application 
### on a case report of BCMA CAR-T cell treatment in multiple myeloma

### maturity prediction during BCMA CAR-T treatment
plot_BCMA_CAR<-ggplot(new_output, aes(x = dpi, y = p)) +
  geom_point(shape=16, size=3) + #none-filled circles
  theme_bw() +
  labs(#title = "one-way ANOVA with Tukey correction",
       x = "days post BCMA CAR-T cell treatment",
       y = "predicted DPI") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  ylim(40,90)+
  scale_x_continuous(breaks = seq(0,160, by=20))
print(plot_BCMA_CAR) 
ggsave(filename = "plot_BCMA_CAR_timeline.png", plot = plot_BCMA_CAR)
## output for Figure 6E
