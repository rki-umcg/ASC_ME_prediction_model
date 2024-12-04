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
# this scripts contains the code to generate plots in Extended Figure 4
# 

### This R-script requires the following packages:
library(dplyr)
library(ggplot2)

## contine with R environment from training the RF model
load(".../ASC_ME_prediction_model/env_train_model_final.RData")
setwd(".../ASC_ME_prediction_model")

# Bin the 'predict' values into sections of width 5 and calculate the midpoint of each bin
binned_data <- cbind(data.imputed_mean,predict=predictions_mean$p)
binned_data$Age <- as.numeric(binned_data$Age)
binned_data <- binned_data %>%
  mutate(predict_bin = cut(predict, breaks=seq(5, 125, by=5), right=FALSE),
         predict_mid = as.numeric(as.character(cut(predict, breaks=seq(5, 125, by=5), labels=seq(7.5, 122.5, by=5)))))
binned_data$Age<-as.numeric(binned_data$Age)

# Plot boxplots with the actual data points plotted at their true x-values
for (a in 3:11) {
print(ggplot(binned_data, aes(x=predict, y=binned_data[,a])) +
        geom_boxplot(aes(x=predict_mid, group=predict_bin), outlier.shape = NA, width=4.5, color="red") +  # Box plots at midpoints
        geom_point(alpha=0.5, shape=1, size=1) +  # Points at actual x-values
        theme_bw() +
        labs(x="predict DPI",
             y=paste(colnames(binned_data[a]), "(normalized mean expression)")) +
        #ylim(-1,5)+
        theme(axis.title.y = element_text(size = 20,face = "bold"),
              axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
              axis.title.x = element_text(size = 20,face = "bold"),
              axis.text.x = element_text(size = 20,face = "bold", colour = "black", angle = 0)))  # Rotate x-axis labels if necessary
ggsave(filename = paste("plot_",colnames(binned_data[a]),".by.predict_low.png", sep = ""))
}

## End of project