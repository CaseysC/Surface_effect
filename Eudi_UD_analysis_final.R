#########################
# Author: Celine Caseys
# Orcid: 0000-0003-4187-9018
# Dept. of Plant sciences, University of California Davis
# July 2024

# Associated to the manuscript:
# "Leaf abaxial and adaxial surfaces differentially affect 
#  the interaction of Botrytis cinerea across several eudicots"

## Warning:
# This was coded for and functional with R version:
# R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
# Platform: x86_64-apple-darwin15.6.0 (64-bit)

##########################

library(tidyverse)
library(ggpubr)

#### This is the code for the analysis of the infection of 16 eudicot species with Botrytis cinerea

# Loading the raw data
Eudi_UD <- read.table(file="Eudi16_72h_all1.txt", header=T)
NStomata <- read.table(file="Stomata_forPlot.txt", header=T)
Stom_col <- read.table(file="Eudicot_Stomata_counts_col.txt", header=T)

# Calculate the summary statistics across replicates
Summary_spsurIso <- summarise(group_by(Eudi_UD, SpSurfIso), 
                              mean.sp=mean(Lesion_mm2), 
                              median.sp=median(Lesion_mm2), 
                              sd.sp=sd(Lesion_mm2), 
                              min.sp=min(Lesion_mm2), 
                              max.sp=max(Lesion_mm2), 
                              length.sp=length(Lesion_mm2), 
                              cv.sp=sd.sp/mean.sp)

Eudi_UD <- merge(Eudi_UD, Summary_spsurIso, by="SpSurfIso")

# Create dataframe for each species
EudicotC_data_Sp <- split(Eudi_UD, Eudi_UD$Species)
EudicotC_data_Sp  <- as.array(EudicotC_data_Sp)

######################################
### Data Cleaning 
## For details on the methodology check Caseys et al. 2021 DOI: 10.1093/g3journal/jkab175

## In short:
## To partition biological from technical failures, the lesion area distribution was analyzed 
## for each species and empirical thresholds were fixed (see below, median and size thresholds).
## A lesion below the size threshold threshold was considered a technical error only if the median 
## of lesion area for a plant genotype - strain pair was larger than the threshold. 
## The rationale is the following: when most lesions are of small size, the likelihood of biological 
## reasons for such small lesion areas is high, while when the majority of lesion areas are large, 
## the likelihood of technical error is high.
######################################


# Plotting the distribution of lesion area for each species
# These distributions inform on the size threshold

plot_list <- list()
pdf(file="LSsize_density_mm2_Eudi16.pdf")
for (i in c(1:16)) {
  plot_list[[i]] <- ggplot(EudicotC_data_Sp[[i]], aes(Lesion_mm2)) + 
    geom_density() + 
    xlim(0,50)+ 
    ggtitle(dimnames(EudicotC_data_Sp)[[1]][i])
  print(plot_list[[i]])
}
dev.off()


# Plotting the distribution of median of lesion area for each species
# These distributions inform on the median threshold

plot_list <- list()
pdf(file="Median_density_50_Eudi16.pdf")
for (i in c(1:16)) {
  plot_list[[i]] <- ggplot(EudicotC_data_Sp[[i]], aes(median.sp)) +
    geom_density() + 
    xlim(0,50) + 
    ggtitle(dimnames(EudicotC_data_Sp)[[1]][i])
  print(plot_list[[i]])
}
dev.off()


splt_Data_cleaning_filterOutliers <- list()
Failed_lesions <- list()


# Arabidopsis [1]
EudicotC_data_Sp[[1]] <- mutate(EudicotC_data_Sp[[1]],
                                median_logic=median.sp > 15,
                                threshold_logic=Lesion_mm2<15,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[1]] <- EudicotC_data_Sp[[1]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[1]] <- EudicotC_data_Sp[[1]] %>% filter(TBC_param=="TRUETRUE")



# Basil [2]
EudicotC_data_Sp[[2]] <- mutate(EudicotC_data_Sp[[2]],
                                  median_logic=median.sp > 15,
                                  threshold_logic=Lesion_mm2<13,
                                  TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                  count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[2]] <- EudicotC_data_Sp[[2]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[2]] <- EudicotC_data_Sp[[2]] %>% filter(TBC_param=="TRUETRUE")


# Chard [3]

EudicotC_data_Sp[[3]] <- mutate(EudicotC_data_Sp[[3]],
                                median_logic=median.sp > 13,
                                threshold_logic=Lesion_mm2<12,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[3]] <- EudicotC_data_Sp[[3]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[3]] <- EudicotC_data_Sp[[3]] %>% filter(TBC_param=="TRUETRUE")

# C.intybus [4]

EudicotC_data_Sp[[4]] <- mutate(EudicotC_data_Sp[[4]],
                                median_logic=median.sp > 17,
                                threshold_logic=Lesion_mm2<17,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[4]] <- EudicotC_data_Sp[[4]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[4]] <- EudicotC_data_Sp[[4]] %>% filter(TBC_param=="TRUETRUE")

# Cowpea (Vigna) [5]

EudicotC_data_Sp[[5]] <- mutate(EudicotC_data_Sp[[5]],
                                median_logic=median.sp > 17,
                                threshold_logic=Lesion_mm2<17,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[5]] <- EudicotC_data_Sp[[5]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[5]] <- EudicotC_data_Sp[[5]] %>% filter(TBC_param=="TRUETRUE")


# Cucumber [6]

EudicotC_data_Sp[[6]] <- mutate(EudicotC_data_Sp[[6]],
                                median_logic=median.sp > 15,
                                threshold_logic=Lesion_mm2<13,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[6]] <- EudicotC_data_Sp[[6]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[6]] <- EudicotC_data_Sp[[6]] %>% filter(TBC_param=="TRUETRUE")


# Eggplant [7]


EudicotC_data_Sp[[7]] <- mutate(EudicotC_data_Sp[[7]],
                                median_logic=median.sp > 16,
                                threshold_logic=Lesion_mm2<20,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[7]] <- EudicotC_data_Sp[[7]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[7]] <- EudicotC_data_Sp[[7]] %>% filter(TBC_param=="TRUETRUE")

# Kale [8]

EudicotC_data_Sp[[8]] <- mutate(EudicotC_data_Sp[[8]],
                                median_logic=median.sp > 14,
                                threshold_logic=Lesion_mm2<14,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[8]] <- EudicotC_data_Sp[[8]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[8]] <- EudicotC_data_Sp[[8]] %>% filter(TBC_param=="TRUETRUE")

# Lettuce [9]

EudicotC_data_Sp[[9]] <- mutate(EudicotC_data_Sp[[9]],
                                median_logic=median.sp > 13,
                                threshold_logic=Lesion_mm2<13,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[9]] <- EudicotC_data_Sp[[9]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[9]] <- EudicotC_data_Sp[[9]] %>% filter(TBC_param=="TRUETRUE")

# Mint [10]

EudicotC_data_Sp[[10]] <- mutate(EudicotC_data_Sp[[10]],
                                median_logic=median.sp > 10,
                                threshold_logic=Lesion_mm2<11,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[10]] <- EudicotC_data_Sp[[10]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[10]] <- EudicotC_data_Sp[[10]] %>% filter(TBC_param=="TRUETRUE")


# Parsley [11]

EudicotC_data_Sp[[11]] <- mutate(EudicotC_data_Sp[[11]],
                                median_logic=median.sp > 15,
                                threshold_logic=Lesion_mm2<16,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[11]] <- EudicotC_data_Sp[[11]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[11]] <- EudicotC_data_Sp[[11]] %>% filter(TBC_param=="TRUETRUE")

# Pepper [12]

EudicotC_data_Sp[[12]] <- mutate(EudicotC_data_Sp[[12]],
                                median_logic=median.sp > 10,
                                threshold_logic=Lesion_mm2<10,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[12]] <- EudicotC_data_Sp[[12]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[12]] <- EudicotC_data_Sp[[12]] %>% filter(TBC_param=="TRUETRUE")

# Phaseolus [13]

EudicotC_data_Sp[[13]] <- mutate(EudicotC_data_Sp[[13]],
                                median_logic=median.sp > 12,
                                threshold_logic=Lesion_mm2<12,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[13]] <- EudicotC_data_Sp[[13]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[13]] <- EudicotC_data_Sp[[13]] %>% filter(TBC_param=="TRUETRUE")

# Raphanus [14]

EudicotC_data_Sp[[14]] <- mutate(EudicotC_data_Sp[[14]],
                                median_logic=median.sp > 13,
                                threshold_logic=Lesion_mm2<13,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[14]] <- EudicotC_data_Sp[[14]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[14]] <- EudicotC_data_Sp[[14]] %>% filter(TBC_param=="TRUETRUE")

# Spinach [15]
EudicotC_data_Sp[[15]] <- mutate(EudicotC_data_Sp[[15]],
                                median_logic=median.sp > 17,
                                threshold_logic=Lesion_mm2<17,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[15]] <- EudicotC_data_Sp[[15]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[15]] <- EudicotC_data_Sp[[15]] %>% filter(TBC_param=="TRUETRUE")

# Tomato [16]

EudicotC_data_Sp[[16]] <- mutate(EudicotC_data_Sp[[16]],
                                median_logic=median.sp > 16,
                                threshold_logic=Lesion_mm2<15,
                                TBC_param=as.factor(paste(median_logic, threshold_logic, sep="")),
                                count_TBC=length(Lesion_mm2[TBC_param=="TRUETRUE"]))

# Select all data points that are not 'TRUETRUE'
splt_Data_cleaning_filterOutliers[[16]] <- EudicotC_data_Sp[[16]] %>% filter(TBC_param!="TRUETRUE")

# Select all data points that are 'TRUETRUE', e.g. the data removed from the dataset
Failed_lesions[[16]] <- EudicotC_data_Sp[[16]] %>% filter(TBC_param=="TRUETRUE")

#############

Eudicot_UD_clean <- rbind(splt_Data_cleaning_filterOutliers[[1]], splt_Data_cleaning_filterOutliers[[2]], splt_Data_cleaning_filterOutliers[[3]], splt_Data_cleaning_filterOutliers[[4]], splt_Data_cleaning_filterOutliers[[5]], splt_Data_cleaning_filterOutliers[[6]], splt_Data_cleaning_filterOutliers[[7]], 
                          splt_Data_cleaning_filterOutliers[[8]], splt_Data_cleaning_filterOutliers[[9]], splt_Data_cleaning_filterOutliers[[10]], splt_Data_cleaning_filterOutliers[[11]], splt_Data_cleaning_filterOutliers[[12]], splt_Data_cleaning_filterOutliers[[13]], splt_Data_cleaning_filterOutliers[[14]], splt_Data_cleaning_filterOutliers[[15]], splt_Data_cleaning_filterOutliers[[16]])
write.table(Eudicot_UD_clean, file="Eudicot16_UD_clean.txt")

Names <- c("Brapa", "Cendivia", "Cintybus", "Glycine", "Helianthus", "Lactuca", "Solanum")

#FailedLesions_UD <- rbind(Failed_lesions[[1]], Failed_lesions[[2]], Failed_lesions[[3]], Failed_lesions[[4]],Failed_lesions[[5]],Failed_lesions[[6]],Failed_lesions[[7]],
#                                  Failed_lesions[[8]], Failed_lesions[[9]], Failed_lesions[[10]], Failed_lesions[[11]], Failed_lesions[[12]], Failed_lesions[[13]], Failed_lesions[[14]], Failed_lesions[[15]], Failed_lesions[[16]])
#write.table(FailedLesions_UD, file="FailedLesions_Eudicot15_Ara_UD.txt")

##################################################

Eudi_UD_clean <- read.table(file="Eudicot16_UD_clean.txt", header=T)

Summary_spsurIso_clean <- summarise(group_by(Eudi_UD_clean, SpSurfIso), 
                             mean.sp=mean(Lesion_mm2))

Summary_spsurIso_clean = separate(Summary_spsurIso_clean, col=SpSurfIso, into=c("Species", "Surface", "Isolate"), sep="\\_")
Summary_spsurIso_clean$Species <- as.factor(Summary_spsurIso_clean$Species)
Summary_spsurIso_clean$Surface <- as.factor(Summary_spsurIso_clean$Surface)
Summary_spsurIso_clean$Isolate <- as.factor(Summary_spsurIso_clean$Isolate)

Surface_diff_data <- pivot_wider(Summary_spsurIso_clean, names_from = "Surface", values_from = "mean.sp")
Surface_diff_data <- mutate(Surface_diff_data, Surface_effect=Down-Up )


Summary_spsur_clean <- summarise(group_by(Eudi_UD_clean, Sp_Surface), 
                                    mean.sp=mean(Lesion_mm2), 
                                    median.sp=median(Lesion_mm2), 
                                    sd.sp=sd(Lesion_mm2), 
                                    min.sp=min(Lesion_mm2), 
                                    max.sp=max(Lesion_mm2), 
                                    length.sp=length(Lesion_mm2), 
                                    cv.sp=sd.sp/mean.sp)



#write.table(Summary_spsurIso_clean, file="Eudicot15_72_Clean_summary.txt")



Eudi_UD_clean_Phylum <- split(Eudi_UD_clean, Eudi_UD_clean$Phylum)


# plot_list <- list()
# pdf(file="Boxplot_Eudi15.pdf")
# for (i in c(1:3)) {
#   plot_list[[i]] <- ggplot(data=Eudi_UD_clean_Phylum[[i]], aes(x=Order_Plot, y=Lesion_mm2, fill=Species)) +
#     geom_boxplot()+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   print(plot_list[[i]])
# }
# dev.off()


EudicotC_data_Phy <- split(Eudi_UD, Eudi_UD$Phylum)
list2env(EudicotC_data_Phy, envir = .GlobalEnv)

Asterids <- droplevels.data.frame(Asterids)
Rosids <- droplevels.data.frame(Rosids)
Core <- droplevels.data.frame(Core)

Eudi_UD$Order_Plot_final <- Eudi_UD$Order_Plot
 

# levels(Eudi_UD$Order_Plot_final)
#  "Asterids_Apiales_Par_Down"     "Asterids_Apiales_Par_Up"       "Asterids_Asterales_Cin_Down"   "Asterids_Asterales_Cin_Up"    
#  "Asterids_Asterales_Lac_Down"   "Asterids_Asterales_Lac_Up"     "Asterids_Lamiales_Bas_Down"    "Asterids_Lamiales_Bas_Up"     
#  "Asterids_Lamiales_Min_Down"    "Asterids_Lamiales_Min_Up"      "Asterids_Solanales_Egg_Down"   "Asterids_Solanales_Egg_Up"    
#  "Asterids_Solanales_Pep_Down"   "Asterids_Solanales_Pep_Up"     "Asterids_Solanales_Tom_Down"   "Asterids_Solanales_Tom_Up"    
#  "Core_Caryophyllales_Char_Down" "Core_Caryophyllales_Char_Up"   "Core_Caryophyllales_Spi_Down"  "Core_Caryophyllales_Spi_Up"   
#  "Rosids_Brassicales_Ara_Down"   "Rosids_Brassicales_Ara_Up"     "Rosids_Brassicales_Kale_Down"  "Rosids_Brassicales_Kale_Up"   
#  "Rosids_Brassicales_Rap_Down"   "Rosids_Brassicales_Rap_Up"     "Rosids_Cucurbitales_Cuc_Down"  "Rosids_Cucurbitales_Cuc_Up"   
#  "Rosids_Fabales_Cp_Down"        "Rosids_Fabales_Cp_Up"          "Rosids_Fabales_Pha_Down"       "Rosids_Fabales_Pha_Up"        
#  "Rosids_Fabales_Rap_Up"  

# Ordering the species by phylogenetic order
levels(Eudi_UD$Order_Plot_final) <-  c("03_Par_Down", "03_Par_Up", "04_Cin_Down","04_Cin_Up", "05_Lac_Down", "05_Lac_Up",  "06_Bas_Down", "06_Bas_Up", "07_Min_Down" ,  
                                     "07_Min_Up", "08_Egg_Down", "08_Egg_Up", "10_Pep_Down","10_Pep_Up", "09_Tom_Down", "09_Tom_Up","01_Char_Down", "01_Char_Up",  
                                     "02_Spi_Down", "02_Spi_Up", "13_Ara_Down","13_Ara_Up", "11_Kale_Down",  "11_Kale_Up", "12_Rap_Down", "12_Rap_Up","14_Cuc_Down",
                                    "14_Cuc_Up", "15_Cp_Down", "15_Cp_Up", "16_Pha_Down", "16_Pha_Up","12_Rap_Up")

Eudi_UD$Order_Plot_final <- as.factor(as.character(Eudi_UD$Order_Plot_final))


## Figure S1A
ggplot(data=Eudi_UD, aes(x=Order_Plot_final, y=Lesion_mm2, color=Species)) +
  geom_boxplot(width=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data=Asterids, aes(x=Order_Plot, y=Lesion_mm2, color=Species)) +
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Asterids")
 

ggplot(data=Rosids, aes(x=Order_Plot, y=Lesion_mm2, color=Species)) +
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Rosids")

ggplot(data=Core, aes(x=Order_Plot, y=Lesion_mm2, color=Species)) +
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Caryophyllales")


#################################
#Statistical models 
################################
library(MASS)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)

Eudi_UD_clean_Species <- split(Eudi_UD_clean, Eudi_UD_clean$Species)


UD_model <- lmer(Lesion_mm2 ~Species*Isolate*Surface + (1|Tray), data=Eudi_UD_clean)
model_anova<- Anova(UD_model)


Anova_species_list <- list()

NamesPlotN <- c("13_Ara","06_Bas", "01_Char", "04_Cin", "15_Cp", "14_Cuc", "08_Egg","11_Kale", "05_Lac", "07_Min","03_Par","10_Pep","16_Pha", "12_Rap","02_Spi","09_Tom") 


for (i in c(1:16)) {
Lm1 <- lm(Lesion_mm2 ~ Isolate*Surface + Tray, data=Eudi_UD_clean_Species[[i]])
Alm1 <- anova(Lm1)
Alm1$PercVar <- (Alm1$`Sum Sq`*100)/sum(Alm1$`Sum Sq`)
Alm1$Species <- rep(NamesPlotN[i], 5)
Alm1$Param <- rownames(Alm1)
Anova_species_list[[i]] <- as.data.frame(Alm1)
}


ModelRes <- rbind(Anova_species_list[[1]], Anova_species_list[[2]])
for (i in c(3:15)) {ModelRes <- rbind(ModelRes, Anova_species_list[[i]])}

#write.table(ModelRes ,file="ModelRes_Eudi16_UD_Clean.txt",sep="\t")
#ModelRes <- read.table(file="ModelRes_Eudi16_UD_Clean.txt", header = T)

ModelRes$PercVar1 <- round(ModelRes$PercVar, digits=1)
ModelRes$Param <- as.factor(ModelRes$Param)
levels(ModelRes$Param) <- c("1_Isolate", "3_Isolate:Surface", "5_Residuals", "2_Surface", "4_Tray" )
ModelRes$Param1 <- as.factor(as.character(ModelRes$Param))


## Figure S1B
ggplot(ModelRes, aes(Species, PercVar, fill=Param1)) +
  geom_col()+
  scale_fill_manual(values=c("#f1b6da", "#7b3294", "#92c5de","#bababa","#ffffbf"))+
  geom_text(aes(label=PercVar1), position=position_stack(vjust=0.5), size=3)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##############
## Least-square means
##############

library(emmeans)

Eudi_UD_clean_Sp <- split(Eudi_UD_clean, Eudi_UD_clean$Species)
Eudi_UD_clean_Sp  <- as.array(Eudi_UD_clean_Sp)

Sp_names <- names(Eudi_UD_clean_Sp)

emmeans_list <- list()

Lesion_LS <- as.data.frame(matrix(ncol=7))
colnames(Lesion_LS) <- c("Surface", "emmean", "SE" , "df", "lower.CL", "upper.CL","species" )

for (i in c(1:16)) {
  Species_model <- lmer(Lesion_mm2 ~ Surface + (1|Isolate) + (1|Tray), data=Eudi_UD_clean_Sp[[i]]) 
  tmp <- emmeans(Species_model, ~Surface, lmer.df='satterthwaite')
  tmp <- as.data.frame(print(tmp))
  tmp$species <-  rep(Sp_names[[i]],2)
  emmeans_list[[i]] <- tmp
  Lesion_LS <- rbind(Lesion_LS, tmp)
}

Lesion_LSs <-  Lesion_LS[-1,c(1,2,7)]
Surface_diff_LSM <- pivot_wider(Lesion_LSs, names_from = "Surface", values_from = "emmean")



LesionSp_LS <- data.frame()
  Species_model <- lmer(Lesion_mm2 ~ Species*Surface + (1|Isolate) + (1|Tray), data=Eudi_UD_clean) 
  tmp <- emmeans(Species_model, ~Species, lmer.df='satterthwaite')
  LesionSp_LS <- as.data.frame(print(tmp))


colnames(LesionSp_LS) <- c( "species", "emmean_sp", "SE_sp", "df_Sp", "lower.CL_sp", "upper.CL_sp")



Surface_diff_LSM$diff <- Surface_diff_LSM$Down-Surface_diff_LSM$Up

Surface_diff_LSM <- merge(Surface_diff_LSM, LesionSp_LS, by="species")

Surface_diff_LSM$species_plot <- c("13_Ara","06_Bas", "01_Char", "04_Cin", "15_Cp", "14_Cuc", "08_Egg","11_Kale", "05_Lac", "07_Min","03_Par","10_Pep","16_Pha", "12_Rap","02_Spi","09_Tom") 

Surface_diff_LSM$Std_Diff <- (Surface_diff_LSM$Down-Surface_diff_LSM$Up)/((Surface_diff_LSM$Down+Surface_diff_LSM$Up)/2)

Surface_diff_LSM$Leaf_Thickness <- c("NA","226", "371", "232", "NA", "293", "310","244", "120", "260","182","277","386", "427","323","324")
Surface_diff_LSM$Leaf_Thickness <- as.numeric(Surface_diff_LSM$Leaf_Thickness)


NStomata <- read.table(file="Stomata_forPlot.txt", header=T)

Surface_diff_LSM <- merge(Surface_diff_LSM, NStomata, by="species_plot")

Eudi16_color <- c("#8c510a", "#9e823e", "#c4674b", "#e69705", "#dba502", "#5f8a2c", "#0d663c", "#bd0026", "#e31a1c", "#f56f53", "#16ab94", "#35978f", "#0c6e66", "#e7298a", "#c51b7d", "#9e4f79") 
#                   Chard       Spinach  Parsley    Cin         Lac         Basil   Mint        Eggplant   TOMATO     Pepper      Kale      Radish    Arabido     Cucumber    Cowpea    Bean



ggplot(Surface_diff_LSM, aes(species_plot, Std_Diff, color=species_plot)) +
  geom_point(size=2.5, show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Standardized Surface Differential")

  
  
ggplot(Surface_diff_LSM, aes(species_plot, emmean_sp, color=species_plot)) +
    geom_point(size=2.5, show.legend = FALSE)+
    scale_color_manual(values=Eudi16_color)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ylab("Lesion Area [mm2]")
  
  
ggplot(Surface_diff_LSM, aes(species_plot, Stom_Diff, color=species_plot)) +
    geom_point(size=2.5, show.legend = FALSE)+
    scale_color_manual(values=Eudi16_color)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ylab("Standardized Stomata Differential")
  



#######################
# Mean centering lesion data
###########
colnames(LesionSp_LS)[1] <- c("Species")

Eudi_UD_clean1 <- merge(Eudi_UD_clean, LesionSp_LS, by="Species", all.x=T)
Eudi_UD_clean1$Lesion_mm2_centered <- Eudi_UD_clean1$Lesion_mm2-Eudi_UD_clean1$emmean_sp


Eudi_UD_clean_Sp <- split(Eudi_UD_clean1, Eudi_UD_clean1$Species)
Eudi_UD_clean_Sp  <- as.array(Eudi_UD_clean_Sp)

Sp_names <- names(Eudi_UD_clean_Sp)

emmeans_list <- list()

Lesion_LS_centered <- as.data.frame(matrix(ncol=7))
colnames(Lesion_LS_centered) <- c("Surface", "emmean", "SE" , "df", "lower.CL", "upper.CL","species" )


for (i in c(1:16)) {
  Species_model <- lmer(Lesion_mm2_centered ~ Surface + (1|Isolate) + (1|Tray), data=Eudi_UD_clean_Sp[[i]]) 
  tmp <- emmeans(Species_model, ~Surface, lmer.df='satterthwaite')
  tmp <- as.data.frame(print(tmp))
  tmp$species <-  rep(Sp_names[[i]],2)
  emmeans_list[[i]] <- tmp
  Lesion_LS_centered <- rbind(Lesion_LS_centered, tmp)
}


Lesion_LS_centered <-    Lesion_LS_centered[-1,]


Species_names_plot <- Surface_diff_LSM[, c(1,2)]
Lesion_LS_centered <- merge(Lesion_LS_centered, Species_names_plot, by="species")
Lesion_LS_centered$species_plot <- as.factor(as.character(Lesion_LS_centered$species_plot))

Lesion_LS_centered$SpSurf <- paste(Lesion_LS_centered$species_plot, Lesion_LS_centered$Surface, sep="_")
Lesion_LS_centered$SpSurf <- as.factor(Lesion_LS_centered$SpSurf)

d <-  ggplot(Lesion_LS_centered, aes(SpSurf, emmean, color=species_plot, shape=Surface)) +
  geom_point(size=3, show.legend = FALSE)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.8, position = position_dodge(0.1), , show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Mean-centered Lesion Area")



Stom_col <- read.table(file="Eudicot_Stomata_counts_col.txt", header=T)
Stom_col$SpSurf <- paste(Stom_col$species_plot, Stom_col$Surface, sep="_")

Mean_Stom_spe <- summarise(group_by(Stom_col, species_plot), 
                              mean.stom=mean(N_stomata))
Stom_col <- merge(Stom_col, Mean_Stom_spe, by="species_plot", all.x=T)
Stom_col$N_Stom_centered <- Stom_col$N_stomata-Stom_col$mean.stom

std_E <- function(x) sd(x)/sqrt(length(x))

Summary_Stom_col <- summarise(group_by(Stom_col, SpSurf), 
                              mean.stom=mean(N_Stom_centered), 
                              median.stom=median(N_Stom_centered), 
                              sd.stom=sd(N_Stom_centered), 
                              se.stom=std_E(N_Stom_centered),
                              min.stom=min(N_Stom_centered), 
                              max.stom=max(N_Stom_centered), 
                              length.stom=length(N_Stom_centered))

Summary_Stom_col$SpSurf1 <- Summary_Stom_col$SpSurf
Summary_Stom_col <- separate(Summary_Stom_col, col=SpSurf1, into = c("N", "Sp", "Surface"), sep="\\_")


Summary_Stom_col$species_plot <- paste(Summary_Stom_col$N, Summary_Stom_col$Sp, sep="_")

e <- ggplot(Summary_Stom_col, aes(SpSurf, mean.stom, color=species_plot, shape=Surface)) +
  geom_point(size=3, show.legend = FALSE)+
  geom_errorbar(aes(ymin=mean.stom-se.stom, ymax=mean.stom+se.stom), width=.8, position = position_dodge(0.1), , show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Mean-centered Stomata Density")


## Figure 2
ggarrange(d, e, nrow=2, ncol=1)


Anova_species_stom <- list()

Stom_col_species <- split(Stom_col, Stom_col$Species_Full)

for (i in c(1:16)) {
  Lm1 <- lm(N_stomata ~ Surface, data=Stom_col_species[[i]])
  Alm1 <- anova(Lm1)
  Alm1$Param <- rownames(Alm1)
  Alm1$PercVar <- (Alm1$`Sum Sq`*100)/sum(Alm1$`Sum Sq`)
  Alm1$Param <- rownames(Alm1)
  Alm1$Species <- rep(names(Stom_col_species)[i], dim(Alm1)[1])
  Anova_species_stom[[i]] <- as.data.frame(Alm1)
}

ModelStom <- rbind(Anova_species_stom[[1]], Anova_species_stom[[2]])
for (i in c(3:16)) {ModelStom <- rbind(ModelStom, Anova_species_stom[[i]])}

#write.table(ModelStom ,file="LMModelStomata.txt",sep="\t"


##################
# For R>4.0 so on windows computer
## Not working even on the windows computer

### Supplemental
p001 <- ggplot(Surface_diff_LSM, aes(emmean_sp, diff, color=species_plot)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point(size=2.5, show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  xlab("Lesion area")+
  ylab("Surface differential")+
  annotate("text", x=65, y=35, label="R2adj=0.1103, p-val=0.113")


Lm001 <-  lm(diff~emmean_sp, data=Surface_diff_LSM)
summary(Lm001)


p002 <-ggplot(Surface_diff_LSM, aes(emmean_sp, Leaf_Thickness, color=species_plot)) +
  geom_smooth(method = "lm", se=T, color="black", formula = y ~ x) +
  #stat_poly_eq(use_label(c("adj.R2", "p"))) +
  geom_point(size=2.5, show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  xlab("Lesion Area")+
  ylab("Leaf thickness")+
  annotate("text", x=65, y=400, label="R2adj=0.01595, p-val=0.2978")

Lm002 <-  lm(Leaf_Thickness~emmean_sp, data=Surface_diff_LSM)
summary(Lm002)


p003 <- ggplot(Surface_diff_LSM,aes(diff, Leaf_Thickness, color=species_plot)) +
  geom_smooth(method = "lm", se=T, color="black", formula = y ~ x) +
  #stat_poly_eq(use_label(c("adj.R2", "p"))) +
  geom_point(size=2.5, show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  xlab("Surface differential")+
  ylab("Leaf thickness")+
  annotate("text", x=15, y=400, label="R2adj=-0.07804, p-val=0.7239")

Lm003 <-  lm(Leaf_Thickness~diff, data=Surface_diff_LSM)
summary(Lm003)


p004 <- ggplot(Surface_diff_LSM, aes(Up, X_stom_adaxial, color=species_plot)) +
  geom_smooth(method = "lm", se=T, color="black", formula = y ~ x) +
  #stat_poly_eq(use_label(c("adj.R2", "p"))) +
  geom_point(size=2.5, show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  xlab("Adaxial Lesion area")+
  ylab("Adaxial stomata density")+
  annotate("text", x=65, y=30, label="R2adj=-0.07121, p-val=0.9586")

Lm004 <-  lm(X_stom_adaxial~Up, data=Surface_diff_LSM)
summary(Lm004)


p005 <- ggplot(Surface_diff_LSM, aes(Up, Leaf_Thickness, color=species_plot)) +
  geom_smooth(method = "lm", se=T, color="black", formula = y ~ x) +
  #stat_poly_eq(use_label(c("adj.R2", "p"))) +
  geom_point(size=2.5, show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  xlab("Abaxial Lesion area")+
  ylab("Leaf thickness")+
  annotate("text", x=65, y=400, label="R2adj=0.04278, p-val=0.241")

Lm005 <-  lm(Leaf_Thickness~Up, data=Surface_diff_LSM)
summary(Lm005)



p006 <- ggplot(Surface_diff_LSM, aes(Down, Leaf_Thickness, color=species_plot)) +
  geom_smooth(method = "lm", se=T, color="black", formula = y ~ x) +
  #stat_poly_eq(use_label(c("adj.R2", "p"))) +
  geom_point(size=2.5, show.legend = F)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  xlab("Abaxial Lesion area")+
  ylab("Leaf thickness")+
  annotate("text", x=80, y=400, label="R2adj=-0.01107, p-val=0.3713")

Lm006 <-  lm(Leaf_Thickness~Down, data=Surface_diff_LSM)
summary(Lm006)


p007 <- ggplot(Surface_diff_LSM, aes(Down, stom_abaxial, color=species_plot)) +
  geom_smooth(method = "lm", se=T, color="black", formula = y ~ x) +
  #stat_poly_eq(use_label(c("adj.R2", "p"))) +
  geom_point(size=2.5, show.legend = FALSE)+
  scale_color_manual(values=Eudi16_color)+
  theme_minimal()+
  xlab("Abaxial Lesion area")+
  ylab("Abaxial stomata density")+
  annotate("text", x=80, y=100, label="R2adj=-0.05011, p-val=0.6023")

Lm007 <-  lm(stom_abaxial~Down, data=Surface_diff_LSM)
summary(Lm007)



ggarrange(p001, p002, p004, p003, p007, p005, p006, nrow=4, ncol=2)
