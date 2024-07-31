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

#### This is the code for the analysis of glucosinolates and camalexin concentration 
#### in infected and uninfected leaves of A.thaliana

## Glucosinolates

GSL_Data <- read.table(file="GSL_Leaf_Data1.txt", header=T)

## Camalexin

Cam_data  <- read.table(file="Cam_Leaf_Data1.txt", header=T)

GSL_Data$treatSurf <-  as.factor(paste(GSL_Data$Surface, GSL_Data$Treatment))
levels(GSL_Data$treatSurf) <- c("Abaxial-Uninfected", "Abaxial-Infected", "Adaxial-Uninfected", "Adaxial-Infected")
GSL_Data$Key <- as.factor(paste(GSL_Data$treatSurf , GSL_Data$Genotype, sep="-"))

GSL_Data <-  mutate(GSL_Data,
                    ## Concentrations were standardized to the leaf area
                          QLLf_3MSO= QL_3MSO/RemainLFcm,
                          QLLf_4MSO= QL_4MSO/RemainLFcm,
                          QLLf_5MSO= QL_5MSO/RemainLFcm,
                          QLLf_6MSO= QL_6MSO/RemainLFcm,
                          QLLf_7MSO= QL_7MSO/RemainLFcm,
                          QLLf_8MSO= QL_8MSO/RemainLFcm,
                          QLLf_4MT= QL_4MT/RemainLFcm,
                          QLLf_I3M= QL_I3M/RemainLFcm,
                          QLLf_4MI3M= QL_4MI3M/RemainLFcm,
                          QLLf_4OHI3M= QL_4OHI3M/RemainLFcm,
                    ## Short chain aliphatic glucosinolates
                          SC_GSL=QLLf_3MSO+QLLf_4MSO+QLLf_4MT,
                    ## Long chain aliphatic glucosinolates
                          LC_GSL=QLLf_7MSO+QLLf_8MSO,
                    ## Total aliphatic glucosinolates
                          aliphatic=QLLf_3MSO+QLLf_4MSO+QLLf_5MSO+QLLf_6MSO+QLLf_7MSO+QLLf_8MSO+QLLf_4MT,
                    ## Total indolic glucosinolates
                          indolic=QLLf_I3M+QLLf_4MI3M+QLLf_4OHI3M,
                    ## Proportions long/short chain aliphatic glucosinolates
                          SCProp=SC_GSL/(LC_GSL+SC_GSL),
                          LCProp=LC_GSL/(LC_GSL+SC_GSL),
                          SLratio=SC_GSL/LC_GSL,
                    ## Level of methylation of 4C glcosinolates
                          Stress=QLLf_4MSO/QLLf_4MT,
                          AlInd_ratio=aliphatic/indolic
)

GSL_Leaf_DataQ <- GSL_Data[, c(80:99)]

GSL_Leaf_Data_noCtrl <- GSL_Data[which(GSL_Data$Isolate!='Ctrl'),]
GSL_Leaf_Data_noCtrl$GenoS <- paste(GSL_Leaf_Data_noCtrl$Geno1, GSL_Leaf_Data_noCtrl$Surface1, sep="-")
GSL_Leaf_DataCtrl <- GSL_Data[which(GSL_Data$Isolate=='Ctrl'),]


###################################
##### Raw data plots
##### Figure S4 Data
###################################


ggplot(GSL_Data, aes(Genotype, QLLf_3MSO*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("3MSO, nmol/cm2")+
  theme_bw()

ggplot(GSL_Data, aes(Genotype, QLLf_4MSO*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("4MSO, nmol/cm2")+
  theme_bw()

ggplot(GSL_Data, aes(Genotype, QLLf_5MSO*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("5MSO, nmol/cm2")+
  theme_bw()

ggplot(GSL_Data, aes(Genotype, QLLf_6MSO*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("6MSO, nmol/cm2")+
  theme_bw()

ggplot(GSL_Data, aes(Genotype, QLLf_7MSO*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("7MSO, nmol/cm2")+
  theme_bw()

ggplot(GSL_Data, aes(Genotype, QLLf_4MT*1000, color=treatSurf))+
  geom_boxplot()+
ylab("4MT, nmol/cm2")+
  theme_bw()

ggplot(GSL_Data, aes(Genotype, QLLf_8MSO*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("8MSO, nmol/cm2")+
  theme_bw()


ggplot(GSL_Data, aes(Genotype, aliphatic*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("Aliphatic glucosinolates, nmol/cm2")+
  theme_bw()


ggplot(GSL_Data, aes(Genotype, indolic*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("Indolic glucosinolates, nmol/cm2")+
  theme_bw()

ggplot(GSL_Data, aes(Genotype, SC_GSL*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("Short chain glucosinolates, nmol/cm2")+
  theme_bw()


ggplot(GSL_Data, aes(Genotype, LC_GSL*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("Long chain glucosinolates, nmol/cm2")+
  theme_bw()


ggplot(GSL_Data, aes(Genotype, SCProp*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("Shor/Long glucosinolates ratio, nmol/cm2")+
  theme_bw()

ggplot(GSL_Data, aes(Genotype, LCProp*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("Shor/Long glucosinolates ratio, nmol/cm2")+
  theme_bw()


ggplot(GSL_Data, aes(Genotype, Stress, color=treatSurf))+
  geom_boxplot()+
  ylab("Long chain glucosinolates, nmol/cm2")+
  theme_bw()


ggplot(GSL_Data, aes(Genotype, QLLf_I3M*1000, color=treatSurf))+
  geom_boxplot()+
  ylab("I3M, nmol/cm2")+
  theme_bw()
 




###################################

##### Correlation to lesion area
##### Figures S6-9
##########################
GSL_Data$Lesion_mm2 <- (GSL_Data$Lesion.Size/(89.2^2))*100
Cam_data$Lesion_mm2 <- (Cam_data$Lesion.Size/(89.2^2))*100

GSL_Data_noCtrl <- GSL_Data[which(GSL_Data$Treatment!="Control"),]
Cam_Data_noCtrl <- Cam_data[which(Cam_data$Treatment!="Control"),]
Cam_Data_noCtrl <- Cam_data[which(Cam_data$Cam_LF>0),]

GSL_geno_split <-  split(GSL_Data_noCtrl,GSL_Data_noCtrl$Genotype) 
Cam_geno_split <-  split(Cam_Data_noCtrl,Cam_Data_noCtrl$Genotype) 

list2env(GSL_geno_split, env=.GlobalEnv)

### Remove the controls as they are not infected so CAN'T have lesions

## I3M (Figure S8)

### Legend: red=Abaxial Green=Adaxial
##########################

I6 <- ggplot(Col0, aes(Lesion_mm2,QL_I3M)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlim(0, 160)+
  ylim(0, 0.03)+
  ggtitle("Col-0 **")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("I3M [nmol/cm2]")


Lm001 <-  lm(QL_I3M~Lesion_mm2, data=Col0)
summary(Lm001)

I1 <- ggplot(MYB28, aes(Lesion_mm2,QL_I3M)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.03)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("MYB28")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("I3M [nmol/cm2]")


Lm001 <-  lm(QL_I3M~Lesion_mm2, data=MYB28)
summary(Lm001)

I2 <- ggplot(m2829, aes(Lesion_mm2,QL_I3M)) +
  geom_point(aes(color=Surface), show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.03)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("myb28/29")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("I3M [nmol/cm2]")

Lm001 <-  lm(QL_I3M~Lesion_mm2, data=m2829)
summary(Lm001)


I3 <- ggplot(tgg12, aes(Lesion_mm2,QL_I3M)) +
  geom_point(aes(color=Surface))+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.03)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("tgg1/2 *")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("I3M [nmol/cm2]")

Lm001 <-  lm(QL_I3M~Lesion_mm2, data=tgg12)
summary(Lm001)

I4 <- ggplot(m3451, aes(Lesion_mm2,QL_I3M)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.03)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("myb34/51")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("I3M [nmol/cm2]")

Lm001 <-  lm(QL_I3M~Lesion_mm2, data=m3451)
summary(Lm001)


I5 <- ggplot(Pad3, aes(Lesion_mm2,QL_I3M)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.03)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("pad3 **")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("I3M [nmol/cm2]")

Lm001 <-  lm(QL_I3M~Lesion_mm2, data=Pad3)
summary(Lm001)


ggarrange(I1, I2, I3, I4, I5, I6, 
          labels = c("A)", "B)", "C)", "D)", "E)", "F)"),
          ncol = 2, nrow = 3)

##########################
## 4MI3M (Figure S9)

M6 <-ggplot(Col0, aes(Lesion_mm2,QL_4MI3M)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlim(0, 160)+
  ylim(0, 0.01)+
  ggtitle("Col-0")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("4MI3M [nmol/cm2]")


Lm001 <-  lm(QL_4MI3M~Lesion_mm2, data=Col0)
summary(Lm001)

M1 <- ggplot(MYB28, aes(Lesion_mm2,QL_4MI3M)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.01)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("MYB28")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("4MI3M [nmol/cm2]")


Lm001 <-  lm(QL_4MI3M~Lesion_mm2, data=MYB28)
summary(Lm001)

M2 <- ggplot(m2829, aes(Lesion_mm2,QL_4MI3M)) +
  geom_point(aes(color=Surface), show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.01)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("myb28/29")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("4MI3M [nmol/cm2]")

Lm001 <-  lm(QL_4MI3M~Lesion_mm2, data=m2829)
summary(Lm001)


M3 <- ggplot(tgg12, aes(Lesion_mm2,QL_4MI3M)) +
  geom_point(aes(color=Surface), , show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.01)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("tgg1/2")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("4MI3M [nmol/cm2]")

Lm001 <-  lm(QL_4MI3M~Lesion_mm2, data=tgg12)
summary(Lm001)

M4 <- ggplot(m3451, aes(Lesion_mm2,QL_4MI3M)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.01)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("myb34/51")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("4MI3M [nmol/cm2]")

Lm001 <-  lm(QL_4MI3M~Lesion_mm2, data=m3451)
summary(Lm001)


M5 <- ggplot(Pad3, aes(Lesion_mm2,QL_4MI3M)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.01)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("pad3")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("4MI3M [nmol/cm2]")

Lm001 <-  lm(QL_4MI3M~Lesion_mm2, data=Pad3)
summary(Lm001)


ggarrange(M1, M2, M3, M4, M5, M6, 
          labels = c("A)", "B)", "C)", "D)", "E)", "F)"),
          ncol = 2, nrow = 3)

##########################
# Aliphatic  (Figure S7)

A6 <-ggplot(Col0, aes(Lesion_mm2, aliphatic)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlim(0, 160)+
  ylim(0, 0.07)+
  ggtitle("Col-0 **")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("aliphatic [nmol/cm2]")


Lm001 <-  lm(aliphatic~Lesion_mm2, data=Col0)
summary(Lm001)

A1 <- ggplot(MYB28, aes(Lesion_mm2,aliphatic)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.07)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("MYB28")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("aliphatic [nmol/cm2]")


Lm001 <-  lm(aliphatic~Lesion_mm2, data=MYB28)
summary(Lm001)

A3 <- ggplot(Cyp79, aes(Lesion_mm2,aliphatic)) +
  geom_point(aes(color=Surface), show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.07)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("cyp79b2/b3 ***")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("aliphatic [nmol/cm2]")

Lm001 <-  lm(aliphatic~Lesion_mm2, data=Cyp79)
summary(Lm001)


A2 <- ggplot(tgg12, aes(Lesion_mm2,aliphatic)) +
  geom_point(aes(color=Surface), , show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.07)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("tgg1/2")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("aliphatic [nmol/cm2]")

Lm001 <-  lm(aliphatic~Lesion_mm2, data=tgg12)
summary(Lm001)

A4 <- ggplot(m3451, aes(Lesion_mm2,aliphatic)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.07)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("myb34/51 **")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("aliphatic [nmol/cm2]")

Lm001 <-  lm(aliphatic~Lesion_mm2, data=m3451)
summary(Lm001)


A5 <- ggplot(Pad3, aes(Lesion_mm2,aliphatic)) +
  geom_point(aes(color=Surface),show.legend = FALSE)+
  geom_smooth(method = "lm", se=T, formula = y ~ x, show.legend = FALSE) +
  xlim(0, 160)+
  ylim(0, 0.07)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("pad3")+
  xlab("Lesion area at 72HPI [mm2]")+
  ylab("aliphatic [nmol/cm2]")


Lm001 <-  lm(aliphatic~Lesion_mm2, data=Pad3)
summary(Lm001)


ggarrange(A1, A2, A3, A4, A5, A6, 
          labels = c("A)", "B)", "C)", "D)", "E)", "F)"),
          ncol = 2, nrow = 3)


#### Linear mixed modeling, calculation of Least-square means, and mean-centering
#### Figure 4
# 
# library(MASS)
# library(lme4)
# library(lmerTest)
# library(car)
# library(emmeans)
# 
# GSL_Data$GenoTreat <-  as.factor(paste( GSL_Data$Geno1, GSL_Data$Treatment, sep="-"))
# 
# GSL_Data1 <- GSL_Data[complete.cases(GSL_Data),]
# GSL_Data1 <- droplevels.data.frame(GSL_Data1)
# 
# # Calculating the LS-means for each genotypes infected/uninfected
# ## Botrytis strains, the experimental tray for the infection, 
# ## and the HPLC plates are considered and random factors
# 
# lm_4M <-  lmer(QLLf_4MI3M ~ GenoTreat + (1|Iso1) + (1|Exp_Tray.x) +(1|GSL_Plates_Run), data=GSL_Data1)
# anova(lm_4M)
# summary(lm_4M)
# Emmean20 <- emmeans(lm_4M, ~ GenoTreat, lmer.df='satterthwaite')
# Emmean_4M <- as.data.frame(print(Emmean20))
# 
# colnames(Emmean_4M) <- paste(colnames(Emmean_4M), "4MI3M", sep="_")
# colnames(Emmean_4M)[1] <- c("GenoTreat")
# 
# 
# lm_3M <-  lmer(QLLf_I3M ~ GenoTreat + (1|Iso1) + (1|Exp_Tray.x) +(1|GSL_Plates_Run), data=GSL_Data1)
# anova(lm_3M)
# summary(lm_3M)
# Emmean3M <- emmeans(lm_3M, ~ GenoTreat, lmer.df='satterthwaite')
# Emmean_3M <- as.data.frame(print(Emmean3M))
# 
# colnames(Emmean_3M) <- paste(colnames(Emmean_3M), "I3M", sep="_")
# colnames(Emmean_3M)[1] <- c("GenoTreat")
# 
# 
# lm_in <-  lmer(indolic ~ GenoTreat + (1|Iso1) + (1|Exp_Tray.x) +(1|GSL_Plates_Run), data=GSL_Data1)
# anova(lm_in)
# summary(lm_in)
# Emmean_in <- emmeans(lm_in, ~ GenoTreat, lmer.df='satterthwaite')
# Emmean_in <- as.data.frame(print(Emmean_in))
# 
# colnames(Emmean_in) <- paste(colnames(Emmean_in), "Ind", sep="_")
# colnames(Emmean_in)[1] <- c("GenoTreat")
# 
# lm_al <-  lmer(aliphatic ~ GenoTreat + (1|Iso1) + (1|Exp_Tray.x) +(1|GSL_Plates_Run), data=GSL_Data1)
# anova(lm_al)
# summary(lm_al)
# Emmean_al <- emmeans(lm_al, ~ GenoTreat, lmer.df='satterthwaite')
# Emmean_al <- as.data.frame(print(Emmean_al))
# 
# colnames(Emmean_al) <- paste(colnames(Emmean_al), "Aliph", sep="_")
# colnames(Emmean_al)[1] <- c("GenoTreat")
# 
# GSL_Data <- merge(GSL_Data, Emmean_4M, by="GenoTreat")
# GSL_Data <- merge(GSL_Data, Emmean_3M, by="GenoTreat")
# GSL_Data <- merge(GSL_Data, Emmean_in, by="GenoTreat")
# GSL_Data <- merge(GSL_Data, Emmean_al, by="GenoTreat")
# 
# GSL_Data$centered_QLLf_4MI3M <- GSL_Data$QL_4MI3M - GSL_Data$emmean_4MI3M
# GSL_Data$centered_QLLf_I3M <- GSL_Data$QL_I3M - GSL_Data$emmean_I3M
# GSL_Data$centered_QLLf_aliph <- GSL_Data$aliphatic - GSL_Data$emmean_Aliph
# GSL_Data$centered_QLLf_ind <- GSL_Data$indolic - GSL_Data$emmean_Ind
# 
# GSL_Data_Treat <-  split(GSL_Data, GSL_Data$Treatment)
# GSL_Data_Treatment <- as.data.frame(GSL_Data_Treat[2])
# 
# colnames(GSL_Data_Treatment) <- colnames(GSL_Data)
# 
# GSL_Data_Treatment$Geno1 <- sub("Col0", "7_Col0", GSL_Data_Treatment$Geno1)
# GSL_Data_Treatment$Geno1 <- sub("Cyp79", "4_Cyp79", GSL_Data_Treatment$Geno1)
# GSL_Data_Treatment$Geno1 <- sub("m2829", "2_m2829", GSL_Data_Treatment$Geno1)
# GSL_Data_Treatment$Geno1 <- sub("m3451", "5_m3451", GSL_Data_Treatment$Geno1)   
# GSL_Data_Treatment$Geno1 <- sub("MYB28", "1_MYB28", GSL_Data_Treatment$Geno1)
# GSL_Data_Treatment$Geno1 <- sub("Pad3", "6_Pad3", GSL_Data_Treatment$Geno1)
# GSL_Data_Treatment$Geno1 <- sub("tgg12", "3_tgg12", GSL_Data_Treatment$Geno1)
# 
# 
# GSL_Data_Geno <- split(GSL_Data_Treatment , GSL_Data_Treatment $Geno1)
# GSL_Data_Geno  <- as.array(GSL_Data_Geno)
# 
# Sp_names <- names(GSL_Data_Geno)
# 
# emmeans_list <- list()
# 
# ############## 
# ## 4MI3M
# #############
# 
# ## Calculating the LS-mean for each genotype x surface 
# 
# GSL_centered <- as.data.frame(matrix(ncol=7))
# colnames(GSL_centered) <- c("Surface", "emmean", "SE" , "df", "lower.CL", "upper.CL","Genotype" )
# 
# anova_Surface <- as.data.frame(matrix(ncol=4))
# colnames(anova_Surface) <- c("Chisq", "Df", "Pr(>Chisq)","Surface")
# 
# 
# for (i in c(1:7)) {
#   Species_model <- lmer(centered_QLLf_4MI3M ~ Surface + (1|Isolate) + (1|Exp_Tray.x) , data=GSL_Data_Geno[[i]]) 
#   anova_tmp <- Anova(Species_model)
#   anova_tmp <- as.data.frame(anova_tmp)
#   anova_tmp$Surface <- Sp_names[[i]]
#   tmp <- emmeans(Species_model, ~Surface, lmer.df='satterthwaite')
#   tmp <- as.data.frame(print(tmp))
#   tmp$Genotype <-  rep(Sp_names[[i]],2)
#   emmeans_list[[i]] <- tmp
#   GSL_centered <- rbind(GSL_centered, tmp)
#   anova_Surface <- rbind(anova_Surface,anova_tmp)
# }
# 
# GSL_centered$Genotype <- as.factor(GSL_centered$Genotype)
# GSL_centered$plot_name <- paste(GSL_centered$Genotype, GSL_centered$Surface, sep="_")
# GSL_centered <- GSL_centered[-1,]
# 
# levels(GSL_centered$Genotype)
# 
# 
# P4m <- ggplot(GSL_centered, aes(plot_name, emmean, shape=Surface, color=Genotype)) +
#   geom_point(size=3, show.legend = FALSE)+
#   geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.8, position = position_dodge(0.1), show.legend = FALSE)+
#   theme_minimal()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   ylab("Induced 4MI3M [nmol/cm2]")
# 
# 
# 
# ############## 
# ## I3M
# #############
# 
# ## Calculating the LS-mean for each genotype x surface 
# 
# I3M_centered <- as.data.frame(matrix(ncol=7))
# colnames(I3M_centered) <- c("Surface", "emmean", "SE" , "df", "lower.CL", "upper.CL","Genotype" )
# 
# I3M_anova_Surface <- as.data.frame(matrix(ncol=4))
# colnames(I3M_anova_Surface) <- c("Chisq", "Df", "Pr(>Chisq)","Surface")
# 
# 
# for (i in c(1:7)) {
#   Species_model <- lmer(centered_QLLf_I3M ~ Surface + (1|Isolate) + (1|Exp_Tray.x) , data=GSL_Data_Geno[[i]]) 
#   anova_tmp <- Anova(Species_model)
#   anova_tmp <- as.data.frame(anova_tmp)
#   anova_tmp$Surface <- Sp_names[[i]]
#   tmp <- emmeans(Species_model, ~Surface, lmer.df='satterthwaite')
#   tmp <- as.data.frame(print(tmp))
#   tmp$Genotype <-  rep(Sp_names[[i]],2)
#   emmeans_list[[i]] <- tmp
#   I3M_centered <- rbind(I3M_centered, tmp)
#   I3M_anova_Surface  <- rbind(I3M_anova_Surface ,anova_tmp)
# }
# 
# I3M_centered$plot_name <- paste(I3M_centered$Genotype, I3M_centered$Surface, sep="_")
# I3M_centered <- I3M_centered[-1,]
# 
# p3m <- ggplot(I3M_centered, aes(plot_name, emmean, shape=Surface, color=Genotype)) +
#   geom_point(size=3, show.legend = FALSE)+
#   geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.8, position = position_dodge(0.1), show.legend = FALSE)+
#   theme_minimal()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   ylab("Induced I3M [nmol/cm2]")
# 
# 
# 
# ############## 
# ## Alphatic
# #############
# 
# ## Calculating the LS-mean for each genotype x surface 
# 
# Ali_centered <- as.data.frame(matrix(ncol=7))
# colnames(Ali_centered) <- c("Surface", "emmean", "SE" , "df", "lower.CL", "upper.CL","Genotype" )
# 
# Ali_anova_Surface <- as.data.frame(matrix(ncol=4))
# colnames(Ali_anova_Surface) <- c("Chisq", "Df", "Pr(>Chisq)","Surface")
# 
# 
# for (i in c(1:7)) {
#   Species_model <- lmer(centered_QLLf_aliph ~ Surface + (1|Isolate) + (1|Exp_Tray.x)  , data=GSL_Data_Geno[[i]]) 
#   anova_tmp <- Anova(Species_model)
#   anova_tmp <- as.data.frame(anova_tmp)
#   anova_tmp$Surface <- Sp_names[[i]]
#   tmp <- emmeans(Species_model, ~Surface, lmer.df='satterthwaite')
#   tmp <- as.data.frame(print(tmp))
#   tmp$Genotype <-  rep(Sp_names[[i]],2)
#   emmeans_list[[i]] <- tmp
#   Ali_centered <- rbind(Ali_centered, tmp)
#   Ali_anova_Surface  <- rbind(Ali_anova_Surface ,anova_tmp)
# }
# 
# Ali_centered$plot_name <- paste(Ali_centered$Genotype, Ali_centered$Surface, sep="_")
# Ali_centered <- Ali_centered[-1,]
# 
# pal <- ggplot(Ali_centered, aes(plot_name, emmean, shape=Surface, color=Genotype)) +
#   geom_point(size=3, show.legend = FALSE)+
#   geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.8, position = position_dodge(0.1), show.legend = FALSE)+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   ylab("Induced Aliphatic [nmol/cm2]")
# #theme(axis.title.x=element_blank(),
# #     axis.text.x=element_blank(),
# #    axis.ticks.x=element_blank())+
# 
# ################
# ##Camalexin
# ################
# 
# ## Calculating the LS-mean for each genotype x surface 
# 
# Cam_data_Treatment <-  split(Cam_data, Cam_data$Treatment)
# Cam_data_Treatment <- as.data.frame(Cam_data_Treatment[2])
# 
# colnames(Cam_data_Treatment) <- colnames(Cam_data)
# Cam_data_Treatment$Genotype <- sub("Col0", "7_Col0", Cam_data_Treatment$Genotype)
# Cam_data_Treatment$Genotype <- sub("Cyp79", "4_Cyp79", Cam_data_Treatment$Genotype)
# Cam_data_Treatment$Genotype <- sub("m2829", "2_m2829", Cam_data_Treatment$Genotype)
# Cam_data_Treatment$Genotype <- sub("m3451", "5_m3451", Cam_data_Treatment$Genotype)   
# Cam_data_Treatment$Genotype <- sub("MYB28", "1_MYB28", Cam_data_Treatment$Genotype)
# Cam_data_Treatment$Genotype <- sub("Pad3", "6_Pad3", Cam_data_Treatment$Genotype)
# Cam_data_Treatment$Genotype <- sub("tgg12", "3_tgg12", Cam_data_Treatment$Genotype)
# 
# Cam_data_Geno <- split(Cam_data_Treatment , Cam_data_Treatment$Genotype)
# Cam_data_Geno  <- as.array(Cam_data_Geno)
# 
# Sp_names <- names(Cam_data_Geno)
# 
# emmeans_list <- list()
# 
# Cam_centered <- as.data.frame(matrix(ncol=7))
# colnames(Cam_centered) <- c("Surface.x", "emmean", "SE" , "df", "lower.CL", "upper.CL","Genotype" )
# 
# Cam_anova_Surface <- as.data.frame(matrix(ncol=4))
# colnames(Cam_anova_Surface) <- c("Chisq", "Df", "Pr(>Chisq)","Surface.x")
# 
# 
# for (i in c(1,2,3,7)) {
#   Species_model <- lmer(Cam_LP ~ Surface.x + (1|Isolate) + (1|Exp_Tray.x) , data=Cam_data_Geno[[i]]) 
#   anova_tmp <- Anova(Species_model)
#   anova_tmp <- as.data.frame(anova_tmp)
#   anova_tmp$Surface.x <- Sp_names[[i]]
#   tmp <- emmeans(Species_model, ~Surface.x, lmer.df='satterthwaite')
#   tmp <- as.data.frame(print(tmp))
#   tmp$Genotype <-  rep(Sp_names[[i]],2)
#   emmeans_list[[i]] <- tmp
#   Cam_centered <- rbind(Cam_centered, tmp)
#   Cam_anova_Surface  <- rbind(Cam_anova_Surface, anova_tmp)
# }
# 
# Cam_centered$plot_name <- paste(Cam_centered$Genotype, Cam_centered$Surface, sep="_")
# Cam_centered <- Cam_centered[-1,]
# 
# pcam <- ggplot(Cam_centered, aes(plot_name, emmean, shape=Surface.x, color=Genotype)) +
#   geom_point(size=3, show.legend = FALSE)+
#   geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.4, position = position_dodge(0.1), show.legend = FALSE)+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   ylab("Induced Camalexin [nmol/cm2]")
# 
# library(ggpubr)
# 
# ggarrange(pcam, p3m, P4m, pal, nrow=4, ncol=1)

###########################################################################
## Chemotype modeling 
## Figure S5

Data_modeling <-  merge(GSL_Data, Cam_data, by="Leaf", all.x=T)

#model1 <- lm(LEScm.x ~ Iso1*Surface.x*Camalexin*SC_GSL*LC_GSL*QLLf_I3M*QLLf_4MI3M*QLLf_4OHI3M, data=Data_modeling)
#amodel1 <- anova(model1)

#model2 <- lm(LEScm.x ~ as.factor(Iso1)*as.factor(Surface.x*Camalexin + Iso1*Surface.x*QLLf_I3M + aliphatic + QLLf_4MI3M + QLLf_4OHI3M, data=Data_modeling)
#amodel2 <- anova(model2)
             
#model3 <- lm(LEScm.x ~ Iso1*Surface.x*Camalexin*aliphatic*indolic, data=Data_modeling)
#amodel3 <- anova(model3)
             
             
#model4 <- lm(LEScm.x ~ Iso1*Surface.x*Camalexin+aliphatic+indolic, data=Data_modeling)
#amodel4 <- anova(model4)
             
#write.table( amodel4,file="Chemotype_LM_model4.txt",sep="\t")
             
Model_pie <- read.table(file="Chemotype_LM_model_pie.txt", header=T)
             
library(plotly)
             
plot_ly(Model_pie, labels = ~Var, values = ~Perc, type = 'pie', 
                          textposition = 'inside',
                          textinfo = 'label+percent',
                          insidetextfont = list(color = '#FFFFFF'),
                          hoverinfo = 'text',
                          text = ~Var,
                          #The 'pull' attribute can also be used to create space between the sectors
                          showlegend = FALSE) %>%
               layout(title = 'Modeling on lesion area',
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
             
  
             
       