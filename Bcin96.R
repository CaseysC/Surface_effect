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
library(lme4)
library(lmerTest)
library(car)
library(emmeans)


source('geom_split_violin.R')

#### This is the code for the analysis of 96 Botrytis strains infected on Col-0

############################
#### Loading the raw data
############################

LSData_48h <-  read.table(file="Lesion_C0l0_48h.txt", header=T)
LSData_48h <- mutate(LSData_48h, Scaled.LSmm=(Lesion.Size/(8.91^2)),Scaled.LScm=(Lesion.Size/(89.1^2)) )
levels(LSData_48h$Genotype) <- c("Abaxial", "Adaxial")
LSData_48h$IsoGen <- as.factor(paste(LSData_48h$Isolate, LSData_48h$Genotype, sep="-"))

Data_72h <- read.table(file="Col0_UpDown72h.txt", header=T)
Data_72h <- mutate(Data_72h, Scaled.LSmm=(Lesion.Size/(8.91^2)),Scaled.LScm=(Lesion.Size/(89.1^2)) )
levels(Data_72h$Genotype) <- c("Abaxial", "Adaxial")
Data_72h$IsoGen <- as.factor(paste(Data_72h$Isolate, Data_72h$Genotype, sep="-"))

Data_96h <- read.table(file="UpDown_Exp1_Results96.txt", header=T)
Data_96h <- mutate(Data_96h, Scaled.LSmm=(Lesion.Size/(8.91^2)),Scaled.LScm=(Lesion.Size/(89.1^2)) )
levels(Data_96h$Genotype) <- c("Abaxial", "Adaxial")
Data_96h$IsoGen <- as.factor(paste(Data_96h$Isolate, Data_96h$Genotype, sep="-"))


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
#######################################

Summary_isogens_96 <- summarise(group_by(Data_96h, IsoGen),  
                                mean_iso=mean(Scaled.LScm), 
                                median_iso=median(Scaled.LScm), 
                                sd_iso=sd(Scaled.LScm), 
                                length_iso=length(Scaled.LScm))

Summary_isogens_96$IsoGen1 <- Summary_isogens_96$IsoGen
Summary_isogens_96 = separate(Summary_isogens_96, col=IsoGen1, into=c("Sname", "Isolate", "Surface"), sep="\\-")


ggplot(Summary_isogens_96, aes(mean_iso)) +
  geom_density() + 
  xlim(0,0.5) 

ggplot(Summary_isogens_96, aes(median_iso)) +
  geom_density()+ 
  xlim(0,0.5) 

UD_TreatClean_96 <- merge(Data_96h, Summary_isogens_96, by="IsoGen")

Summary_isogens_96$IsoGen <- paste(Summary_isogens_96$Isolate, Summary_isogens_96$Surface, sep="-")
UD_TreatClean_72 <- merge(Data_72h, Summary_isogens_96, by="IsoGen")


Data_cleaning <- mutate(UD_TreatClean_96, 
                        median_logic= median_iso > 0.2,
                        threshold_logic= Scaled.LScm <0.15, 
                        TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                        count_TBC=length(Scaled.LScm[TBC_param=="TRUETRUE"]))

Data_cleaning72 <- mutate(UD_TreatClean_72, 
                          median_logic= median_iso > 0.2,
                          threshold_logic= Scaled.LScm <0.15, 
                          TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                          count_TBC=length(Scaled.LScm[TBC_param=="TRUETRUE"]))

Data_72_filterOutliers <- Data_cleaning72 %>% filter(TBC_param!="TRUETRUE")
Failed_lesion_72 <- Data_cleaning72 %>% filter(TBC_param=="TRUETRUE")

Data_96_filterOutliers <- Data_cleaning %>% filter(TBC_param!="TRUETRUE")
Failed_lesion_96 <- Data_cleaning %>% filter(TBC_param=="TRUETRUE")

##############################################################


LSData_48hm<- lmer(Scaled.LSmm ~ Isolate + Genotype + (1|Image) + (1|Plant_tray/Plant), data=LSData_48h)
UpDown48_model_Anova<- anova(LSData_48hm)
H48_Isolate_lsm <- emmeans(LSData_48hm, ~Isolate, lmer.df='satterthwaite')
H48_Isolate_lsm <- as.data.frame(print(H48_Isolate_lsm))
H48_Isolate_lsm <- separate(H48_Isolate_lsm, Isolate, into=c("Sname", "Isolate"), sep="-")
#write.table(H48_Isolate_lsm, file= "H48_lsmean.txt")

H48_Geno_lsm <- emmeans(LSData_48hm, ~Genotype, lmer.df='satterthwaite')
H48_Geno_lsm <- as.data.frame(print(H48_Geno_lsm))

Sgen48 <- summarise(group_by(LSData_48h, Genotype) ,
                    mean_geno=mean(Scaled.LSmm), 
                    median_geno=median(Scaled.LSmm), 
                    sd_geno=sd(Scaled.LSmm), 
                    length_genoo=length(Scaled.LSmm))

IG_LSData_48hm<- lmer(Scaled.LSmm ~ IsoGen + (1|Image) + (1|Plant_tray/Plant), data=LSData_48h)
IG48_model_Anova<- anova(IG_LSData_48hm)
H48_IsoGen_lsm <- emmeans(IG_LSData_48hm, ~IsoGen, lmer.df='satterthwaite')
H48_IsoGen_lsm <- as.data.frame(print(H48_IsoGen_lsm))
#write.table(H48_IsoGen_lsm, file= "IsoGen_H48_lsmean.txt")


H48_IsoGen_lsm$IsoGen1 <- H48_IsoGen_lsm$IsoGen
H48_IsoGen_lsm <- separate(H48_IsoGen_lsm, IsoGen, into=c("Sname", "Isolate", "Surface"), sep="-")
H48_IsoGen_lsm_col <- pivot_wider(H48_IsoGen_lsm, id_cols = Isolate, values_from=emmean:upper.CL, names_from =Surface)
H48_IsoGen_lsm_col <- mutate(H48_IsoGen_lsm_col, Lsm48_diff = emmean_Abaxial-emmean_Adaxial)

All_H48_Iso_lsm <- merge(H48_Isolate_lsm,H48_IsoGen_lsm_col, By="Isolate")

LmData_48hm <- lm(Scaled.LSmm ~ Isolate * Genotype + Image + Plant_tray/Plant, data=LSData_48h)
Lm48_model_Anova<- anova(LmData_48hm)

aov48lmtable <- broom::tidy(Lm48_model_Anova)
aov48lmtable <- as.data.frame(aov48lmtable)
aov48lmtable <-  mutate(aov48lmtable,
                        heritabilityPerc=sumsq*100/sum(aov48lmtable$sumsq),
                        heritability=sumsq/sum(aov48lmtable$sumsq))
aov48lmtable$Time <- rep(48, 7)
aov48lmtable$Cat <- c("a", "a", "b", "b", "a", "b", "b")

######
#72H Data
######


LSData_72hm<- lmer(Scaled.LSmm ~ Isolate.x + Genotype + (1|Image) + (1|Plant_tray/Plant), data=Data_72_filterOutliers)
UpDown72_model_Anova<- anova(LSData_72hm)
H72_Isolate_lsm <- emmeans(LSData_72hm, ~Isolate.x, lmer.df='satterthwaite')
H72_Isolate_lsm <- as.data.frame(print(H72_Isolate_lsm))
#write.table(H72_Isolate_lsm, file= "H72_lsmean.txt")
colnames(H72_Isolate_lsm)[1] <- c("Isolate")

Sgen72 <- summarise(group_by(Data_72_filterOutliers, Genotype) ,
                    mean_geno=mean(Scaled.LSmm), 
                    median_geno=median(Scaled.LSmm), 
                    sd_geno=sd(Scaled.LSmm), 
                    length_genoo=length(Scaled.LSmm))

IG_LSData_72hm<- lmer(Scaled.LSmm ~ IsoGen + (1|Image) + (1|Plant_tray/Plant), data=Data_72_filterOutliers)
IG72_model_Anova<- anova(IG_LSData_72hm)
H72_IsoGen_lsm <- emmeans(IG_LSData_72hm, ~IsoGen, lmer.df='satterthwaite')
H72_IsoGen_lsm <- as.data.frame(print(H72_IsoGen_lsm))
#write.table(H72_IsoGen_lsm, file= "IsoGen_H72_lsmean.txt")
colnames(H72_IsoGen_lsm)

H72_IsoGen_lsm$IsoGen1 <- H72_IsoGen_lsm$IsoGen
H72_IsoGen_lsm <- separate(H72_IsoGen_lsm, IsoGen, into=c("Isolate", "Surface"), sep="-")
H72_IsoGen_lsm_col <- pivot_wider(H72_IsoGen_lsm, id_cols = Isolate, values_from=emmean:upper.CL, names_from =Surface)
H72_IsoGen_lsm_col <- mutate(H72_IsoGen_lsm_col, Lsm72_diff = emmean_Abaxial-emmean_Adaxial)

All_H72_Iso_lsm <- merge(H72_Isolate_lsm,H72_IsoGen_lsm_col, By="Isolate")

LmData_72hm <- lm(Scaled.LSmm ~ Isolate.x * Genotype + Image + Plant_tray/Plant, data=Data_72_filterOutliers)
Lm72_model_Anova<- anova(LmData_72hm)

aov72lmtable <- broom::tidy(Lm72_model_Anova)
aov72lmtable <- as.data.frame(aov72lmtable)
aov72lmtable <-  mutate(aov72lmtable,
                        heritabilityPerc=sumsq*100/sum(aov72lmtable$sumsq),
                        heritability=sumsq/sum(aov72lmtable$sumsq))
aov72lmtable$Time <- rep(72, 7)
aov72lmtable$Cat <- c("a", "a", "b", "b", "a", "b", "b")

######
# 96H Data
######

LSData_96hm<- lmer(Scaled.LSmm ~ Isolate.x + Genotype + (1|Image) + (1|Plant_tray/Plant), data=Data_96_filterOutliers)
UpDown96_model_Anova<- anova(LSData_96hm)
H96_Isolate_lsm <- emmeans(LSData_96hm, ~Isolate.x, lmer.df='satterthwaite')
H96_Isolate_lsm <- as.data.frame(print(H96_Isolate_lsm))
H96_Isolate_lsm <- separate(H96_Isolate_lsm, Isolate.x, into=c("Sname", "Isolate"), sep="-")
#write.table(H96_Isolate_lsm, file= "H96_lsmean.txt")

IG_LSData_96hm<- lmer(Scaled.LSmm ~ IsoGen + (1|Image) + (1|Plant_tray/Plant), data=Data_96_filterOutliers)
IG96_model_Anova<- anova(IG_LSData_96hm)
H96_IsoGen_lsm <- emmeans(IG_LSData_96hm, ~IsoGen, lmer.df='satterthwaite')
H96_IsoGen_lsm <- as.data.frame(print(H96_IsoGen_lsm))
#write.table(H96_IsoGen_lsm, file= "IsoGen_H96_lsmean.txt")

Sgen96 <- summarise(group_by(Data_96_filterOutliers, Genotype) ,
                    mean_geno=mean(Scaled.LSmm), 
                    median_geno=median(Scaled.LSmm), 
                    sd_geno=sd(Scaled.LSmm), 
                    length_genoo=length(Scaled.LSmm))

H96_IsoGen_lsm$IsoGen1 <- H96_IsoGen_lsm$IsoGen
H96_IsoGen_lsm <- separate(H96_IsoGen_lsm, IsoGen, into=c("Sname", "Isolate", "Surface"), sep="-")
H96_IsoGen_lsm_col <- pivot_wider(H96_IsoGen_lsm, id_cols = Isolate, values_from=emmean:upper.CL, names_from =Surface)
H96_IsoGen_lsm_col <- mutate(H96_IsoGen_lsm_col, Lsm48_diff = emmean_Abaxial-emmean_Adaxial)

All_H96_Iso_lsm <- merge(H96_Isolate_lsm,H96_IsoGen_lsm_col, By="Isolate")

LmData_96hm <- lm(Scaled.LSmm ~ Isolate.x * Genotype + Image + Plant_tray/Plant, data=Data_96_filterOutliers)
Lm96_model_Anova<- anova(LmData_96hm)
aov96lmtable <- broom::tidy(Lm96_model_Anova)
aov96lmtable <- as.data.frame(aov96lmtable)
aov96lmtable <-  mutate(aov96lmtable,
                        heritabilityPerc=sumsq*100/sum(aov96lmtable$sumsq),
                        heritability=sumsq/sum(aov96lmtable$sumsq))
aov96lmtable$Time <- rep(96, 7)
aov96lmtable$Cat <- c("a", "a", "b", "b", "a", "b", "b")

############

Heritability <- rbind(aov48lmtable, aov72lmtable, aov96lmtable)
Heritability$term <- as.factor(Heritability$term)
levels(Heritability$term)[c(5,6)] <- c("Isolate", "Isolate:Genotype" )
Heritability$Cat <- as.factor(Heritability$Cat)

library(viridis)

PalCol <- viridis_pal(alpha = 0.8, begin = 0.2, end = 0.9, direction = 1,
                      option = "A")(7)


## Figure 5B
ggplot(data=Heritability, aes(Time, heritabilityPerc, color=term)) +
  geom_line(aes(linetype=Cat))+
  geom_point() +
  scale_color_manual(values=PalCol)+
  ylab("Percentage of the variance in lesion area")+
  xlab("Hours post inoculation")+
  theme_bw()



########
# Merging all
########

All_H48_Iso_lsm$time <- as.numeric(rep("48", 99))
colnames(All_H48_Iso_lsm)[18] <- c("Lsm_diff")

All_H72_Iso_lsm$time <- as.numeric(rep("72", 99))
colnames(All_H72_Iso_lsm)[17] <- c("Lsm_diff")

All_H96_Iso_lsm$time <- as.numeric(rep("96", 99))
colnames(All_H96_Iso_lsm)[18] <- c("Lsm_diff")


All_lsm_497296 <- rbind(All_H48_Iso_lsm[,-2], All_H72_Iso_lsm, All_H96_Iso_lsm[,-2])

#removing the controls

All_lsm_497296 <- All_lsm_497296[which(All_lsm_497296$Isolate!='ctrl'),]
All_lsm_497296 <- All_lsm_497296[which(All_lsm_497296$Isolate!='Ctrl'),]
All_lsm_497296$Time <- as.factor(All_lsm_497296$time)

All_lsm_497296$chosen <- rep("Shadow", 291)
rownames(All_lsm_497296) <- c(1:291)
#All_lsm_497296_1 <-   All_lsm_497296[order(All_lsm_497296$Isolate),]

All_lsm_497296$chosen[c(44,47,51,52,141,144,148,149, 238, 241, 245, 246)] <- c("ab")
All_lsm_497296$chosen <- as.factor(All_lsm_497296$chosen)

All_lsm_497296 <- All_lsm_497296[which(All_lsm_497296$Isolate!='01_05_17'),]


### Iso * Surface LS-means Data
## Figure 5C
ggplot(data=All_lsm_497296, aes(emmean, Lsm_diff)) +
  facet_grid(.~chosen)+
  #facet_wrap(Isolate~.)+
  #geom_line(aes(group= Isolate))+
  geom_line(aes( group= Isolate,linetype=chosen))+
  scale_linetype_manual(values=c("solid", "dotted"))+
  geom_point(aes(color=Time)) +
  ylab("Surface differential in Lesion area [mm2]")+
  xlab("Lesion area [mm2]")+
  theme_bw()


### Emmean based half violin plots
### Figure 5A

colnames(All_lsm_497296)
S1 <- unique(All_lsm_497296[, c(1, 19, 7,8)])
S1L <- pivot_longer(S1, emmean_Abaxial:emmean_Adaxial, names_to="Surface", values_to="emmean")



## Figure 5A
ggplot(data=S1L, aes(x=Surface, y=emmean)) + 
  facet_grid(.~Time)+
  #geom_point()+
  geom_split_violin(fill="grey")+
  geom_boxplot(width=0.05)+ 
  ylab("Lesion area [mm2]")+
  theme_bw()

