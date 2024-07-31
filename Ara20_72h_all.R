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
library(MASS)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(ggpubr)

#### This is the code for the analysis of the infection of 20 Arabidopsis thaliana genotypes with Botrytis cinerea

# Loading the raw data
G24I10 <- read.table(file="G24I10_72H.txt", header=T)

G24I10$IsoGeno <- paste(G24I10$Isolate, G24I10$Genotype, sep="_")


# Converting pixels to cm2 or mm2
G24I10 <- mutate(G24I10, 
                 Scaled.LSmm=(Lesion.Size/(8.79^2)),
                 Scaled.LScm=(Lesion.Size/(87.9^2)),
                 Leaf.Size.Cm2= Leaf.Size/(87.9^2))


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

Summary_isogens_all <- summarise(group_by(G24I10, IsoGeno),  
                                 mean_iso=mean(Scaled.LScm), 
                                 median_iso=median(Scaled.LScm), 
                                 sd_iso=sd(Scaled.LScm), 
                                 length_iso=length(Scaled.LScm))

Summary_isogens_all$IsoGen1 <- Summary_isogens_all$IsoGeno
Summary_isogens_all = separate(Summary_isogens_all, col=IsoGen1, into=c("Isolate", "Genotype"), sep="\\_")


# Plotting the distribution of lesion area 
# This distribution informs on the size threshold

ggplot(Summary_isogens_all, aes(mean_iso)) +
  geom_density() + 
  xlim(0,0.5) 

# Plotting the distribution of median of lesion area 
# This distribution informs on the median threshold

ggplot(Summary_isogens_all, aes(mean_iso)) +
  geom_density()+ 
  xlim(0,0.5) 

UD_ALL_TreatClean <- merge(G24I10, Summary_isogens_all, by="IsoGeno")

splt_Data_cleaning <- split(UD_ALL_TreatClean, UD_ALL_TreatClean$Genotype.x)
splt_Data_cleaning <- as.array(splt_Data_cleaning)


Data_cleaning <- mutate(UD_ALL_TreatClean, median_logic= median_iso > 0.2,
                        threshold_logic= Scaled.LScm <0.15, 
                        TBC_param=as.factor(paste(median_logic,threshold_logic, sep="")),
                        count_TBC=length(Scaled.LScm[TBC_param=="TRUETRUE"]))

Data_cleaning_filterOutliers <- Data_cleaning %>% filter(TBC_param!="TRUETRUE")
Failed_lesion <- Data_cleaning %>% filter(TBC_param=="TRUETRUE")



######################

## Removing columns that are not useful for this analysis
G24I10S <- Data_cleaning_filterOutliers[, -c(22:152)]
colnames(G24I10S)[c(3,4)] <- c("Isolate", "Genotype")

# Cheching the contribution to the variation in lesion area
lm1 <-  lm(Scaled.LSmm~Isolate*Genotype*Surface, data=G24I10S)
alm1 <- anova(lm1)
alm1$PercVar <- (alm1$`Sum Sq`*100)/sum(alm1$`Sum Sq`)
alm1$names1 <- rownames(alm1)
alm1$names <-   paste(rownames(alm1),round(alm1$PercVar, 2), "%", sep=" ")
pie(alm1$PercVar, labels=alm1$names1, main="percentage of variance", col=rainbow(length(alm1$names1)))


#### Analysis at the genotype level
# Calculating the emmean for Col0
lm1 <-  lmer(Scaled.LSmm~Genotype + Surface + (1|Isolate) + (1|Tray), data=G24I10S)
anova(lm1)
summary(lm1)
Emmean20 <- emmeans(lm1, ~Genotype, lmer.df='satterthwaite')
Emmean20 <- as.data.frame(print(Emmean20))

#write.table( Emmean20,file="Genotype_Emmeans.txt",sep="\t")
Cat <- read.table(file="Genotype_Emmeans_Cat.txt", header=T)

#### Mean-centering the data

# Centering by substracting the emmean of Col0
G24I10S$Centered.LSmm <- G24I10S$Scaled.LSmm-32.1


G24I10S_Geno <- split(G24I10S, G24I10S$Genotype)
G24I10S_Geno  <- as.array(G24I10S_Geno)

Sp_names <- names(G24I10S_Geno)

emmeans_list <- list()

Lesion_LS_centered <- as.data.frame(matrix(ncol=7))
colnames(Lesion_LS_centered) <- c("Surface", "emmean", "SE" , "df", "lower.CL", "upper.CL","Genotype" )

anova_Surface <- as.data.frame(matrix(ncol=4))
colnames(anova_Surface) <- c("Chisq", "Df", "Pr(>Chisq)","Surface")

# Analysis at the genotype x surface level

for (i in c(1:20)) {
  Species_model <- lmer(Centered.LSmm ~ Surface + (1|Isolate) + (1|Tray), data=G24I10S_Geno[[i]]) 
  anova_tmp <- Anova(Species_model)
  anova_tmp <- as.data.frame(anova_tmp)
  anova_tmp$Surface <- Sp_names[[i]]
  tmp <- emmeans(Species_model, ~Surface, lmer.df='satterthwaite')
  tmp <- as.data.frame(print(tmp))
  tmp$Genotype <-  rep(Sp_names[[i]],2)
  emmeans_list[[i]] <- tmp
  Lesion_LS_centered <- rbind(Lesion_LS_centered, tmp)
  anova_Surface <- rbind(anova_Surface,anova_tmp)
}

#write.table(anova_Surface ,file="lmer_anova_Surface_significance.txt",sep="\t")

Lesion_LS_centered <-    Lesion_LS_centered[-1,]

Lesion_LS_centered$SpSurf <- paste(Lesion_LS_centered$Genotype, Lesion_LS_centered$Surface, sep="_")
Lesion_LS_centered$SpSurf <- as.factor(Lesion_LS_centered$SpSurf)
Lesion_LS_centered <- merge(Lesion_LS_centered, Cat[,1:3], by="Genotype")
Lesion_LS_centered$plot_name <- paste(Lesion_LS_centered$Plot_Cat, Lesion_LS_centered$Full_Genotype,Lesion_LS_centered$Surface,  sep="_")


#Figure 3A

ggplot(Lesion_LS_centered, aes(plot_name, emmean, shape=Surface, color=Plot_Cat)) +
  geom_point(size=3, show.legend = FALSE)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.8, position = position_dodge(0.1), show.legend = FALSE)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Mean-centered Lesion Area")


########################
### Analysis of variance
########################

Spl_G24I10 <- split(G24I10S, G24I10S$Genotype)
list2env(Spl_G24I10, envir=.GlobalEnv)

names(Spl_G24I10)

# Creating 19 datasets, with 1 mutant and Col-0

Spl_G24I10_pairwise <- as.list(array(dim=19))

ilist <- c(1:6, 8:20)

for (j in c(1:19)){
    i = ilist[j]
    names(Spl_G24I10_pairwise)[j] <- paste(names(Spl_G24I10[i]), names(Spl_G24I10[7]), sep="_")
    temp1 <- as.data.frame(Spl_G24I10[i])
    colnames(temp1) <- colnames(G24I10S)
    Col0_data <- as.data.frame(Spl_G24I10[7])
    colnames(Col0_data ) <- colnames(G24I10S)
    Spl_G24I10_pairwise[[j]] <- rbind(temp1, Col0_data)
    rm(temp1)
  }  

### Calculating pairwise ANOVA

Anova_Results_pairwise <- as.data.frame(matrix(ncol=32, nrow=19))
for (j in c(1:19))  {
  Anova_Results_pairwise[j,1] <-  names(Spl_G24I10_pairwise[j])
  temp1 <- droplevels.data.frame( Spl_G24I10_pairwise[[j]])
  lmt <- lm(formula=Scaled.LScm ~ Genotype*Surface*Isolate, data=temp1)
  lmt1 <- summary(lmt)
  aovlmt <-  anova(lmt)
  aovlmt <- broom::tidy(aovlmt)
  Anova_Results_pairwise[j,2] <- aovlmt$term[1]
  Anova_Results_pairwise[j,3] <- aovlmt$sumsq[1]
  Anova_Results_pairwise[j,4] <- aovlmt$meansq[1]
  Anova_Results_pairwise[j,5] <- aovlmt$p.value[1]
  Anova_Results_pairwise[j,6] <- aovlmt$term[2]
  Anova_Results_pairwise[j,7] <- aovlmt$sumsq[2]
  Anova_Results_pairwise[j,8] <- aovlmt$meansq[2]
  Anova_Results_pairwise[j,9] <- aovlmt$p.value[2]
  Anova_Results_pairwise[j,10] <- aovlmt$term[3]
  Anova_Results_pairwise[j,11] <- aovlmt$sumsq[3]
  Anova_Results_pairwise[j,12] <- aovlmt$meansq[3]
  Anova_Results_pairwise[j,13] <- aovlmt$p.value[3]
  Anova_Results_pairwise[j,14] <- aovlmt$term[4]
  Anova_Results_pairwise[j,15] <- aovlmt$sumsq[4]
  Anova_Results_pairwise[j,16] <- aovlmt$meansq[4]
  Anova_Results_pairwise[j,17] <- aovlmt$p.value[4]
  Anova_Results_pairwise[j,18] <- aovlmt$term[5]
  Anova_Results_pairwise[j,19] <- aovlmt$sumsq[5]
  Anova_Results_pairwise[j,20] <- aovlmt$meansq[5]
  Anova_Results_pairwise[j,21] <- aovlmt$p.value[5]
  Anova_Results_pairwise[j,22] <- aovlmt$term[6]
  Anova_Results_pairwise[j,23] <- aovlmt$sumsq[6]
  Anova_Results_pairwise[j,24] <- aovlmt$meansq[6]
  Anova_Results_pairwise[j,25] <- aovlmt$p.value[6]
  Anova_Results_pairwise[j,26] <- aovlmt$term[7]
  Anova_Results_pairwise[j,27] <- aovlmt$sumsq[7]
  Anova_Results_pairwise[j,28] <- aovlmt$meansq[7]
  Anova_Results_pairwise[j,29] <- aovlmt$p.value[7]
  Anova_Results_pairwise[j,30] <- aovlmt$term[8]
  Anova_Results_pairwise[j,31] <- aovlmt$sumsq[8]
  Anova_Results_pairwise[j,32] <- aovlmt$meansq[8]
  rm(temp1, lmt, lmt1, aovlmt)
}

colnames(Anova_Results_pairwise) <- c("GenotypePair", "term1", "SSqr1", "meansq1", "pval1", "term2", "SSqr2", "meansq2", "pval2","term3", "SSqr3", "meansq3", "pval3", "term4", "SSqr4", "meansq4", "pval4", "term5", "SSqr5", "meansq5", "pval5", "term6", "SSqr6", "meansq6", "pval6","term7", "SSqr7", "meansq7", "pval7", "term8", "SSqr8", "meansq8")

##### Analysis of variance od the surface within genotype

Anova_Results1 <- as.data.frame(matrix(ncol=12, nrow=20))
#Genotype
#Pval

for (i in c(1:20)) {
  Anova_Results1[i, 1]  <- names(Spl_G24I10)[i]
  
  temp <- droplevels.data.frame( as.data.frame(Spl_G24I10[i]))
  colnames(temp) <- colnames(G24I10S)
  Aov_sum <- as.data.frame(t(unlist(summary(aov(Scaled.LSmm ~ Surface, data=temp)))))
  lmt <- lm(formula=Scaled.LScm ~ Surface, data=temp)
  Anova_Results1[i, 9] <- confint(lmt, level=0.90)[2]
  Anova_Results1[i, 10] <- confint(lmt,  level=0.90)[4]
  lmt <- summary(lmt)

  
  Anova_Results1[i, 2]  <- Aov_sum$`Pr(>F)1`
  Anova_Results1[i, 3]  <- Aov_sum$`Sum Sq1`
  Anova_Results1[i, 4]  <- Aov_sum$`Sum Sq2`
  Anova_Results1[i, 5] <- lmt$coefficients[1]
  Anova_Results1[i, 6] <- lmt$coefficients[2]
  Anova_Results1[i, 7] <- lmt$coefficients[3]
  Anova_Results1[i, 8] <- lmt$coefficients[4]
  Anova_Results1[i, 11] <- lmt$r.squared
  Anova_Results1[i, 12] <- lmt$adj.r.squared
  
  rm(temp, lmt, Aov_sum)
}

colnames(Anova_Results1) <- c("Genotype", "Pval", "SSquare_surface", "Ssquare_res", "intercept", "slope", "Int_error", "slope_error", "CI95min", "CI95max","Rsq", "adj_Rsq")
Anova_Results1$fdr <- p.adjust(Anova_Results1$Pval, method = "fdr", n = length(Anova_Results1$Pval))


Anova_Results1$Sign <- NA

for (i in c(1:dim(Anova_Results1)[1])){
  if (Anova_Results1$Pval[i]<0.05) {
    Anova_Results1$Sign[i] <- c("Signif")
  } else {
    Anova_Results1$Sign[i] <- c("UnSignif")
  }
}

write.table(Anova_Results1, file="Anova_Results_singleMut.txt",sep="\t")


########################################################
# Analysis across time
# Compiling and filtering data at 48 and 96h for the cleaned dataset at 72h
########################################################

## Identifying every leaf in each tray across time
## to remove the lesions that were cleaned at 72h 

D96 <- read.table(file="G24I10_96h.txt", header=T)
D96 <-  separate(D96, Image, into=c( "Tray", "IMG", "N_IMG", "trash"), sep="_")
D96 <-  D96[, -c(6,8)]
D96$Key <- paste(D96$Tray, D96$Isolate, D96$Genotype, sep="_")
D96S <- D96[, -c(6:148)]
D96S <- mutate(D96S, Scaled.LSmm=(Lesion.Size/(8.892^2)),Scaled.LScm=(Lesion.Size/(88.92^2)) )
colnames(D96S)[c(1:12, 14:15)] <- paste("D96", colnames(D96S)[c(1:12, 14:15)], sep="_")



D48 <- read.table(file="G24I10_UD_48h.txt", header=T)
D48 <-  separate(D48, Image, into=c( "Tray", "IMG", "N_IMG", "trash"), sep="_")
D48 <-  D48[, -c(6,8)]
D48$Key <- paste(D48$Tray, D48$Isolate, D48$Genotype, sep="_")
D48S <- D48[, -c(6:148)]
D48S <- mutate(D48S, Scaled.LSmm=(Lesion.Size/(8.892^2)),Scaled.LScm=(Lesion.Size/(88.92^2)) )
colnames(D48S)[c(1:12, 14:15)] <- paste("D48", colnames(D48S)[c(1:12, 14:15)], sep="_")


G24I10S <- Data_cleaning_filterOutliers[, -c(22:152)]
colnames(G24I10S)[c(3,4)] <- c("Isolate", "Genotype")
G24I10S$Key <- paste(G24I10S$Tray, G24I10S$IsoGeno, G24I10S$Surf, sep="_")


G24I10T3 <- merge(G24I10S, D48S, by="Key", all.x=T)
G24I10T3 <- merge(G24I10T3, D96S, by="Key", all.x=T)

rm(D48, D48S, D96, D96S)
rm(Data_cleaning, Data_cleaning_filterOutliers, Failed_lesion, G24I10, splt_Data_cleaning, Summary_isogens_all, UD_ALL_TreatClean)

###################
### Merged Data 3Times and calculate durface effect and relative growth rate

G24I10T3S <- G24I10T3[, c(1,2, 4:8, 10,11, 24:31, 49:56, 63:70)]
G24I10T3SS <- G24I10T3S[, c(1:9, 17, 25, 33)]
colnames(G24I10T3SS)[10:12] <- c("D72", "D48", "D96")
rm(G24I10T3S, G24I10T3)

### Calculating the differential in lesion area across leaf surfaces (dSurf)
# The linear growth rate is also calculated (linRGR)

G24I10T3S_splt <- split(G24I10T3SS, G24I10T3SS$Surface)

tmp1 <- as.data.frame(G24I10T3S_splt[1])
tmp1$Key <- paste(tmp1$Abaxial.Tray, tmp1$Abaxial.IsoGeno)
tmp2 <- as.data.frame(G24I10T3S_splt[2])
tmp2$Key <- paste(tmp2$Adaxial.Tray, tmp2$Adaxial.IsoGeno)

G24I10T3Surf <- merge(tmp1, tmp2, by="Key")
rm(tmp1, tmp2)

G24I10T3Surf <- mutate(G24I10T3Surf,
                       dSurf48 = Abaxial.D48-Adaxial.D48,
                       dSurf72 = Abaxial.D72-Adaxial.D72,
                       dSurf96 = Abaxial.D96-Adaxial.D96,
                       linRGR0_ab= (Abaxial.D48)/6,
                       linRGR0_ad= (Adaxial.D48)/6,
                       linRGR1_ab= (Abaxial.D72-Abaxial.D48)/24,
                       linRGR1_ad= (Adaxial.D72-Adaxial.D48)/24,
                       linRGR2_ab= (Abaxial.D96-Abaxial.D72)/24,
                       linRGR2_ad= (Adaxial.D96-Adaxial.D72)/24)

G24I10T3_dSurf_long <- pivot_longer(G24I10T3Surf[, c(1, 3:5, 26:28)], dSurf48:dSurf96, names_to="Time", values_to="dSurf")

G24I10T3_dSurf_long$IsoGenTime <- paste(G24I10T3_dSurf_long$Abaxial.IsoGeno, G24I10T3_dSurf_long$Time, sep="_")
G24I10T3_dSurf_sum <- summarise(group_by(G24I10T3_dSurf_long, IsoGenTime),
                                mean_dSurf=mean(dSurf,na.rm = T), 
                                median_dSurf=median(dSurf,na.rm = T), 
                                sd_dSurf=sd(dSurf,na.rm = T), 
                                length_dSurf=length(dSurf))

G24I10T3_dSurf_sum$IsoGenTime <- as.factor(gsub('_dSurf', '_D', G24I10T3_dSurf_sum$IsoGenTime, fixed=F))


G24I10T3_dSurf_long$GenTime <- paste(G24I10T3_dSurf_long$Abaxial.Genotype, G24I10T3_dSurf_long$Time, sep="_")
G24I10T3_dSurf_GT <- summarise(group_by(G24I10T3_dSurf_long, GenTime),
                               mean_dSurf=mean(dSurf,na.rm = T), 
                               median_dSurf=median(dSurf,na.rm = T), 
                               sd_dSurf=sd(dSurf,na.rm = T), 
                               length_dSurf=length(dSurf))
G24I10T3_dSurf_GT$GenTime <- as.factor(gsub('_dSurf', '_D', G24I10T3_dSurf_GT$GenTime, fixed=F))


G24I10T3S_long <- pivot_longer(G24I10T3SS, D72:D96, names_to="Time", values_to="Lsmm")
G24I10T3S_long$group <- paste(G24I10T3S_long$IsoGeno, G24I10T3S_long$Surf, G24I10T3S_long$Time, sep="_")
G24I10T3S_long$IsoGenTime <- paste(G24I10T3S_long$IsoGeno, G24I10T3S_long$Time, sep="_")
G24I10T3S_long$GenTime <- paste(G24I10T3S_long$Genotype, G24I10T3S_long$Time, sep="_")

G24I10_IGT_sum <- summarise(group_by(G24I10T3S_long, IsoGenTime),
                            mean_lesion=mean(Lsmm,na.rm = T), 
                            median_lesion=median(Lsmm,na.rm = T), 
                            sd_lesion=sd(Lsmm,na.rm = T), 
                            length_lesion=length(Lsmm))

G24I10_GT_sum <- summarise(group_by(G24I10T3S_long, GenTime),
                           mean_lesion=mean(Lsmm,na.rm = T), 
                           median_lesion=median(Lsmm,na.rm = T), 
                           sd_lesion=sd(Lsmm,na.rm = T), 
                           length_lesion=length(Lsmm))


Data_Plot_Lesion <- merge(G24I10_IGT_sum,G24I10T3_dSurf_sum, by="IsoGenTime")
Data_Plot_Lesion$IsoGenTime1 <- Data_Plot_Lesion$IsoGenTime
Data_Plot_Lesion <- separate(Data_Plot_Lesion, IsoGenTime, into=c("Isolate", "Genotype", "Time"), sep="_")
Data_Plot_Lesion$Isolate <- as.factor(Data_Plot_Lesion$Isolate)
Data_Plot_Lesion$Genotype <- as.factor(Data_Plot_Lesion$Genotype)
Data_Plot_Lesion$Time <- as.factor(Data_Plot_Lesion$Time)


Data_Plot_Lesion_GT <- merge(G24I10_GT_sum,G24I10T3_dSurf_GT, by="GenTime")
Data_Plot_Lesion_GT$GenTime1 <- Data_Plot_Lesion_GT$GenTime
Data_Plot_Lesion_GT<- separate(Data_Plot_Lesion_GT, GenTime, into=c("Genotype", "Time"), sep="_")
Data_Plot_Lesion_GT$Genotype <- as.factor(Data_Plot_Lesion_GT$Genotype)
Data_Plot_Lesion_GT$Time <- as.factor(Data_Plot_Lesion_GT$Time)

Data_Plot_Lesion_GTS <- Data_Plot_Lesion_GT[, c(1,2,11, 3,7)]

Lesion_Data <- pivot_wider(Data_Plot_Lesion_GTS, id_cols = Genotype, values_from=mean_lesion:mean_dSurf, names_from =Time)

#write.table( Lesion_Data,file="Data_Plot_Lesion_GT.txt",sep="\t")


###############
## Dotted line plot
###############

palette1 <- c("#a6dba0", "#ae017e", "#315717", "#df65b0",  "#88419d" , "#d16eff",
              "#050505","#ba6aae", "#FDBE85", "#737373",  "#02adbd", "#9ECAE1", 
              "#762a83", "#5aae61", "#6BAED6", "#3182BD", "#dd3497",
              "#00441b", "#4d9221", "#FD8D3C", "#F8766D", "#00BA38", "#619CFF")

Data_Plot_Lesion_GT <- read.table(file="Data_Plot_Lesion_GT.txt", header=T)
Data_Plot_Lesion_GT1 <- Data_Plot_Lesion_GT[, c(1:4)]
colnames(Data_Plot_Lesion_GT1) <- c("Genotype", "D48", "D72", "D96")
Data_Plot_Lesion_GT2 <- Data_Plot_Lesion_GT[, c(1, 5:7)]
colnames(Data_Plot_Lesion_GT2) <- c("Genotype", "D48", "D72", "D96")

tmp1 <- pivot_longer(Data_Plot_Lesion_GT1, D48:D96, names_to="Time_mean", values_to="Mean_lesion")
tmp1$key <- paste(tmp1$Genotype, tmp1$Time_mean, sep="_")
tmp2 <- pivot_longer(Data_Plot_Lesion_GT2, D48:D96, names_to="Time_surf", values_to="Surface_Diff")
tmp2$key <- paste(tmp2$Genotype, tmp2$Time_surf, sep="_")

Data_Plot_Lesion_GTL <- merge(tmp1, tmp2, by="key")
rm(tmp1, tmp2)


library(directlabels)

palette1 <- c("#a6dba0", "#ae017e", "#315717", "#df65b0",  "#88419d" , "#d16eff",
              "#050505","#ba6aae", "#FDBE85", "#737373",  "#02adbd", "#9ECAE1", 
              "#762a83", "#5aae61", "#6BAED6", "#3182BD", "#dd3497",
              "#00441b", "#4d9221", "#FD8D3C", "#F8766D", "#00BA38", "#619CFF")

PaletteN <- c("#800000", "#3b2e7d", "#a10b0b","#682085", "#ff0066",   
              "#bb21cc",  "#211f20", "#f7488e", "#45105c","#666364",
              "#01657d", "#14caf7",  "#2a6ed4", "#ff6600", "#2c89a0",
              "#52b8d1",  "#52329c", "#bf5915", "#db7530",  "#2a4fb5")

PaletteP <- c("#f76454", "#479c44", "#279ff5")



levels(Data_Plot_Lesion_GTL$Genotype.x) <- c("anac055", "cyp79b2/b3", "coi1-16","cyp71A12","cyp81D8",  "cyp82C2", "col-0" , "fox5","gGP1",   
                                             "Ler", "myb28/29", "myb29", "myb34/51", "npr1", "MYB28", "MYB29", "pad3", "pad4",
                                             "tga3-2","tgg1/2")


#Figure 3B
ggplot(data=Data_Plot_Lesion_GTL, aes(Mean_lesion, Surface_Diff)) +
  geom_point(aes(color=Time_mean), show.legend = FALSE) +
  scale_fill_manual(values=PaletteP)+
  ylab("Surface differential in Lesion area [mm2]")+
  xlab("Lesion area [mm2]")+
  theme_bw()


q<- ggplot(data=Data_Plot_Lesion_GTL, aes(Mean_lesion, Surface_Diff)) + 
  geom_line(aes( group= Genotype.x, color=Genotype.x), linetype=2)+
  scale_color_manual(values=PaletteN)+
  theme_bw()+
  ylab("Surface differential in Lesion area [mm2]")+
  xlab("Lesion area [mm2]")

direct.label(q,"last.points")

