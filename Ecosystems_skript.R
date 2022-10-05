#####################################################################################################################
#               Variability in tree-ring width and NDVI responses to climate at a landscape level                   #
#####################################################################################################################

#Jiří Mašek a*, Jan Tumajer a, Jelena Lange a, Ryszard Kaczka a, Petr Fišer a, Václav Treml a
#a Department of Physical Geography and Geoecology, Faculty of Science, Charles University, Albertov 6, 128 43 Prague, Czech Republic
#* Corresponding author: jiri.masek@natur.cuni.cz (Jiří Mašek)

### Necessary R packages
library(dplR)
library(treeclim)
library(tidyr)
library(devtools)
library(dendRolAB)

# to instal the last package, run this comand:
#devtools::install_github("AllanBuras/dendRolAB")

### set working directory
setwd("C:/Users/jirka/Desktop/Data/")

###############################################################
#                 1. Detrending of NDVI data                  #
###############################################################

NDVI_orig <- read.table("NDVI_orig.txt", check.names=FALSE, row.names = 1, , header=TRUE, sep="\t", na.strings="NA", dec=",")

NDVI_det<- data.frame(matrix(ncol = ncol(NDVI_orig), nrow = nrow(NDVI_orig))); colnames(NDVI_det)<- colnames(NDVI_orig); rownames(NDVI_det)<- rownames(NDVI_orig)

for (a in c(1:ncol(NDVI_orig))) {
  
  model<- lm(NDVI_orig[,a] ~as.numeric(rownames(NDVI_orig)))
  
  NDVI_det[,a] <- model$residuals+mean(NDVI_orig[,a])
}

NDVI<- NDVI_det[rownames(NDVI_det)<2018,]


###############################################################
#          2. Calculating of TRI chronologies                #
###############################################################

List_loc<- read.table("List_loc.txt", header=TRUE, sep="\t", na.strings="NA", dec=",", strip.white=TRUE, encoding = "UTF-8")

Chron<- data.frame(matrix(ncol = 40, nrow = 245)); colnames(Chron)<- List_loc$Code; rownames(Chron)<- c(1778:2022)

for (i in c(1:40)){
  
  serie <- read.rwl(paste("C:/Users/jirka/Desktop/Data/", List_loc[i, "Code"], ".rwl", sep = ""))
  
  detrendovane.serie <- detrend(serie, method = "Spline", nyrs = 30)
  chronologie <- chron(detrendovane.serie, biweight = T)
  
  for (b in rownames(chronologie)) {
    
    Chron[rownames(Chron[b,]==rownames(chronologie[b,])), colnames(Chron) == List_loc[i, "Code"]] <- chronologie[b,1]
    
  }
}

TRI<- Chron[rownames(Chron)<2018 & rownames(Chron)>1984,]


###############################################################
#              	     3. PCA calculation      	              #
###############################################################


TRI_PISY<- TRI[,1:20]
TRI_PCAB<- TRI[,21:40]

NDVI_PISY<- NDVI[,1:20]
NDVI_PCAB<- NDVI[,21:40]

PCGA_TRI_PISY<-pcga(TRI_PISY)

PCGA_TRI_PCAB<-pcga(TRI_PCAB)

PCGA_NDVI_PISY<-pcga(NDVI_PISY)

PCGA_NDVI_PCAB<-pcga(NDVI_PCAB)

#################################
# ANOVA of PC2 and PC3 loadings
#################################

# data preparation
ROT_TRI_PISY<-as.data.frame(print(PCGA_TRI_PISY$pca$rotation)[,c(2:3)])
ROT_TRI_PCAB<-as.data.frame(print(PCGA_TRI_PCAB$pca$rotation)[,c(2:3)])
ROT_NDVI_PISY<-as.data.frame(print(PCGA_NDVI_PISY$pca$rotation)[,c(2:3)])
ROT_NDVI_PCAB<-as.data.frame(print(PCGA_NDVI_PCAB$pca$rotation)[,c(2:3)])

ROT_TRI_PISY$TOPO<- substring(rownames(ROT_TRI_PISY),3,4)
ROT_TRI_PCAB$TOPO<- substring(rownames(ROT_TRI_PCAB),3,4)
ROT_NDVI_PISY$TOPO<- substring(rownames(ROT_NDVI_PISY),3,4)
ROT_NDVI_PCAB$TOPO<- substring(rownames(ROT_NDVI_PCAB),3,4)

# ANOVA for PC2

PC2_aov<- data.frame(matrix(ncol = 6, nrow = 4)); rownames(PC2_aov)<- c("PISY_TRI", "PCAB_TRI", "PISY_NDVI", "PCAB_NDVI"); colnames(PC2_aov)<- c("PL_NS", "SS_NS", "VA_NS", "SS_PL", "VA_PL", "VA_SS")

A_TRI_PISY<- aov(ROT_TRI_PISY$PC2 ~ ROT_TRI_PISY$TOPO)
T_TRI_PISY<- TukeyHSD(A_TRI_PISY)
PC2_aov[1,]<- T_TRI_PISY[["ROT_TRI_PISY$TOPO"]]

A_TRI_PCAB<- aov(ROT_TRI_PCAB$PC2 ~ ROT_TRI_PCAB$TOPO)
T_TRI_PCAB<- TukeyHSD(A_TRI_PCAB)
PC2_aov[2,]<- T_TRI_PCAB[["ROT_TRI_PCAB$TOPO"]]

A_NDVI_PISY<- aov(ROT_NDVI_PISY$PC2 ~ ROT_NDVI_PISY$TOPO)
T_NDVI_PISY<- TukeyHSD(A_NDVI_PISY)
PC2_aov[3,]<- T_NDVI_PISY[["ROT_NDVI_PISY$TOPO"]]

A_NDVI_PCAB<- aov(ROT_NDVI_PCAB$PC2 ~ ROT_NDVI_PCAB$TOPO)
T_NDVI_PCAB<- TukeyHSD(A_NDVI_PCAB)
PC2_aov[4,]<- T_NDVI_PCAB[["ROT_NDVI_PCAB$TOPO"]]


# ANOVA for PC3

PC3_aov<- data.frame(matrix(ncol = 6, nrow = 4)); rownames(PC3_aov)<- c("PISY_TRI", "PCAB_TRI", "PISY_NDVI", "PCAB_NDVI"); colnames(PC3_aov)<- c("PL_NS", "SS_NS", "VA_NS", "SS_PL", "VA_PL", "VA_SS")

A_TRI_PISY<- aov(ROT_TRI_PISY$PC3 ~ ROT_TRI_PISY$TOPO)
T_TRI_PISY<- TukeyHSD(A_TRI_PISY)
PC3_aov[1,]<- T_TRI_PISY[["ROT_TRI_PISY$TOPO"]]

A_TRI_PCAB<- aov(ROT_TRI_PCAB$PC3 ~ ROT_TRI_PCAB$TOPO)
T_TRI_PCAB<- TukeyHSD(A_TRI_PCAB)
PC3_aov[2,]<- T_TRI_PCAB[["ROT_TRI_PCAB$TOPO"]]

A_NDVI_PISY<- aov(ROT_NDVI_PISY$PC3 ~ ROT_NDVI_PISY$TOPO)
T_NDVI_PISY<- TukeyHSD(A_NDVI_PISY)
PC3_aov[3,]<- T_NDVI_PISY[["ROT_NDVI_PISY$TOPO"]]

A_NDVI_PCAB<- aov(ROT_NDVI_PCAB$PC3 ~ ROT_NDVI_PCAB$TOPO)
T_NDVI_PCAB<- TukeyHSD(A_NDVI_PCAB)
PC3_aov[4,]<- T_NDVI_PCAB[["ROT_NDVI_PCAB$TOPO"]]


###############################################################
#              4. Calculation of climatic signal              #
###############################################################

### Loading of climadata
PISY_Temp <- read.table("PISY_temp.txt", header=F, na.strings="NA", dec=",", strip.white=TRUE)
PISY_SPEI <- read.table("PISY_SPEI.txt", header=F, na.strings="NA", dec=",", strip.white=TRUE)
PISY_clima <- list(temperature=PISY_Temp, SPEI=PISY_SPEI)

PCAB_Temp <- read.table("PCAB_temp.txt", header=F, na.strings="NA", dec=",", strip.white=TRUE)
PCAB_SPEI <- read.table("PCAB_SPEI.txt", header=F, na.strings="NA", dec=",", strip.white=TRUE)
PCAB_clima <- list(temperature=PCAB_Temp, SPEI=PCAB_SPEI)


### Climatic signal for PISY TRI

PISY_TRI_Cor<- as.data.frame(matrix(ncol = 20, nrow = 32)); colnames(PISY_TRI_Cor)<- colnames(TRI_PISY)
PISY_TRI_Sig<- as.data.frame(matrix(ncol = 20, nrow = 32)); colnames(PISY_TRI_Sig)<- colnames(TRI_PISY)

for (i in c(1:ncol(TRI_PISY))) {
  
  LOK<- as.data.frame(colnames(TRI_PISY))
  
  chron<- as.data.frame((TRI_PISY[,i])); colnames(chron)<- LOK[i,]; rownames(chron)<- rownames(TRI_PISY)
  
  korelace <- dcc(chron, PISY_clima, selection= .range(-6:9), method="correlation")
  
  PISY_TRI_Cor[,colnames(PISY_TRI_Cor) == colnames(chron)] <- korelace$coef$coef 
  PISY_TRI_Sig[,colnames(PISY_TRI_Sig) == colnames(chron)] <- korelace$coef$significant
  
}

rownames(PISY_TRI_Cor)<- paste(korelace$coef$varname ,korelace$coef$month)
rownames(PISY_TRI_Sig)<- paste(korelace$coef$varname ,korelace$coef$month)


### Climatic signal for PCAB TRI

PCAB_TRI_Cor<- as.data.frame(matrix(ncol = 20, nrow = 32)); colnames(PCAB_TRI_Cor)<- colnames(TRI_PCAB)
PCAB_TRI_Sig<- as.data.frame(matrix(ncol = 20, nrow = 32)); colnames(PCAB_TRI_Sig)<- colnames(TRI_PCAB)

for (i in c(1:ncol(TRI_PCAB))) {
  
  LOK<- as.data.frame(colnames(TRI_PCAB))
  
  chron<- as.data.frame((TRI_PCAB[,i])); colnames(chron)<- LOK[i,]; rownames(chron)<- rownames(TRI_PCAB)
  
  korelace <- dcc(chron, PCAB_clima, selection= .range(-6:9), method="correlation")
  
  PCAB_TRI_Cor[,colnames(PCAB_TRI_Cor) == colnames(chron)] <- korelace$coef$coef 
  PCAB_TRI_Sig[,colnames(PCAB_TRI_Sig) == colnames(chron)] <- korelace$coef$significant
  
}

rownames(PCAB_TRI_Cor)<- paste(korelace$coef$varname ,korelace$coef$month)
rownames(PCAB_TRI_Sig)<- paste(korelace$coef$varname ,korelace$coef$month)

### Climatic signal for PISY NDVI

PISY_NDVI_Cor<- as.data.frame(matrix(ncol = 20, nrow = 32)); colnames(PISY_NDVI_Cor)<- colnames(NDVI_PISY)
PISY_NDVI_Sig<- as.data.frame(matrix(ncol = 20, nrow = 32)); colnames(PISY_NDVI_Sig)<- colnames(NDVI_PISY)

for (i in c(1:ncol(NDVI_PISY))) {
  
  LOK<- as.data.frame(colnames(NDVI_PISY))
  
  chron<- as.data.frame((NDVI_PISY[,i])); colnames(chron)<- LOK[i,]; rownames(chron)<- rownames(NDVI_PISY)
  
  korelace <- dcc(chron, PISY_clima, selection= .range(-6:9), method="correlation")
  
  PISY_NDVI_Cor[,colnames(PISY_NDVI_Cor) == colnames(chron)] <- korelace$coef$coef 
  PISY_NDVI_Sig[,colnames(PISY_NDVI_Sig) == colnames(chron)] <- korelace$coef$significant
  
}

rownames(PISY_NDVI_Cor)<- paste(korelace$coef$varname ,korelace$coef$month)
rownames(PISY_NDVI_Sig)<- paste(korelace$coef$varname ,korelace$coef$month)


### Climatic signal for PCAB NDVI

PCAB_NDVI_Cor<- as.data.frame(matrix(ncol = 20, nrow = 32)); colnames(PCAB_NDVI_Cor)<- colnames(NDVI_PCAB)
PCAB_NDVI_Sig<- as.data.frame(matrix(ncol = 20, nrow = 32)); colnames(PCAB_NDVI_Sig)<- colnames(NDVI_PCAB)

for (i in c(1:ncol(NDVI_PCAB))) {
  
  LOK<- as.data.frame(colnames(NDVI_PCAB))
  
  chron<- as.data.frame((NDVI_PCAB[,i])); colnames(chron)<- LOK[i,]; rownames(chron)<- rownames(NDVI_PCAB)
  
  korelace <- dcc(chron, PCAB_clima, selection= .range(-6:9), method="correlation")
  
  PCAB_NDVI_Cor[,colnames(PCAB_NDVI_Cor) == colnames(chron)] <- korelace$coef$coef 
  PCAB_NDVI_Sig[,colnames(PCAB_NDVI_Sig) == colnames(chron)] <- korelace$coef$significant
  
}

rownames(PCAB_NDVI_Cor)<- paste(korelace$coef$varname ,korelace$coef$month)
rownames(PCAB_NDVI_Sig)<- paste(korelace$coef$varname ,korelace$coef$month)

###############################################################
#             5. ANOVA between site categories                #
###############################################################

### Reorganization of PISY_TRI climatic signal for ANOVA

PISY_TRI_Cor<- as.data.frame(t(PISY_TRI_Cor))
PISY_TRI_Cor$TOPO<- substring(rownames(PISY_TRI_Cor), 3)

PISY_TRI_ANOVAs<- data.frame(matrix(nrow = ncol(PISY_TRI_Cor)-1, ncol = 1)); rownames(PISY_TRI_ANOVAs)<- colnames(PISY_TRI_Cor)[1:32]; colnames(PISY_TRI_ANOVAs)<- "P val"

PISY_TRI_PHT<- data.frame(matrix(nrow = 6, ncol = nrow(PISY_TRI_ANOVAs))); colnames(PISY_TRI_PHT)<- rownames(PISY_TRI_ANOVAs)


### ANOVA calculation

for (i in c(1:ncol(PISY_TRI_Cor))) {
  
  #########  basic ANOVA   ######
  ANOVA<- aov(PISY_TRI_Cor[,i] ~ PISY_TRI_Cor$TOPO)
  
  AOV_sum<- summary(ANOVA)
  
  Pval<- AOV_sum[[1]][["Pr(>F)"]]
  
  PISY_TRI_ANOVAs[i,"P val"]<- Pval[1]
  
  ##########   post hoc tests   #########
  
  PostHoc<- TukeyHSD(ANOVA)
  PostHoc<- as.data.frame(PostHoc[["PISY_TRI_Cor$TOPO"]])
  
  PISY_TRI_PHT[,i]<- PostHoc$`p adj`
  
}

rownames(PISY_TRI_PHT)<- rownames(PostHoc)

### Reorganization of PCAB_TRI climatic signal for ANOVA

PCAB_TRI_Cor<- as.data.frame(t(PCAB_TRI_Cor))
PCAB_TRI_Cor$TOPO<- substring(rownames(PCAB_TRI_Cor), 3)

PCAB_TRI_ANOVAs<- data.frame(matrix(nrow = ncol(PCAB_TRI_Cor)-1, ncol = 1)); rownames(PCAB_TRI_ANOVAs)<- colnames(PCAB_TRI_Cor)[1:32]; colnames(PCAB_TRI_ANOVAs)<- "P val"

PCAB_TRI_PHT<- data.frame(matrix(nrow = 6, ncol = nrow(PCAB_TRI_ANOVAs))); colnames(PCAB_TRI_PHT)<- rownames(PCAB_TRI_ANOVAs)


### ANOVA calculation

for (i in c(1:ncol(PCAB_TRI_Cor))) {
  
  #########  basic ANOVA   ######
  ANOVA<- aov(PCAB_TRI_Cor[,i] ~ PCAB_TRI_Cor$TOPO)
  
  AOV_sum<- summary(ANOVA)
  
  Pval<- AOV_sum[[1]][["Pr(>F)"]]
  
  PCAB_TRI_ANOVAs[i,"P val"]<- Pval[1]
  
  ##########   post hoc tests   #########
  
  PostHoc<- TukeyHSD(ANOVA)
  PostHoc<- as.data.frame(PostHoc[["PCAB_TRI_Cor$TOPO"]])
  
  PCAB_TRI_PHT[,i]<- PostHoc$`p adj`
  
}

rownames(PCAB_TRI_PHT)<- rownames(PostHoc)

### Reorganization of PISY_NDVI climatic signal for ANOVA

PISY_NDVI_Cor<- as.data.frame(t(PISY_NDVI_Cor))
PISY_NDVI_Cor$TOPO<- substring(rownames(PISY_NDVI_Cor), 3)

PISY_NDVI_ANOVAs<- data.frame(matrix(nrow = ncol(PISY_NDVI_Cor)-1, ncol = 1)); rownames(PISY_NDVI_ANOVAs)<- colnames(PISY_NDVI_Cor)[1:32]; colnames(PISY_NDVI_ANOVAs)<- "P val"

PISY_NDVI_PHT<- data.frame(matrix(nrow = 6, ncol = nrow(PISY_NDVI_ANOVAs))); colnames(PISY_NDVI_PHT)<- rownames(PISY_NDVI_ANOVAs)


### ANOVA calculation

for (i in c(1:ncol(PISY_NDVI_Cor))) {
  
  #########  basic ANOVA   ######
  ANOVA<- aov(PISY_NDVI_Cor[,i] ~ PISY_NDVI_Cor$TOPO)
  
  AOV_sum<- summary(ANOVA)
  
  Pval<- AOV_sum[[1]][["Pr(>F)"]]
  
  PISY_NDVI_ANOVAs[i,"P val"]<- Pval[1]
  
  ##########   post hoc tests   #########
  
  PostHoc<- TukeyHSD(ANOVA)
  PostHoc<- as.data.frame(PostHoc[["PISY_NDVI_Cor$TOPO"]])
  
  PISY_NDVI_PHT[,i]<- PostHoc$`p adj`
  
}

rownames(PISY_NDVI_PHT)<- rownames(PostHoc)

### Reorganization of PCAB_NDVI climatic signal for ANOVA

PCAB_NDVI_Cor<- as.data.frame(t(PCAB_NDVI_Cor))
PCAB_NDVI_Cor$TOPO<- substring(rownames(PCAB_NDVI_Cor), 3)

PCAB_NDVI_ANOVAs<- data.frame(matrix(nrow = ncol(PCAB_NDVI_Cor)-1, ncol = 1)); rownames(PCAB_NDVI_ANOVAs)<- colnames(PCAB_NDVI_Cor)[1:32]; colnames(PCAB_NDVI_ANOVAs)<- "P val"

PCAB_NDVI_PHT<- data.frame(matrix(nrow = 6, ncol = nrow(PCAB_NDVI_ANOVAs))); colnames(PCAB_NDVI_PHT)<- rownames(PCAB_NDVI_ANOVAs)


### ANOVA calculation

for (i in c(1:ncol(PCAB_NDVI_Cor))) {
  
  #########  basic ANOVA   ######
  ANOVA<- aov(PCAB_NDVI_Cor[,i] ~ PCAB_NDVI_Cor$TOPO)
  
  AOV_sum<- summary(ANOVA)
  
  Pval<- AOV_sum[[1]][["Pr(>F)"]]
  
  PCAB_NDVI_ANOVAs[i,"P val"]<- Pval[1]
  
  ##########   post hoc tests   #########
  
  PostHoc<- TukeyHSD(ANOVA)
  PostHoc<- as.data.frame(PostHoc[["PCAB_NDVI_Cor$TOPO"]])
  
  PCAB_NDVI_PHT[,i]<- PostHoc$`p adj`
  
}

rownames(PCAB_NDVI_PHT)<- rownames(PostHoc)


###############################################################
#       6. T-test between TRI and NDVI climatic signal        #
###############################################################

### Reorganization of climatic signal for t test

# PISY TRI
PISY_TRI_Cor$LOK<- rownames(PISY_TRI_Cor)
PISY_TRI_Cor1<- gather(PISY_TRI_Cor, "Month_Clima", "COR", 1:32)
PISY_TRI_Cor1$VAR<- "TRI"
PISY_TRI_Cor1$SPECIES<- "PISY"

# PCAB TRI
PCAB_TRI_Cor$LOK<- rownames(PCAB_TRI_Cor)
PCAB_TRI_Cor1<- gather(PCAB_TRI_Cor, "Month_Clima", "COR", 1:32)
PCAB_TRI_Cor1$VAR<- "TRI"
PCAB_TRI_Cor1$SPECIES<- "PCAB"

# PISY NDVI
PISY_NDVI_Cor$LOK<- rownames(PISY_NDVI_Cor)
PISY_NDVI_Cor1<- gather(PISY_NDVI_Cor, "Month_Clima", "COR", 1:32)
PISY_NDVI_Cor1$VAR<- "NDVI"
PISY_NDVI_Cor1$SPECIES<- "PISY"

# PCAB NDVI
PCAB_NDVI_Cor$LOK<- rownames(PCAB_NDVI_Cor)
PCAB_NDVI_Cor1<- gather(PCAB_NDVI_Cor, "Month_Clima", "COR", 1:32)
PCAB_NDVI_Cor1$VAR<- "NDVI"
PCAB_NDVI_Cor1$SPECIES<- "PCAB"

# gathering everything together

T_test_data<- rbind(PISY_TRI_Cor1, PCAB_TRI_Cor1, PISY_NDVI_Cor1, PCAB_NDVI_Cor1)

T_test_data$CLIMA<- lapply(strsplit(as.character(T_test_data$Month_Clima), "\\ "), "[",1)
T_test_data$Month<- lapply(strsplit(as.character(T_test_data$Month_Clima), "\\ "), "[",2)
T_test_data<- T_test_data[,-3]
T_test_data$Code<- paste(T_test_data$CLIMA, T_test_data$Month, T_test_data$SPECIES, T_test_data$TOPO, sep = "_")

### T-test itself

Levels<- as.data.frame(unique(T_test_data$Code)); colnames(Levels)<- "Code"

T_test<- data.frame(matrix(nrow = nrow(Levels), ncol = 1)); rownames(T_test)<- Levels$Code; colnames(T_test)<- "P val"

for (i in c(1:nrow(Levels))) {
  
  Set<- subset(T_test_data, subset = T_test_data$Code == Levels[i,] )
  
  t_test<- t.test(Set$COR~Set$VAR)
  
  T_sum<- summary(t_test)
  
  Pval<- t_test[["p.value"]]
  
  T_test[i,"P val"]<- Pval[1]
  
}


###############################################################
#               6. Correlation of TRI and NDVI                #
###############################################################

COR_TRI_NDVI<- data.frame(matrix(ncol = 2, nrow = ncol(TRI))); rownames(COR_TRI_NDVI)<- colnames(TRI); colnames(COR_TRI_NDVI)<- c("COR", "SIG")

for (i in c(1:40)) {
  
  COR_t<- cor.test(TRI[,i], NDVI[,i])
  
  COR_TRI_NDVI[i,"COR"]<- COR_t$estimate
  COR_TRI_NDVI[i,"SIG"]<- COR_t$p.value
  
}


