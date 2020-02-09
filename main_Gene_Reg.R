########################################################################################################
## title: "Construction of Gene Regulatory Network for tissue-specific genes"
## author: "Madara Hetti-Arachchilage; mhettiar 'at' illinois.edu"
## date: "February 3, 2020"
########################################################################################################

##### Importing libraries and custom scripts 
# Loading R source code which has defined custom functions needed for later steps
# Loading R/Biconductor packages (checks if the requested package is already installed in 
# the environment. If not, the package is installed from the appropriate source and then 
# loaded into the workspace

source("./get_CRE_Enrichment.R")
source("./get_Coexpression_Network.R")
source("./get_Gene_Reg_Network.R")
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(qvalue)
library(psych)
########################################################################################################

##### Step 1.1: Calculate Average Frequency of TFBS motifs in the genome ##### 
# The input files required for this step is a text file containing positions of cis-regulatory motifs present in 
# promoters (e.g. 1KB upstream of TSS) of all genes in the interested genome. This file should contain
# these columns; "Sequence.ID", "TFBS.ID", "Posistion", "Strand",	"Binding.sequence", "Family".
# This step save a output file containing Average Frequency of TFBS motifs present in the genome
 

TFBS_Genome_File <- read.table("./Data/PlantPAN1KBGenesTFIDs.txt", 
                         header=TRUE, sep="\t", fill = TRUE)
AveragesFile <- AvgFrequency_TFBS_Occurence(TFBS_Genome_File)
write.table(AveragesFile,"./Data/AvgFrequencyTFIDOccurencePlantPan1KB.txt",
            quote=F,row.names = F, sep="\t")
########################################################################################################

##### Step 1.2: Calculate Fold change for occurence of TFBS motifs in promoters of given list of genes  
# Input 1: a text file containing positions of cis-regulatory motifs present in  
# promoters (e.g. 1KB upstream of TSS) of all genes in the interested genome. This file should contain
# these columns; "Sequence.ID", "TFBS.ID", "Posistion", "Strand",	"Binding.sequence", "Family". 
# Input 2: contains Average Frequency of TFBS motifs present in the genome
# Input 3: contains a list of interested genes 

AveragesFile <- read.table("Data/AvgFrequencyTFIDOccurencePlantPan1KB.txt", 
                               header=TRUE, sep="\t", fill = TRUE)
GeneList <- read.table("Data/GeneList.txt", header= TRUE, sep="\t", as.is = T)
TFBS_Family <- read.csv("Data/TFBS_Family_TFID.csv", header = T, as.is = T)
TFBS_Family <- TFBS_Family[!duplicated(TFBS_Family$Matrix.ID),]

# Enriched CREs for given gene list based on FC > 1.1
enrich_CRE <- data.frame(matrix(data = NA, ncol = 3))
colnames(enrich_CRE) <- c("Gene_ID", "TFBS_ID", "FC")
for (i in 1:length(GeneList$Gene_ID)){
print(i)
t <- enrich_TFBS_Family(TFBS_Genome_File, GeneList$Gene_ID[i], AveragesFile, 1.1)
enrich_CRE <- rbind(enrich_CRE, t)
}

enrich_CRE_Family <- merge(enrich_CRE, TFBS_Family, by.x = "TFBS_ID", by.y = "Matrix.ID")
write.csv(enrich_CRE_Family, "Output/Enrichment/CRE_Enrichment_FC_1.1_GeneList.csv", row.names = F)
################################################################################################

##### Step 1.3 Create reg_Assoc file containing all TF-target gene relationship
input1 <- "Data/PlantPAN1KBGenesTFIDs_bothStrand_PlantTFDBFamilyUpdated.txt"
input2 <- "Data/CRE_Enrichment_FC_1.1_GeneList.csv"
input3 <- "Data/GeneList.txt"
output <- "Output/Enrichment/GeneList_reg_intr_CRE_Enrichment_FC1.1.txt"

Motif <- read.table(input1, header=T, sep="\t") 
Enriched_CRE <- read.csv(input2, header = T, as.is = T)
TF_Family <- read.table(input3, header = T, sep= "\t", as.is = T)
SP <- read.csv(input4, header = T, as.is = T)

uniq <- unique(c(TF_Family$Gene_ID, SP$Gene_ID))
network <- data.frame()
reg_Assoc <- get_reg_Assoc(Enriched_CRE, uniq)
write.table(reg_Assoc, output, sep="\t", quote=F, row.names = F)

########################################################################################################

##### Step 2.1 create a count file for the given geneList with individual samples including replicates
### Read input files
count_file <- "Data/TPM.csv"
input1 <- "Data/GeneList.txt"
count <- read.csv(count_file, header=T, sep=",", row.names = 1) 
GeneList <- read.table(input1, header = T, sep= "\t", as.is = T)

df <- subset(count, row.names(count) %in% uniq)
# Define Stages
stage_order <- c("S1", "S2", "S3", "S4")
# Define Tissues
tissue <- c("tis1", "tis2", "tis3", "tis4", "tis5")

# Generate stem expression profiles by each stage for expressed genes (5 TPM cutoff)
for (i in 1:length(stage_order)){
  
  edata <- df[,grepl(stage_order[i], colnames(df))]
  edata <- edata[,grepl(tissue[1], colnames(edata))]
  cc <- edata[rowSums(edata >= 5) > 1,]
  edata <- data.frame(as.character(rownames(cc)), cc)
  colnames(edata)[1] <- "GeneID"
  
  output1 <- paste("Output/Coexpression/GeneList_TPM_tis1_", stage_order[i], ".txt", sep= "")
  write.table(edata, output1, sep="\t", quote=F, row.names = F)
  rm(edata, cc, output1)
}
########################################################################################################

#### Step 2.2 create a peasrson correlation coexpression network
stage <- c("S1", "S2", "S3", "S4")

for (i in 1:5){
  print(i)
  output <- paste("Output/Coexpression/Cor_TPM", stage[i], ".cor", sep="")
  input <- paste("Output/Coexpression/GeneList_TPM_tis1_", stage[i], ".txt", sep="")
  # cpm values of all genes for allcondirions
  Data <- read.table(input, header = T, sep = "\t")
  row.names(Data) <- Data[, 1]
  Cor_Network <- cor.test.mk.R(Data)
  write.table(
    Cor_Network,
    file = output,
    quote = F,
    sep = ",",
    row.names = F
  )
}
########################################################################################################

#### Step 2.3 Calculate Mutual Rank based on pearson correlation coexpression network
require(data.table)
stage_order <- c("S1", "S2", "S3", "S4")

for (x in 1:5){
  output <- paste("Output/Coexpression/Corr_MRRank_", stage[x], ".cor", sep="")
  input <-paste("Output/Coexpression/Cor_TPM", stage[x], ".cor", sep="")
  MR_Network<- mutual_Rank(input)
  write.table(MR_Network,file=output,quote=F,sep=",", row.names = FALSE)
}
########################################################################################################

#### Step 3.1 Merge Coexpression network with gene regulatory network
# Correlation files
input1 <- list.files(path = paste(getwd(), "Output/Coexpression/", sep=""), 
                     pattern = "Corr_MRRank_*", full.names = T)

# regulatory interaction files
input2 <- "Enrichment/GeneList_reg_intr_CRE_Enrichment_FC1.1.txt"
final_reg <- read.csv(file = input2, header = T, sep= "\t")

# Output files
output1 <- data.frame()
output1 <- gsub("Output/Coexpression/", "Output/GRN/", input1, fixed=TRUE)
output1 <- gsub(".cor", ".txt", output1, fixed=TRUE)
output1 <- gsub("Corr_MRRank", "GRN_0.7_0.05", output1, fixed=TRUE)

for (i in 1:length(input1)){
  
  Cor1 <- read.table(input1[i], header = T, sep= ",", as.is = T)
  colnames(Cor1) <- c("Gene1", "Gene2", "MR", "Corr", "P-Value", "PCC1", "PCC2", "MR_Scale")
  Cor_TS_All <- subset(Cor1, Cor1$Corr >= 0.7 & Cor1$`P-Value` <= 0.05 & Cor1$MR_Scale >= 0.7)
  colnames(final_reg) <- c("Gene1", "rg", "Gene2")
  
  temp1 <- dplyr::inner_join(final_reg, Cor_TS_All, by=c("Gene1" = "Gene1", "Gene2" = "Gene2"))
  temp2 <- dplyr::inner_join(final_reg, Cor_TS_All, by=c("Gene2" = "Gene1", "Gene1" = "Gene2"))
  Cor_TS_Final <- rbind(temp1,temp2)
  Cor_TS_Final <- Cor_TS_Final[!(duplicated(Cor_TS_Final)| duplicated(Cor_TS_Final, fromLast=TRUE)),]
  Cor_TS_Final <- Cor_TS_Final[which(Cor_TS_Final$Gene1 != Cor_TS_Final$Gene2),]
  Cor_TS_Final <- Cor_TS_Final[complete.cases(Cor_TS_Final), ]
  Cor_TS_Final$rg <- "rg"
  Cor_TS_Final$rg <- paste(Cor_TS_Final$rg, round(as.numeric(Cor_TS_Final$Corr),digits = 1), sep = "")
  Cor_TS_Final <- Cor_TS_Final[,c(1,2,3,5,6,9)]
  write.table(Cor_TS_Final, output1[i], sep="\t", quote=F, row.names = F)
}

########################################################################################################
