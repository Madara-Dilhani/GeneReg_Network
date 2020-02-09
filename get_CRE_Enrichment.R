########################################################################################################
## title: "Construction of Gene Regulatory Network for tissue-specific genes"
## author: "Madara Hetti-Arachchilage; mhettiar 'at' illinois.edu", Kavya Kannan
## date: "February 3, 2020"
########################################################################################################
AvgFrequency_TFBS_Occurence <- function(TFBS_Pos_Genome){

TempIDs <- TFBS_Pos_Genome$TFBS.ID
TempIDs <- unique(TempIDs)
AveragesFile <- data.frame()

for(i in 1:length(TempIDs)){
  TempFile <- TFBS_Pos_Genome[TFBS_Pos_Genome$TFBS.ID == TempIDs[i],]
  TF <- count(TempFile, TempFile$Sequence.ID)
  TempAvg <- mean(TF$n)
  Temp_n <- cbind(as.character(TempIDs[i]),TempAvg, TempFile$Family)
  AveragesFile <- rbind(AveragesFile, Temp_n)
  print(i)
}

colnames(AveragesFile)[1] <- "TFBS_ID"
colnames(AveragesFile)[2] <- "AvgNumBindingSites"
colnames(AveragesFile)[3] <- "TFFamily"
AveragesFile <- AveragesFile[order(AveragesFile$AvgNumBindingSites),]
return(AveragesFile)
}

########################################################################################################
enrich_TFBS_Family <- function(TFBS_Pos_Genome, Gene, AveragesFile, FoldChange){

Enrich <- data.frame(matrix(ncol = 3))
colnames(Enrich) <- c("Gene_ID", "TFBS_ID", "FC")
TFBS_Pos_Gene <- subset(TFBS_Pos_Genome, TFBS_Pos_Genome$Sequence.ID == Gene)
TF <- count(TFBS_Pos_Gene, TFBS_Pos_Gene$TFBS.ID)
colnames(TF) <- c("TFBS_ID", "Count")

for (y in 1:nrow(TF)){
  FC <- (TF[y,2])/(AveragesFile[AveragesFile$TFBS_ID == TF[y,1]$TFBS_ID,][,2])
  Temp_n <- cbind(as.character(Gene),as.character(TF[y,1]$TFBS_ID),as.numeric(FC))
  colnames(Temp_n) <- c("Gene_ID", "TFBS_ID", "FC")
  Enrich <- rbind(Enrich, Temp_n)
}

t <-subset(Enrich, Enrich$FC > FoldChange)
return(t)
}
########################################################################################################