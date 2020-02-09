########################################################################################################
## title: "Construction of Gene Regulatory Network for tissue-specific genes"
## author: "Madara Hetti-Arachchilage; mhettiar 'at' illinois.edu"
## date: "February 3, 2020"
#############################################################################################

get_reg_Assoc <- function(Enriched_CRE, GeneList){
  
  for (i in 1:length(uniq)){
    Enriched_CRE_sub <- subset(Enriched_CRE, Enriched_CRE$Gene_ID == uniq[i])
    Motif_sub <- subset(Motif, Motif$Sequence.ID == uniq[i] & Motif$Family %in% Enriched_CRE_sub$Family)
    TF <- unique(Motif_sub[Motif_sub$Sequence.ID == uniq[i],]$Family)
    tt <- TF_Family[TF_Family$Family %in% TF,]
    df <- cbind(rep(uniq[i], times= nrow(tt)), tt)
    network <- rbind(network,df)
    print(i)
  }
  
  final_reg <- as.data.frame(cbind(as.character(network[,2]), 
                                   rep("interacts", times= nrow(network)), 
                                   as.character(network[,1])))
  return(final_reg)
}
#############################################################################################

get_Gene_Reg_Network <- function(reg_Assoc, Corr_Net, Output_File){

  temp1 <- dplyr::inner_join(reg_Assoc, Corr_Flt, by=c("Gene1" = "Gene1", "Gene2" = "Gene2"))
  temp2 <- dplyr::inner_join(reg_Assoc, Corr_Flt, by=c("Gene2" = "Gene1", "Gene1" = "Gene2"))
  GRN <- rbind(temp1,temp2)
  GRN <- GRN[!(duplicated(GRN)| duplicated(GRN, fromLast=TRUE)),]
  GRN <- GRN[which(GRN$Gene1 != GRN$Gene2),]
  GRN <- GRN[complete.cases(GRN), ]
  GRN$rg <- "rg"
  GRN$rg <- paste(GRN$rg, round(as.numeric(GRN$Corr),digits = 1), sep = "")
  GRN <- GRN[,c(1,2,3,5,6,9)]
  return(GRN)
}

#############################################################################################
