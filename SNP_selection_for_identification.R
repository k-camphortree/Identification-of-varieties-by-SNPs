# argument 1: Input SNP information file
# argument 2: MAF (Minor Allele frequency; 0-50 [%]) 
# argument 3: Minimum distance between markers [bp]
# argument 4: The number of calculation for barcodeing of all samples using random marker sets [times]
# argument 5: Output folder
# argument 6: Minimum number of markers used for identification
# argument 7: Maximum number of markers used for identification
# 
# Input file: Input detected SNP information from Ion Proton and CLC Genomics Workbench
# 1st column: Chromosome
# 2nd column: Region (Physical position)
# 3rd column: Type (SNV, Insertion, Deletion)
# 4th column: Reference (Reference allele)
# 5th columns: Allele (Alternative allele)
# 6th columns: Reference.allele
# 7th columns: Length
# 8th columns: Linkage
# 9th columns: Zygosity
# 10th column: Count
# 11th columns: Coverage
# 12th columns: Frequency
# 13th columns: Probability
# 14th columns: Forward.read.count
# 15th columns: Reverse.read.count
# 16th columns: Forward.reverse.balance
# 17th columns: Average.quality
# 18th columns: Origin.tracks
# 19th columns: Sample.count
# 20th columns: Total.number.of.samples
# 21th columns: Sample.frequency
# 
# Output file:
# (1) Output_marker_sets.csv: Marker information used for random marker selection
# (2) Markertable_with_X_markers.csv: Marker sets which successfully identified all samples
# (3) Markertable_usage_count.csv: Number of marker sets
# (4) Markertable_usage_markers.csv: Usage count of each marker for the marker sets
# 
# Example: Rscript --vanilla --slave Identification-of-varieties-by-SNPs.R $HOME/Input/Input.csv 40 500000 500000 $HOME/Output/ 10 25
#
# 
# 

#Argument
in_file = commandArgs(trailingOnly=TRUE)[1]
MAF = commandArgs(trailingOnly=TRUE)[2]
Distance = commandArgs(trailingOnly=TRUE)[3]
Caltimes = commandArgs(trailingOnly=TRUE)[4]
out_folder = commandArgs(trailingOnly=TRUE)[5]
M_min = commandArgs(trailingOnly=TRUE)[6]
M_max = commandArgs(trailingOnly=TRUE)[7]


input <- read.csv(in_file)
SNP_list <- subset(input, Type=="SNV"&Reference.allele=="No"&Zygosity=="Homozygous"&Sample.frequency>=MAF&Sample.frequency<=(100-MAF))
SNP_list[,2]<- as.numeric(as.character(SNP_list[,2]))
rem <- c()
for (i in 1:nrow(SNP_list)){
  near <- subset(SNP_list, Chromosome==(as.numeric(as.character(SNP_list[i, 1])))&Region>=(as.numeric(as.character(SNP_list[i, 2])))&Region<=(as.numeric(as.character(SNP_list[i, 2]))+Distance))
  rem <- c(rem,rownames(near[-1,]))
}
rem2 <- unique(rem)
SNP_list_without_neighbor <- SNP_list[!(rownames(SNP_list) %in% rem2), ]

sampleno <- strsplit(as.character(SNP_list_without_neighbor$Origin.tracks), ", ")
mat <- matrix(rep(0,96*nrow(SNP_list_without_neighbor)),ncol=96,nrow=nrow(SNP_list_without_neighbor))
for (i in 1:nrow(SNP_list_without_neighbor)){
  mat[i,as.numeric(sampleno[[i]])] <- 1
}

SNP_list_without_neighbor_allele <- cbind(SNP_list_without_neighbor,mat[,1:n_sample])
colnames(SNP_list_without_neighbor_allele)[22:103]<-paste("Sample",colnames(SNP_list_without_neighbor_allele)[22:103],sep="")
SNP_list_without_neighbor_allele_marker <- cbind(paste("M",formatC(1:nrow(SNP_list_without_neighbor),width=2,flag="0"),sep=""),SNP_list_without_neighbor_allele)
colnames(SNP_list_without_neighbor_allele_marker)[1]<- c("Marker")
write.csv(SNP_list_without_neighbor_allele_marker, paste(out_folder,"Output_marker_sets.csv",sep=""), row.names = F)
SNPmatrix <- SNP_list_without_neighbor_allele_marker[,23:ncol(SNP_list_without_neighbor_allele_marker)]
rownames(SNPmatrix) <- SNP_list_without_neighbor_allele_marker[,1]

SNPmatrix_t<-data.frame(t(SNPmatrix))
SNPmatrix_s<-cbind(SNPmatrix_t,colnames(SNPmatrix))
colnames(SNPmatrix_s)<-c(as.character(SNP_list_without_neighbor_allele_marker$Marker),"Line")


library(tidyverse)
usage_count <- data.frame()
for (n_marker in M_min:M_max){
  markertable<-matrix(nrow=1,ncol=n_marker)
  markertable<-markertable[-1,]
  classification_no <- c()
  for (i in 1:Caltimes){
    Selected_matrix <- SNPmatrix_s[,sample(1:nrow(SNPmatrix_s), n_marker, replace = F)]
    Selected_matrix[,n_marker+1] <- tidyr::unite_(Selected_matrix, paste(colnames(Selected_matrix), collapse="_"), colnames(Selected_matrix))
    classification_no <- c(classification_no,length(unique(Selected_matrix[,n_marker+1])))
    if(length(unique(Selected_matrix[,n_marker+1]))==n_sample){
      markertable<-rbind(markertable, matrix(sort(colnames(Selected_matrix)[1:n_marker]),nrow=1))
    }
  }
  markertable<-data.frame(markertable)
  markertable[,n_marker+1] <- tidyr::unite_(markertable, paste(colnames(markertable), collapse="_"), colnames(markertable))
  colnames(markertable)[n_marker+1] <- "Barcode"
  write.csv(markertable, paste(out_folder, "Markertable_with_",m,"_markers.csv",sep=""), row.names = F)
  unique(markertable$Barcode)
  print(max(classification_no))
  if(nrow(markertable)>0){
    marker <- c(as.matrix(markertable[,1:n_marker]))
    markertable2 <- data.frame(cbind(n_marker,names(table(marker)),as.matrix(table(marker))))
    usage_count <- rbind(usage_count,markertable2)
  }
}
write.csv(classification_no, paste(out_folder, "Markertable_usage_count.csv", sep = ""), row.names = F)
write.csv(usage_count, paste(out_folder, "Markertable_usage_markers.csv", sep = ""), row.names = F)
