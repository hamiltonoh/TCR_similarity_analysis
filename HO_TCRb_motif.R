# This script identifies shared motifs in TCRb sequences among all patient samples, then identifies information
# (i.e. antigen specificity, T cell type) about TCRs containing specific motifs from the McPAS-TCR database.
# Of all TCRs in our dataset, we retain clonal TCRs with unambiguous TCR beta (TCRb) sequences.
# This script was developed by Benoit Lehallier.


setwd("/home/hoh3/Hamilton/CSF_scRNA_analysis/proj_CSF_methods")

library(tidyverse)
library(ggplot2)
library(cowplot)


# open data
df=read.csv("data/csf1-24_clonotypes.csv",stringsAsFactors = F)

# remove cells without tcrb
dim(df)
df=df[which(df$CDR3b!=""),]

# keep clones only (frequency > 1)
df=df[which(df$Frequency>1),]

# remove ambiguous TCRs - with ;
df=df[-grep(";",df$CDR3b),]
dim(df)


{
  # Identify all motifs using sliding k-mer approach
  nrow(df)
  for(jjj in 3:8){
    nletters=jjj
    l.list=list()
    # loop over each tcrb
    for(k in 1:nrow(df)){
      tcr_tmp=df$CDR3b[k]
      tcr_tmp.save=tcr_tmp
      # Exclude the 1st 2 and last 2 characters from the analysis
      tcr_tmp=substr(tcr_tmp,start = 3,nchar(tcr_tmp)-2)
      tcr_split_tmp=NULL
      if(nchar(tcr_tmp)==nletters){  #if length of tcr == length motif, then tcr_split_tmp = tcr sequence
        tcr_split_tmp=tcr_tmp
      }
      if(nchar(tcr_tmp)>nletters){   #if length of tcr > length motif, then loop through remaining kmers
        for(i in 1:(nchar(tcr_tmp)-(nletters)+1)){
          tcr_split_tmp=c(tcr_split_tmp,substr(x = tcr_tmp,i,i+nletters-1))  #vector of kmer sequences
        }
      }
      l.list[[k]]<-tcr_split_tmp  #compile kmer sequences in list, organized by tcr. if length of tcr < length motif, will add nothing
    }
    names(l.list)<-df$CDR3b
    assign(paste("l.list.",nletters,sep=""),value = l.list)
  }
  
  head(l.list.7)
  table(unlist(l.list.7))
  
  # Create a frequency table
  Freq.Table.kmers=data.frame(l=3,table(unlist(l.list.3)),stringsAsFactors = F)
  Freq.Table.kmers=rbind(Freq.Table.kmers,data.frame(l=4,table(unlist(l.list.4)),stringsAsFactors = F))
  Freq.Table.kmers=rbind(Freq.Table.kmers,data.frame(l=5,table(unlist(l.list.5)),stringsAsFactors = F))
  Freq.Table.kmers=rbind(Freq.Table.kmers,data.frame(l=6,table(unlist(l.list.6)),stringsAsFactors = F))
  Freq.Table.kmers=rbind(Freq.Table.kmers,data.frame(l=7,table(unlist(l.list.7)),stringsAsFactors = F))
  Freq.Table.kmers=rbind(Freq.Table.kmers,data.frame(l=8,table(unlist(l.list.8)),stringsAsFactors = F))
  Freq.Table.kmers=data.frame(Freq.Table.kmers,stringsAsFactors = F)
  Freq.Table.kmers=Freq.Table.kmers[which(Freq.Table.kmers$Freq != 1),]
  
  colnames(Freq.Table.kmers) <- c("Length", "Motif", "Frequency")
  head(Freq.Table.kmers)
}


# Plot - frequency vs motif length
ggplot(data = Freq.Table.kmers, mapping = aes(x=Length, y=Frequency)) +
  geom_point(size=0.5) +
  scale_y_continuous(trans='log10', breaks=c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256)) +
  geom_jitter(width = 0.05, height = 0.01, size = 0.5) +
  ylab("Frequency") + 
  xlab("Length") +
  ggtitle("CDR3b Motifs") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white", color="white"),
        panel.background = element_rect(fill="white", color="black", size = 1))


# Identify specific motifs of size 7,8
Freq.Table.kmers <- Freq.Table.kmers[order(Freq.Table.kmers$Length, Freq.Table.kmers$Frequency, decreasing=TRUE),]
dim(Freq.Table.kmers)  #4195
head(Freq.Table.kmers)

Freq.Table.kmers_7 <- Freq.Table.kmers[which(Freq.Table.kmers$Length == 7),]
dim(Freq.Table.kmers_7)  #127 shared motifs of length 7
head(Freq.Table.kmers_7)

Freq.Table.kmers_8 <- Freq.Table.kmers[which(Freq.Table.kmers$Length == 8),]
dim(Freq.Table.kmers_8)  #29 shared motifs of length 8
head(Freq.Table.kmers_8)
  
{ 
  # Search TCR motif in McPAS-TCR database to identify candidate antigen specificity of motif containing TCRs
  # First, identify TCRs in database that contain motifs identified above.
  # http://friedmanlab.weizmann.ac.il.stanford.idm.oclc.org/McPAS-TCR/
  McPAS.TCR=read.csv("data/McPAS-TCR.csv",stringsAsFactors = F)
  head(McPAS.TCR)
  dim(McPAS.TCR)  ##21885 x 30
  colnames(McPAS.TCR)
  
  # focus on human TCRs
  McPAS.TCR=McPAS.TCR[which(McPAS.TCR$Species=="Human"),] 
  head(McPAS.TCR)
  dim(McPAS.TCR)  ##18179 x 30
  
  # Find all relevant data on TCRs containing any shared motif (frequency > 1) in our dataset
  res.TCR=Freq.Table.kmers
  res.TCR$Pathology=rep(NA,nrow(res.TCR))
  res.TCR$Category=rep(NA,nrow(res.TCR))
  res.TCR$T.Cell.Type=rep(NA,nrow(res.TCR))
  res.TCR$T.cell.characteristics=rep(NA,nrow(res.TCR))
  for(i in 1:nrow(res.TCR)){
    res.TCR$Category[i]<-paste(McPAS.TCR$Category[grep(res.TCR$Motif[i], McPAS.TCR$CDR3.beta.aa)],collapse = " | ")
    res.TCR$Pathology[i]<-paste(McPAS.TCR$Pathology[grep(res.TCR$Motif[i],McPAS.TCR$CDR3.beta.aa)],collapse = " | ")
    res.TCR$T.Cell.Type[i]<-paste(McPAS.TCR$T.Cell.Type[grep(res.TCR$Motif[i],McPAS.TCR$CDR3.beta.aa)],collapse = " | ")
    res.TCR$T.cell.characteristics[i]<-paste(McPAS.TCR$T.cell.characteristics[grep(res.TCR$Motif[i],McPAS.TCR$CDR3.beta.aa)],collapse = " | ")
  }
  
  #write.csv(res.TCR,"data/csf1-24_tcrb_motifs-in-McPAS.csv",row.names=F)
}  

# Search information about specific motifs
#df.motif <- read.csv("data/csf1-24_tcrb_motifs-in-McPAS.csv", stringsAsFactors = F)
df.motif <- res.TCR
head(df.motif)

#GATNEKL motif
motif.search <- as.character(Freq.Table.kmers_7[2,]$Motif)
motif.search 

# patient data of TCRs with GATNEKL motif
patientdf_GATNEKL <- df[which(str_detect(df$CDR3b, motif.search)),]
patientdf_GATNEKL

# info of GATNEKL motif from McPAS-TCR database
df.motif.search <- df.motif[which(df.motif$Motif == motif.search),]
df.motif.search$Motif
df.motif.search$Pathology
pathogens.tmp <- strsplit(df.motif.search$Pathology, " | ", fixed=T)
table(pathogens.tmp)  #pathogens targeted by T cells with GATNEKL motif

df.motif.search$T.Cell.Type
TcellType.tmp <- strsplit(df.motif.search$T.Cell.Type, " | ", fixed=T)
table(TcellType.tmp)  #T cell types that target above pathogens

  
    
  