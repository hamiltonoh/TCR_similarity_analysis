# This script contains demonstrates T cell receptor (TCR) similarity analysis using Levenshtein Distance.
# Of all TCRs in the dataset, we retain clonal TCRs with unambiguous TCR alpha (TCRa) and beta (TCRb) sequences.
# Downstream analysis of TCR similarity is performed on the full length TCRab sequence.
# Analysis can also be done with only TCRb, which contains the antigen binding domain. 
# For TCRb analysis, retain clonal TCRs with unambiguous beta sequences.

setwd("/home/hoh3/Hamilton/CSF_scRNA_analysis/proj_CSF_methods")

require(RecordLinkage)
require(gplots)
require("qgraph")  


# Load homemade function
mat2edge<-function(x,cutoff){
  if(missing(cutoff)==T){cutoff<-0}
  res=NULL
  for(i in 1:nrow(x)){
    for(k in 1:ncol(x)){
      if(x[i,k]>cutoff){
        res=rbind(res,data.frame(X=rownames(x)[i],Y=colnames(x)[k],weight=x[i,k],stringsAsFactors = F))
      }
    }
    print(i)
  }
  res
}


# Load data containing clonotype info for all samples
df=read.csv("data/csf1-24_clonotypes.csv",stringsAsFactors = F)
head(df)
unique(df$ID)

# Save ID and Diagnosis information for network visualization
ForNetwork1=data.frame(Diagnosis=df$Diagnosis,ID=df$ID,stringsAsFactors = F)
ForNetwork1=ForNetwork1[which(duplicated(paste(ForNetwork1[,1],ForNetwork1[,2],sep="_"))==F),]
colnames(ForNetwork1)=c("x","y")
head(ForNetwork1)


# Filter unwanted TCRs
# keep TCRs with both tcrA and b - i.e remove those with empty values
dim(df)
df=df[which(df$CDR3b!="" & df$CDR3a!=""),]  #for tcrab analysis
#df=df[which(df$CDR3b!=""),]  #for tcrb analysis

# keep clones only (frequency > 1)
df=df[which(df$Frequency>1),]

# remove TCRs with ambiguous sequences (both for tcrab analysis, only cdr3b for tcrb analysis)
df=df[-grep(";",df$CDR3b),]
df=df[-grep(";",df$CDR3a),]
dim(df)

# create clonotype_id.long as unique ID for each clonotype
df$clonotype_id.long=paste(df$CDR3a,df$CDR3b,sep="_")
head(df)

# List of patients sharing identical TCRs
df[which(df$clonotype_id.long %in% df$clonotype_id.long[which(duplicated(df$clonotype_id)==T)]),]  #identical TCRab
df[which(df$CDR3b %in% df$CDR3b[which(duplicated(df$CDR3b)==T)]),]  #identical TCRb


# Calculate levenshtein similarity between all possible pairs of TCRab
res.clonotype_id.long=NULL
for(i in 1:length(df$clonotype_id)){
  res.clonotype_id.long=rbind(res.clonotype_id.long,levenshteinSim(df$clonotype_id.long[i],df$clonotype_id.long))
}

# rename rows/cols
colnames(res.clonotype_id.long)<-df$clonotype_id.long
rownames(res.clonotype_id.long)<-df$clonotype_id.long
head(res.clonotype_id.long[0:5, 0:5])

# plot histogram of levenshtein similarity score distribution
toHist=res.clonotype_id.long
toHist[upper.tri(toHist)]<-NA
diag(toHist)<-NA
hist(toHist,xlab="Similarity between TCRs",main="")
hist(toHist,xlab="Similarity between TCRs",main="",xlim=c(.8,1),ylim = c(0,30),breaks = 100)


# Make heatmap of TCR similarity
res.clonotype_id.long.2=res.clonotype_id.long
res.clonotype_id.long.2[res.clonotype_id.long.2<0.8]<-0  #only visualize pairs with levenshteinSim>0.8
diag(res.clonotype_id.long.2)<-0
res.clonotype_id.long.2=res.clonotype_id.long.2[which(apply(res.clonotype_id.long.2,1,function(x) sum(x>0,na.rm=T))>0),
                                                which(apply(res.clonotype_id.long.2,1,function(x) sum(x>0,na.rm=T))>0)]

toHeatmap=res.clonotype_id.long.2
pairs.breaks <- seq(0, 1, by=0.01)
mycol <- colorpanel(n=length(pairs.breaks)-1,low="darkslateblue",mid="red",high="yellow")
tmp=(heatmap.2(toHeatmap,
               cexRow=.2,cexCol=.2,
               trace="none",
               dendrogram="both",
               breaks=pairs.breaks,
               col=mycol,
               Rowv=T,key=F,
               Colv=T,
               lhei=c(0.2,4),
               lwid=c(.2,3)
))


# TCR node network cutoff 0.8 using qraph
{
  # extract link between TCRs and patients ID for network
  ForNetwork2=data.frame(df$ID,df$clonotype_id.long,stringsAsFactors = F)
  ForNetwork2=ForNetwork2[which(duplicated(paste(ForNetwork2[,1],ForNetwork2[,2],sep="_"))==F),]
  colnames(ForNetwork2)=c("x","y")
  head(ForNetwork2)

  # add similar TCRs to the network
  ForNetwork3=res.clonotype_id.long
  ForNetwork3[upper.tri(ForNetwork3)]<-NA
  diag(ForNetwork3)<-NA
  head(ForNetwork3[,1:5])
  ForNetwork3[is.na(ForNetwork3)]<-0
  
  ForNetwork3=mat2edge(x=ForNetwork3, cutoff=0.8)  #Set levenshteinSim cutoff for visualization
  colnames(ForNetwork3)=c("x","y","weight")
  head(ForNetwork3)
  dim(ForNetwork3)
  
  
  # network
  {
    toQgraph=rbind(ForNetwork1,ForNetwork2)
    toQgraph$weight=5
    toQgraph=rbind(toQgraph,ForNetwork3)

    # customize qgraph 
    {
      # color of disease groups
      l=unique(c(toQgraph[,1],toQgraph[,2]))
      col.tmp=rep(adjustcolor("grey85",.3),length(l))
      col.tmp[which(l %in% c("AD"))]<-adjustcolor("red",.8)
      col.tmp[which(l %in% c("HC"))]<-adjustcolor("blue",.8)
      col.tmp[which(l %in% c("MCI"))]<-adjustcolor("orange",.8)
      col.tmp[which(l %in% c("PD"))]<-adjustcolor("green4",.8)

      # color of subjects nodes
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="AD")]))]<-adjustcolor("red",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="HC")]))]<-adjustcolor("blue",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="MCI")]))]<-adjustcolor("orange",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="PD")]))]<-adjustcolor("green4",.3)

      # size of nodes
      size.tmp=rep(1,length(l))
      size.tmp[which(l %in% c("AD","HC","MCI","PD"))]<-15
      size.tmp[grep("CSF",l)]<-10
      table(size.tmp)
      
      # color of edges
      line.tmp=rep("black",nrow(toQgraph))
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("HC"))]<-"blue"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("MCI"))]<-"orange"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("PD"))]<-"green4"
      
      labels.cex.tmp=l
      labels.cex.tmp=rep(0.00000000001,length(l))
      labels.cex.tmp[which(l %in% c("AD","HC","MCI","PD"))]<-1.5
      labels.cex.tmp[grep("CSF",l)]<-.8
    }
  }
  
  # plot
  toQgraph$weight<-1
  qgraph(toQgraph, color=col.tmp, vsize=size.tmp/2, edge.color=line.tmp, labels=TRUE, label.cex=labels.cex.tmp, directed=F, layout="spring")
}





# done
