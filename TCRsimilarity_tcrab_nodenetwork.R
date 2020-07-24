# Tcell receptor analysis

setwd("/home/hoh3/Hamilton/CSF_scRNA_analysis/proj_CSF_methods")


# open data
df=read.csv("../metadata/all_CSF_TCRs_for_Benoit_2-21-20.csv",stringsAsFactors = F)
head(df)
unique(df$ID)

# subset CSF1-24 (except CSF22)
df <- df[which(df$ID %in% c("CSF1", "CSF2", "CSF3", "CSF4", "CSF5", "CSF6", "CSF7", "CSF8", "CSF9", "CSF10",
                            "CSF11", "CSF12", "CSF13", "CSF14", "CSF15", "CSF16", "CSF17", "CSF18", "CSF19", "CSF20",
                            "CSF21", "CSF22", "CSF23", "CSF24")),]

df <- df[-which(df$ID %in% c("CSF22", "CSF25", "CSF26", "CSF27", "CSF28", "CSF29", "CSF30", "CSF31", "CSF32")),]  ##includes plate-seq


# load homemade functions
{
  bigpar <-
    function(a = 1,
             b = 1,
             brewer.n = 8,
             brewer.name = "Dark2",
             cex.lab = 2,
             cex.main = 2,
             cex.axis = 1.5,
             mar = c(5.1, 5.1, 3.5, 2.1),
             mgp = c(3, 1, 0),
             ...) {
      par(
        mar = mar,
        mgp = mgp,
        cex.lab = cex.lab,
        cex.main = cex.main,
        cex.axis = cex.axis
      )
      par(mfrow = c(a, b), ...)
      palette(RColorBrewer::brewer.pal(brewer.n, brewer.name))
    }
  
  
  
  
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
  
  
}


# Clean up IDs to have CSFxxx as ID.
table(df$ID)
df$ID=gsub("CSFCSF","CSF",paste("CSF",df$ID,sep=""))

# Save ID / Diag information for network visualization
ForNetwork1=data.frame(Diagnosis=df$Diagnosis,ID=df$ID,stringsAsFactors = F)
ForNetwork1=ForNetwork1[which(duplicated(paste(ForNetwork1[,1],ForNetwork1[,2],sep="_"))==F),]
colnames(ForNetwork1)=c("x","y")


# preprocessing

# keep TCRs with both tcrA and b - i.e remove those with empty values
dim(df)
df=df[which(df$CDR3b!="" & df$CDR3a!=""),]
dim(df)

df=df[which(df$CDR3b!=""),] #for just tcrb analysis

# keep clones only (frequency > 1)
dim(df)
df=df[which(df$Frequency>1),]
dim(df)

# remove ambiguous TCRs - with ;
dim(df)
df=df[-grep(";",df$CDR3b),]
df=df[-grep(";",df$CDR3a),]
dim(df)

# create clonotype_id.long as unique ID for each clonotype
df$clonotype_id.long=paste(df$CDR3a,df$CDR3b,sep="_")
head(df)


# List of patients sharing TCR long
df[which(df$clonotype_id.long %in% df$clonotype_id.long[which(duplicated(df$clonotype_id)==T)]),]
df[which(df$CDR3b %in% df$CDR3b[which(duplicated(df$CDR3b)==T)]),]  ##identical TCRB

# export table of clonotype_id.long shared between patients
# toExport=df[which(df$clonotype_id.long %in% df$clonotype_id.long[which(duplicated(df$clonotype_id.long)==T)]),]
# write.csv(toExport,file = "/Users/benoitlehallier/BoxSync/bl/data/Aging/Human/Tcellreceptordata.Olivia.June19/combined_TCRs_for_Benoit.duplicated.clonotypes.csv",row.names = F)

#load package to calculate levenshteinSim
# install.packages("RecordLinkage")
require("RecordLinkage")
# calculate levenshtein distance between clonotype_ids
res.clonotype_id.long=NULL
for(i in 1:length(df$clonotype_id)){
  res.clonotype_id.long=rbind(res.clonotype_id.long,levenshteinSim(df$clonotype_id.long[i],df$clonotype_id.long))
  # print(paste(round((i/length(df$clonotype_id)*100),1),"%"))
}

res.cdr3b=NULL  #for tcrb analysis
for(i in 1:length(df$CDR3b)){
  res.cdr3b=rbind(res.cdr3b,levenshteinSim(df$CDR3b[i],df$CDR3b))
  # print(paste(round((i/length(df$clonotype_id)*100),1),"%"))
}

# rename rows/cols
colnames(res.clonotype_id.long)<-df$clonotype_id.long
rownames(res.clonotype_id.long)<-df$clonotype_id.long
head(res.clonotype_id.long[0:5, 0:5])

colnames(res.cdr3b)<-df$CDR3b
rownames(res.cdr3b)<-df$CDR3b
head(res.cdr3b[0:5, 0:5])

# Visualize levenshteinSim distribution
toHist=res.clonotype_id.long

toHist=res.cdr3b

# set NA to the upper triangle
toHist[upper.tri(toHist)]<-NA

# set NA to the diagonale
diag(toHist)<-NA


#pdf("histograms.TCR.analysis.pdf")
bigpar()
hist(toHist,xlab="Similarity between TCRs",main="")
hist(toHist,xlab="Similarity between TCRs",main="",xlim=c(.8,1),ylim = c(0,30),breaks = 100)
graphics.off()

# Make heatmap of clonotypes with levenshteinSim>0.8

res.clonotype_id.long.2=res.clonotype_id.long
res.clonotype_id.long.2[res.clonotype_id.long.2<0.8]<-0
diag(res.clonotype_id.long.2)<-0
res.clonotype_id.long.2=res.clonotype_id.long.2[which(apply(res.clonotype_id.long.2,1,function(x) sum(x>0,na.rm=T))>0),
                                                which(apply(res.clonotype_id.long.2,1,function(x) sum(x>0,na.rm=T))>0)]
toHeatmap=res.clonotype_id.long.2

res.cdr3b.2=res.cdr3b
res.cdr3b.2[res.cdr3b.2<1]<-0
diag(res.cdr3b.2)<-0
res.cdr3b.2=res.cdr3b.2[which(apply(res.cdr3b.2,1,function(x) sum(x>0,na.rm=T))>0),
                                                which(apply(res.cdr3b.2,1,function(x) sum(x>0,na.rm=T))>0)]
toHeatmap=res.cdr3b.2


pdf("tcr_similarity/final_plots/Heatmap.TCR.pdf")
require(gplots)
pairs.breaks <- seq(0, 1, by=0.01)
mycol <- colorpanel(n=length(pairs.breaks)-1,low="darkslateblue",mid="red",high="yellow")

par(oma=c(8,4,4,8))

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

graphics.off()


# network cutoff 0.8
{
  # extract link between clonotype_id.long and patients ID for network
  ForNetwork2=data.frame(df$ID,df$clonotype_id.long,stringsAsFactors = F)
  ForNetwork2=ForNetwork2[which(duplicated(paste(ForNetwork2[,1],ForNetwork2[,2],sep="_"))==F),]
  colnames(ForNetwork2)=c("x","y")
  head(ForNetwork2)
  # create table for qgraph
  
  # add TCRs similarity to the network
  ForNetwork3=res.clonotype_id.long
  ForNetwork3[upper.tri(ForNetwork3)]<-NA
  diag(ForNetwork3)<-NA
  head(ForNetwork3[,1:5])
  ForNetwork3[is.na(ForNetwork3)]<-0
  
  ForNetwork3=mat2edge(x = ForNetwork3,cutoff =  .8)
  dim(ForNetwork3)
  colnames(ForNetwork3)=c("x","y","weight")
  head(ForNetwork3)
  
  
  
  # network
  {
    toQgraph=rbind(ForNetwork1,ForNetwork2)
    toQgraph$weight=5
    toQgraph=rbind(toQgraph,ForNetwork3)
    head(toQgraph)
    
    # customize qgraph 
    {
      # color of Dx
      l=unique(c(toQgraph[,1],toQgraph[,2]))
      col.tmp=rep(adjustcolor("grey85",.3),length(l))
      size.tmp=rep(1,length(l))
      col.tmp[which(l %in% c("AD"))]<-adjustcolor("red",.8)
      col.tmp[which(l %in% c("HC"))]<-adjustcolor("blue",.8)
      col.tmp[which(l %in% c("MCI"))]<-adjustcolor("orange",.8)
      col.tmp[which(l %in% c("PD"))]<-adjustcolor("green4",.8)
      col.tmp[which(l %in% c("LBD"))]<-adjustcolor("purple",.8)
      
      # color of subjects
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="AD")]))]<-adjustcolor("red",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="HC")]))]<-adjustcolor("blue",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="MCI")]))]<-adjustcolor("orange",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="PD")]))]<-adjustcolor("green4",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="LBD")]))]<-adjustcolor("purple",.3)
      
      # size of nodes
      table(df$Diagnosis)
      size.tmp=rep(1,length(l))
      size.tmp[which(l %in% c("AD","HC","MCI","PD","LBD"))]<-100
      size.tmp[grep("CSF",l)]<-90
      table(size.tmp)
      
      # adjust size.tmp based on frequency
      #for(i in which(size.tmp==1)){
      #  size.tmp[i]=df$Frequency[which(df$clonotype_id.long==l[i])]
      }
      table(size.tmp)
      # Not that R will show warnings when 1 clone is shared by different subjects 
      # --> 2 warnings because 2 clones are shared between subjects
      # size of the dots is not accurate for these two clones
      size.tmp[which(size.tmp<50)]<-1
      size.tmp[which(size.tmp==90)]<-10
      size.tmp[which(size.tmp==100)]<-15
      
      
      # color of edges
      line.tmp=rep("black",nrow(toQgraph))
      # line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD","HC","MCI","PD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("HC"))]<-"blue"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("MCI"))]<-"orange"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("PD"))]<-"green4"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("LBD"))]<-"purple"
      
      labels.cex.tmp=l
      labels.cex.tmp=rep(0.00000000001,length(l))
      labels.cex.tmp[which(l %in% c("AD","HC","MCI","PD","LBD"))]<-1.5
      labels.cex.tmp[grep("CSF",l)]<-.8
    }
    
    require("qgraph")  
    pdf("tcr_similarity/final_plots/csf1-24_Lsim0.8_nodenetwork.pdf")
    toQgraph$weight<-1
    qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,labels=TRUE,label.cex=labels.cex.tmp,directed=F,layout="spring")
    graphics.off()
    
}
  


# network cutoff 0.9
{
  # extract link between clonotype_id.long and patients ID for network
  ForNetwork2=data.frame(df$ID,df$clonotype_id.long,stringsAsFactors = F)
  ForNetwork2=ForNetwork2[which(duplicated(paste(ForNetwork2[,1],ForNetwork2[,2],sep="_"))==F),]
  colnames(ForNetwork2)=c("x","y")
  head(ForNetwork2)
  # create table for qgraph
  
  # add TCRs similarity to the network
  ForNetwork3=res.clonotype_id.long
  ForNetwork3[upper.tri(ForNetwork3)]<-NA
  diag(ForNetwork3)<-NA
  head(ForNetwork3[,1:5])
  ForNetwork3[is.na(ForNetwork3)]<-0
  
  ForNetwork3=mat2edge(x = ForNetwork3,cutoff = .9)
  dim(ForNetwork3)
  colnames(ForNetwork3)=c("x","y","weight")
  head(ForNetwork3)

  
  # network
  {
    toQgraph=rbind(ForNetwork1,ForNetwork2)
    toQgraph$weight=5
    toQgraph=rbind(toQgraph,ForNetwork3)
    head(toQgraph)
    
    # customize qgraph 
    {
      # color of Dx
      l=unique(c(toQgraph[,1],toQgraph[,2]))
      col.tmp=rep(adjustcolor("grey85",.3),length(l))
      size.tmp=rep(1,length(l))
      col.tmp[which(l %in% c("AD"))]<-adjustcolor("red",.8)
      col.tmp[which(l %in% c("HC"))]<-adjustcolor("blue",.8)
      col.tmp[which(l %in% c("MCI"))]<-adjustcolor("orange",.8)
      col.tmp[which(l %in% c("PD"))]<-adjustcolor("green4",.8)
      col.tmp[which(l %in% c("LBD"))]<-adjustcolor("purple",.8)
      
      # color of subjects
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="AD")]))]<-adjustcolor("red",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="HC")]))]<-adjustcolor("blue",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="MCI")]))]<-adjustcolor("orange",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="PD")]))]<-adjustcolor("green4",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="LBD")]))]<-adjustcolor("purple",.3)
      
      # size of nodes
      table(df$Diagnosis)
      size.tmp=rep(1,length(l))
      size.tmp[which(l %in% c("AD","HC","MCI","PD","LBD"))]<-100
      size.tmp[grep("CSF",l)]<-90
      table(size.tmp)
      
      # adjust size.tmp based on frequency
      #for(i in which(size.tmp==1)){
      #  size.tmp[i]=df$Frequency[which(df$clonotype_id.long==l[i])]
      #}
      #table(size.tmp)
      # Not that R will show warnings when 1 clone is shared by different subjects 
      # --> 2 warnings because 2 clones are shared between subjects
      # size of the dots is not accurate for these two clones
      size.tmp[which(size.tmp<50)]<-1
      size.tmp[which(size.tmp==90)]<-10
      size.tmp[which(size.tmp==100)]<-15
      
      
      # color of edges
      line.tmp=rep("black",nrow(toQgraph))
      # line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD","HC","MCI","PD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("HC"))]<-"blue"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("MCI"))]<-"orange"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("PD"))]<-"green4"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("LBD"))]<-"purple"
      
      labels.cex.tmp=l
      labels.cex.tmp=rep(0.00000000001,length(l))
      labels.cex.tmp[which(l %in% c("AD","HC","MCI","PD","LBD"))]<-1.5
      labels.cex.tmp[grep("CSF",l)]<-.8
    }
    
    require("qgraph")  
    pdf("node.network.Threshold.09.pdf")
    toQgraph$weight<-1
    qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,labels=TRUE,label.cex=labels.cex.tmp,directed=F,layout="spring")
    graphics.off()
    
  }
  
}


# network cutoff 1.0
{
  # extract link between clonotype_id.long and patients ID for network
  ForNetwork2=data.frame(df$ID,df$clonotype_id.long,stringsAsFactors = F)
  ForNetwork2=ForNetwork2[which(duplicated(paste(ForNetwork2[,1],ForNetwork2[,2],sep="_"))==F),]
  colnames(ForNetwork2)=c("x","y")
  head(ForNetwork2)
  # create table for qgraph
  
  # add TCRs similarity to the network
  ForNetwork3=res.clonotype_id.long
  ForNetwork3[upper.tri(ForNetwork3)]<-NA
  diag(ForNetwork3)<-NA
  head(ForNetwork3[,1:5])
  ForNetwork3[is.na(ForNetwork3)]<-0
  
  ForNetwork3=mat2edge(x = ForNetwork3,cutoff =  .9999)
  dim(ForNetwork3)
  colnames(ForNetwork3)=c("x","y","weight")
  head(ForNetwork3)
  
  
  # network
  {
    toQgraph=rbind(ForNetwork1,ForNetwork2)
    toQgraph$weight=5
    toQgraph=rbind(toQgraph,ForNetwork3)
    head(toQgraph)
    
    # customize qgraph 
    {
      # color of Dx
      l=unique(c(toQgraph[,1],toQgraph[,2]))
      col.tmp=rep(adjustcolor("grey85",.3),length(l))
      size.tmp=rep(1,length(l))
      col.tmp[which(l %in% c("AD"))]<-adjustcolor("red",.8)
      col.tmp[which(l %in% c("HC"))]<-adjustcolor("blue",.8)
      col.tmp[which(l %in% c("MCI"))]<-adjustcolor("orange",.8)
      col.tmp[which(l %in% c("PD"))]<-adjustcolor("green4",.8)
      col.tmp[which(l %in% c("LBD"))]<-adjustcolor("purple",.8)
      
      # color of subjects
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="AD")]))]<-adjustcolor("red",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="HC")]))]<-adjustcolor("blue",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="MCI")]))]<-adjustcolor("orange",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="PD")]))]<-adjustcolor("green4",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="LBD")]))]<-adjustcolor("purple",.3)
      
      # size of nodes
      table(df$Diagnosis)
      size.tmp=rep(1,length(l))
      size.tmp[which(l %in% c("AD","HC","MCI","PD","LBD"))]<-100
      size.tmp[grep("CSF",l)]<-90
      table(size.tmp)
      
      # adjust size.tmp based on frequency
      #for(i in which(size.tmp==1)){
      #  size.tmp[i]=df$Frequency[which(df$clonotype_id.long==l[i])]
      #}
      #table(size.tmp)
      # Not that R will show warnings when 1 clone is shared by different subjects 
      # --> 2 warnings because 2 clones are shared between subjects
      # size of the dots is not accurate for these two clones
      size.tmp[which(size.tmp<50)]<-1
      size.tmp[which(size.tmp==90)]<-10
      size.tmp[which(size.tmp==100)]<-15
      
      
      # color of edges
      line.tmp=rep("black",nrow(toQgraph))
      # line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD","HC","MCI","PD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("HC"))]<-"blue"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("MCI"))]<-"orange"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("PD"))]<-"green4"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("LBD"))]<-"purple"
      
      labels.cex.tmp=l
      labels.cex.tmp=rep(0.00000000001,length(l))
      labels.cex.tmp[which(l %in% c("AD","HC","MCI","PD","LBD"))]<-1.5
      labels.cex.tmp[grep("CSF",l)]<-.8
    }
    
    require("qgraph")  
    pdf("tcr_similarity/node_network/csf1-24_Lsim1.0_nodenetwork.pdf")
    toQgraph$weight<-1
    qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,labels=TRUE,label.cex=labels.cex.tmp,directed=F,layout="spring")
    graphics.off()
    
  }
  
}



# network cutoff 1.0 cdr3b
{
  # extract link between clonotype_id.long and patients ID for network
  ForNetwork2=data.frame(df$ID,df$CDR3b,stringsAsFactors = F)
  ForNetwork2=ForNetwork2[which(duplicated(paste(ForNetwork2[,1],ForNetwork2[,2],sep="_"))==F),]
  colnames(ForNetwork2)=c("x","y")
  head(ForNetwork2)
  # create table for qgraph
  
  # add TCRs similarity to the network
  ForNetwork3=res.cdr3b
  ForNetwork3[upper.tri(ForNetwork3)]<-NA
  diag(ForNetwork3)<-NA
  head(ForNetwork3[,1:5])
  ForNetwork3[is.na(ForNetwork3)]<-0
  
  ForNetwork3=mat2edge(x = ForNetwork3,cutoff = .999)
  dim(ForNetwork3)
  colnames(ForNetwork3)=c("x","y","weight")
  head(ForNetwork3)
  
  
  
  # network
  {
    toQgraph=rbind(ForNetwork1,ForNetwork2)
    toQgraph$weight=5
    toQgraph=rbind(toQgraph,ForNetwork3)
    head(toQgraph)
    
    # customize qgraph 
    {
      # color of Dx
      l=unique(c(toQgraph[,1],toQgraph[,2]))
      col.tmp=rep(adjustcolor("grey85",.3),length(l))
      size.tmp=rep(1,length(l))
      col.tmp[which(l %in% c("AD"))]<-adjustcolor("red",.8)
      col.tmp[which(l %in% c("HC"))]<-adjustcolor("blue",.8)
      col.tmp[which(l %in% c("MCI"))]<-adjustcolor("orange",.8)
      col.tmp[which(l %in% c("PD"))]<-adjustcolor("green4",.8)
      col.tmp[which(l %in% c("LBD"))]<-adjustcolor("purple",.8)
      
      # color of subjects
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="AD")]))]<-adjustcolor("red",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="HC")]))]<-adjustcolor("blue",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="MCI")]))]<-adjustcolor("orange",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="PD")]))]<-adjustcolor("green4",.3)
      col.tmp[which(l %in% unique(df$ID[which(df$Diagnosis=="LBD")]))]<-adjustcolor("purple",.3)
      
      # size of nodes
      table(df$Diagnosis)
      size.tmp=rep(1,length(l))
      size.tmp[which(l %in% c("AD","HC","MCI","PD","LBD"))]<-100
      size.tmp[grep("CSF",l)]<-90
      table(size.tmp)
      
      # adjust size.tmp based on frequency
      for(i in which(size.tmp==1)){
        size.tmp[i]=df$Frequency[which(df$CDR3b==l[i])]
      }
      table(size.tmp)
      # Not that R will show warnings when 1 clone is shared by different subjects 
      # --> 2 warnings because 2 clones are shared between subjects
      # size of the dots is not accurate for these two clones
      size.tmp[which(size.tmp<50)]<-1
      size.tmp[which(size.tmp==90)]<-10
      size.tmp[which(size.tmp==100)]<-15
      
      
      # color of edges
      line.tmp=rep("black",nrow(toQgraph))
      # line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD","HC","MCI","PD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AD"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("HC"))]<-"blue"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("MCI"))]<-"orange"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("PD"))]<-"green4"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("LBD"))]<-"purple"
      
      labels.cex.tmp=l
      labels.cex.tmp=rep(0.00000000001,length(l))
      labels.cex.tmp[which(l %in% c("AD","HC","MCI","PD","LBD"))]<-1.5
      labels.cex.tmp[grep("CSF",l)]<-.8
    }
    
    require("qgraph")  
    pdf("tcr_similarity/final_plots/csf1-24_tcrb_Lsim1.0_nodenetwork.pdf")
    toQgraph$weight<-1
    qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,labels=TRUE,label.cex=labels.cex.tmp,directed=F,layout="spring")
    graphics.off()
    
  }
  
}
  

  