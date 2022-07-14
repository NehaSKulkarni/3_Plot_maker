options(repos = list(CRAN="http://cran.rstudio.com/"))
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install(c("biomaRt","gridExtra"),update = F)
list.of.packages <- c(c("Hmisc","Tmisc","tidyverse","Rtsne"))
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
library(Hmisc)
library(biomaRt)
library(Tmisc)
library(Rtsne)
library(gridExtra)
library(tidyverse)

#Functions
compute_var_genes=function(dataset,my_mean=1,cv=0.5,png=F,name=NULL){
  dataset=dataset[apply(dataset, MARGIN = 1, function(x) any(x > 0)), ]
  
  means <- rowMeans(dataset)
  vars <- apply(dataset,1,var)
  cv2 <- vars/means^2
  Data=data.frame(means, log(cv2))
  colnames(Data)=c("means","cv")
  
  Data_ordered=Data[with(Data, order(means)), ]
  nls_fit <- nls(cv ~ a+b*(means^c), Data, start = list(a = 1,b=1,c=1))
  i <- order(Data$means)
  
  Data_ordered_filtered=Data_ordered[which(Data_ordered$means>my_mean & Data_ordered$cv>predict(nls_fit)[i]+cv),]
  list2=as.vector(rownames(Data_ordered_filtered))
  if(png==T) {
    png(paste0(name,"_cv",cv,"_mean_",my_mean,"_Var_genes.png",sep=""),width=4000,height=4000,res=250,pointsize=12)
    par(mar = c(6, 6, 4, 4) + 0.1,mgp=c(2,0.65,0),cex=0.8)
    plot(Data_ordered,cex=0.8,pch=21,bg="grey",col=ifelse((Data_ordered$means>my_mean)&(predict(nls_fit)[i]+cv<Data_ordered$cv),"red","blue"))
    lines(Data$means[i], predict(nls_fit)[i], col = "red")
    dev.off()
  } 
  return(dataset[list2,])
}

#make a tsne df from counts matrix
make_tsne <-function(counts,var_genes,perplexity=30){
  set.seed(100)
  tsne <- Rtsne(t(counts[var_genes,]),perplexity = perplexity)
  rownames(tsne$Y) <- colnames(counts) 
  tsne_df <- as.data.frame(tsne$Y)
  colnames(tsne_df) <- c("tsne1","tsne2") 
  tsne_df<-cbind(tsne_df,as.data.frame(t(counts)))
  
  print(identical(rownames(t(counts)),rownames(tsne_df)))
  
  return(tsne_df)
}

#plot markers on a tsne/pca - create with make_pca or make_tsne
plot_markers <- function(df,markers,return_grid=T,type=c("tsne","pca","umap")){
  if(type == "pca"){
    colnames(df)[1:2] <- c("PC1","PC2")
  }
  else if(type == "tsne"){
    colnames(df)[1:2] <- c("tSNE1","tSNE2")
  }
  
  markers <- markers[markers %in% colnames(df)]
  
  if(type == "pca"){
    plot_data_column = function (data, column){
      ggplot(df)+geom_point(aes_string("PC1","PC2",colour= column))
    }
  }
  else if(type == "tsne"){
    plot_data_column = function (data, column){
      ggplot(df)+geom_point(aes_string("tSNE1","tSNE2",colour= column))
    }
  }
  
  else if(type == "umap"){
    plot_data_column = function (data, column){
      ggplot(df)+geom_point(aes_string("UMAP1","UMAP2",colour= column))
    }
  }
  
  myplots <- lapply(markers, plot_data_column, data = pca_df)
  
  # if T return a grid instead of list
  if(return_grid){
    n <- length(myplots)
    nCol <- floor(sqrt(n))
    grid<-do.call("grid.arrange", c(myplots, ncol=2))
    return(grid)
  }
  
  #return list 
  return(myplots)
}

args <- commandArgs(trailingOnly = T)
markers <- unlist(strsplit(args,","))

#load in data
counts <- read.table("Counts_Lanner_Yan_Blackley.txt",header = T)
annotations <- read.csv("Sample_Info.csv",row.names = 1) 
colnames(annotations) <- capitalize(tolower(colnames(annotations)))
annotations$lineage<- sapply(annotations$Lineage,function(x){if( x == "epiblast"){"EPI"} else if(x == "primitive_endoderm"){"PrE"} else if(x == "trophectoderm"){"TE"} else {x}})
rownames(annotations) <- annotations$Sample_name
annotations <- annotations[,-1]
annotations$Timepoint<- sapply(rownames(annotations),function(x){strsplit(x,"_")[[1]][1]})
annotations$Embryo <- as.factor(annotations$Embryo)

# get gene names and mean transcript lengths
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = "www.ensembl.org"))
genes <- counts$EnsemblID
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","transcript_length"),values=genes,mart= mart)
lengths <- aggregate(G_list[,3], list(G_list$ensembl_gene_id), mean)
colnames(lengths) <- c("ensembl_gene_id","length")
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),1:2]
G_list <- merge(G_list,lengths,by = "ensembl_gene_id")
counts <- merge(G_list,counts,by.x = "ensembl_gene_id",by.y = "EnsemblID")
counts$hgnc_symbol <- apply(counts, 1,function(row){ if(row[2] == ""){ row[1]} else{row[2]}})
rownames(counts) <- make.unique(counts$hgnc_symbol)
lengths <- counts[,3,drop = F]
counts <- counts[,colnames(counts) %in% rownames(annotations)]

#log2FPKM
fpkm <- counts2fpkm(as.matrix(counts),lengths$length)
fpkm <- log2(fpkm+1)


markers <- c(colnames(annotations)[-which(colnames(annotations) %in% c("Embryo","Lineage"))],markers)

# plot entire dataset
var_genes <- rownames(compute_var_genes(fpkm))
tsne_df <- make_tsne(fpkm,var_genes,20)
tsne_df <- merge(tsne_df,annotations,by = "row.names")
rownames(tsne_df) <- tsne_df$Row.names
tsne_df <- tsne_df[,-1]
pdf("Lanner_Yan_Blackley_tSNE.pdf",height = length(markers)*2, width = 14)
plot_markers(tsne_df,markers)
dev.off()

# plot Lanner 
lanner_fpkm <- fpkm[,rownames(annotations[annotations$Dataset == "LANNER",])]
var_genes <- rownames(compute_var_genes(lanner_fpkm))
tsne_df <- make_tsne(lanner_fpkm,var_genes,20)
tsne_df <- merge(tsne_df,annotations,by = "row.names")
rownames(tsne_df) <- tsne_df$Row.names
tsne_df <- tsne_df[,-1]
pdf("Lanner_tSNE.pdf",height = length(markers)*2, width = 14)
plot_markers(tsne_df,markers,type = "tsne")
dev.off()

# plot Yan 
yan_fpkm <- fpkm[,rownames(annotations[annotations$Dataset == "YAN",])]
var_genes <- rownames(compute_var_genes(yan_fpkm))
tsne_df <- make_tsne(yan_fpkm,var_genes,10)
tsne_df <- merge(tsne_df,annotations,by = "row.names")
rownames(tsne_df) <- tsne_df$Row.names
tsne_df <- tsne_df[,-1]
pdf("Yan_tSNE.pdf",height = length(markers)*2, width = 14)
plot_markers(tsne_df,markers,type = "tsne")
dev.off()

# plot Niakan 
niakan_fpkm <- fpkm[,rownames(annotations[annotations$Dataset == "NIAKAN",])]
var_genes <- rownames(compute_var_genes(niakan_fpkm))
tsne_df <- make_tsne(niakan_fpkm,var_genes,5)
tsne_df <- merge(tsne_df,annotations,by = "row.names")
rownames(tsne_df) <- tsne_df$Row.names
tsne_df <- tsne_df[,-1]
pdf("Niakan_tSNE.pdf",height = length(markers)*2, width = 14)
plot_markers(tsne_df,markers,type = "tsne")
dev.off()


