# convert raw counts to tpm
#please see use the following script as an example to convert raw counts to tpm
#https://mp.weixin.qq.com/sbiz=MzAxMDkxODM1Ng==&mid=2247508798&idx=1&sn=3fbbc1c173aa188db4e80de257ffb861&chksm=9b4be385ac3c6a9318426b1f2eef05550a330578d0a5060d520e185fbe027f0b961e6f947592&scene=178&cur_album_id=2035227969870168066#rd
#http://www.bio-info-trainee.com/956.html
#https://www.ensembl.org/info/data/ftp/index.html 

# set up 
library(rtracklayer)
library(AnnotationDbi)
library(GenomicFeatures)
library(tidyverse)
library(biomaRt)
library(dplyr)

# human
#export counts data 56602
counts<- read.csv("clin_counts.csv",row.names=1,check.names =F) #
#gene annotation
ensembl <-  useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="uswest.ensembl.org")
attributes <-  listAttributes(ensembl)

#'hgnc_symbol','entrezgene_id'
bioM.res<-getBM(attributes=c('ensembl_gene_id','hgnc_symbol','description','gene_biotype','entrezgene_id'), 
                filters = 'ensembl_gene_id', values = rownames(counts) , mart = ensembl)
dim(bioM.res) 

#combine raw counts data and gene annotation
counts.anno<-cbind(bioM.res[match(rownames(counts), bioM.res$ensembl_gene_id),],counts)
counts.anno[1:3,1:3]

#Human:  transcript lengths form UCSC exon
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#here we define the length of exons as the final gene length
if(!file.exists('gene_length.Rdata')){
  exon_txdb=exons(txdb)
  genes_txdb=genes(txdb)
  
  o = findOverlaps(exon_txdb,genes_txdb)
  t1=exon_txdb[queryHits(o)]
  t2=genes_txdb[subjectHits(o)]
  t1=as.data.frame(t1)
  t1$geneid=mcols(t2)[,1]
  g_l = lapply(split(t1,t1$geneid),function(x){
    # x=split(t1,t1$geneid)[[1]]
    head(x)
    tmp=apply(x,1,function(y){
      y[2]:y[3]
    })
    length(unique(unlist(tmp)))
  })
  head(g_l)
  g_l=data.frame(gene_id=names(g_l),length=as.numeric(g_l))
  
  save(g_l,file = 'gene_length_human.Rdata')
} 

load(file = 'gene_length_human.Rdata')

#gene_length 25750
head(g_l)
colnames(g_l) <- c("entrezgene_id","length")
gene_length <- g_l


#download from https://www.ensembl.org/info/data/ftp/index.html choose 
#Human>>Gene annotation>>Download GTF or GFF3>>File: Homo_sapiens.GRCh38.105.gtf.gz
#Mouse>>Gene annotation>>Download GTF or GFF3>>File: Mus_musculus.GRCm39.105.gtf.gz
#Human: transcript lengths for all the transcripts
#Read the genes from the GTP files
txdb <-makeTxDbFromGRanges(import("/Users/birongzhang/Downloads/reference/Homo_sapiens.GRCh38.105.gtf"))

#Export the transcript lengths for all the transcripts
allTranscripts <- transcriptLengths(txdb)

allGeneIDs <- unique(allTranscripts$gene_id)
allGeneLengths <- as.data.frame(allTranscripts %>%
                                  dplyr::group_by(gene_id) %>%
                                  dplyr::summarize(max.tx_len = max(tx_len)))


gene_length <- allGeneLengths[allGeneLengths$gene_id %in% rownames(counts),] 
colnames(gene_length) <- c("ensembl_gene_id","length")
gene_length[1:3,1:2] 


##combine count_anno and genelength 22922
counts_length <- merge(gene_length, counts.anno, by = "ensembl_gene_id")
counts_length[1:3,1:6]

#create a TPM matrix by dividing each column of the counts matrix by the estimatation of the gene length 
x <- counts_length[,colnames(counts)] /counts_length[,"length"]
#sum to 1 million.
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
#check it again
head(colSums(tpm.mat))
#tpm.mat <- tpm.mat %>% data.frame()  
dim(tpm.mat)

##add annoatation information
tpm.mat$ensembl_gene_id <- counts_length$ensembl_gene_id
tpm.mat$hgnc_symbol<- counts_length$hgnc_symbol
tpm.mat$entrezgene_id <- counts_length$entrezgene_id
write.csv(tpm.mat, "TPM_counts_gtf.csv")


# remove duplicated gene
tpm.mat2<-tpm.mat[!duplicated(tpm.mat$hgnc_symbol),] 
rownames(tpm.mat2) <- tpm.mat2$hgnc_symbol
tpm.mat2[1:3,1:3]
write.csv(tpm.mat, "TPM_counts_gtf_symbol.csv")





# mouse
#export raw counts data 56602
counts<- read.csv("clin_counts.csv",row.names=1,check.names =F) #
#gene annotation
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="uswest.ensembl.org")
#attributes = listAttributes(ensembl)

bioM.res<-getBM(attributes=c('ensembl_gene_id','mgi_symbol','gene_biotype', 'entrezgene_id'), 
                filters = 'ensembl_gene_id', values = rownames(counts) , mart = ensembl)
dim(bioM.res) 

#combine raw counts data and gene annotation
counts.anno<-cbind(bioM.res[match(rownames(counts), bioM.res$ensembl_gene_id),],counts)
counts.anno[1:3,1:3]

#Export the transcript lengths form UCSC
BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#here we define the length of exons as the final gene length
if(!file.exists('gene_length.Rdata')){
  exon_txdb=exons(txdb)
  genes_txdb=genes(txdb)
  
  o = findOverlaps(exon_txdb,genes_txdb)
  t1=exon_txdb[queryHits(o)]
  t2=genes_txdb[subjectHits(o)]
  t1=as.data.frame(t1)
  t1$geneid=mcols(t2)[,1]
  g_l = lapply(split(t1,t1$geneid),function(x){
    # x=split(t1,t1$geneid)[[1]]
    head(x)
    tmp=apply(x,1,function(y){
      y[2]:y[3]
    })
    length(unique(unlist(tmp)))
  })
  head(g_l)
  g_l=data.frame(gene_id=names(g_l),length=as.numeric(g_l))
  
  save(g_l,file = 'gene_length_mouse.Rdata')
} 

load(file = 'gene_length_mouse.Rdata')

#gene_length 24528
head(g_l)
colnames(g_l) <- c("entrezgene_id","length")
gene_length <- g_l


#Mouse
#Read the genes from the GTP files
txdb <-makeTxDbFromGRanges(import("/Users/birongzhang/Downloads/reference/Mus_musculus.GRCm39.105.gtf"))

#Export the transcript lengths for all the transcripts
allTranscripts <- transcriptLengths(txdb)

allGeneIDs <- unique(allTranscripts$gene_id)
allGeneLengths <- as.data.frame(allTranscripts %>%
                                  dplyr::group_by(gene_id) %>%
                                  dplyr::summarize(max.tx_len = max(tx_len)))


gene_length <- allGeneLengths[allGeneLengths$gene_id %in% rownames(counts),] 
colnames(gene_length) <- c("ensembl_gene_id","length")
gene_length[1:3,1:2] 



#combine count_anno and genelength 22922
counts_length <- merge(gene_length, counts.anno, by = "ensembl_gene_id")
counts_length[1:3,1:6]

#create a TPM matrix by dividing each column of the counts matrix by the estimatation of the gene length 
x <- counts_length[,colnames(counts)] /counts_length[,"length"]
#sum to 1 million.
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
#check it again
head(colSums(tpm.mat))
tpm.mat <- tpm.mat %>% data.frame()  
dim(tpm.mat)

#add annoatation information
tpm.mat$ensembl_gene_id <- counts_length$ensembl_gene_id
tpm.mat$mgi_symbol<- counts_length$mgi_symbol
tpm.mat$entrezgene_id <- counts_length$entrezgene_id
write.csv(tpm.mat, "TPM_counts_gtf.csv")

#remove duplicated gene
tpm.mat2<-tpm.mat[!duplicated(tpm.mat$mgi_symbol),] 
rownames(tpm.mat2) <- tpm.mat2$mgi_symbol
tpm.mat2[1:3,1:3]
write.csv(tpm.mat, "TPM_counts_gtf_symbol.csv")


# TPM filter
#Removing genes that are not expressed
table(rowSums(tpm.mat==0) == dim(tpm.mat)[2])
#Which values in tpm.mat are greater than 0
keep.exprs <- tpm.mat > 0
#This produces a logical matrix with TRUEs and FALSEs
head(keep.exprs)
#Summary of how many TRUEs there are in each row
table(rowSums(keep.exprs))
#we would like to keep genes that have at least 2 TRUES in each row of thresh
#keep <- rowSums(keep.exprs) >= round(dim(tpm.mat)[2]*0.1)
#keep <- rowSums(keep.exprs) >= 0
keep <- rowSums(keep.exprs) > 0
summary(keep)
#Subset the rows of countdata to keep the more highly expressed genes
tpm.mat.filtered <- tpm.mat[keep,]
rm(keep.exprs)
tpm.mat.filtered[1:4,1:4]

# TPM boxplot
boxplot(log2(tpm.mat.filtered+1),xlab="",ylab="Log2 transcript per million unnormalised",las=2,outline = F)
#Let's add a blue horizontal line that corresponds to the median log2TPM
abline(h=median(log2(tpm.mat.filtered+1)),col="blue")
title("Boxplots of log2TPM (unnormalised)")

# TPM normalised
BiocManager::install("preprocessCore")
library(preprocessCore)
tpm.mat.norm<- normalize.quantiles(log2(tpm.mat.filtered+1));
colnames(tpm.mat.norm)<- colnames(tpm.mat.filtered);
rownames(tpm.mat.norm)<- rownames(tpm.mat.filtered);
boxplot(tpm.mat.norm,xlab="",ylab="Log2 transcript per million normalised",las=2,outline = F)
abline(h=median(tpm.mat.norm),col="blue")
title("Boxplots of log2TPM (normalised)")

write.table(tpm.mat.filtered,file = "tpm.mat_filtered.txt",quote = FALSE,sep = "\t")
write.table(tpm.mat.norm,file = "tpm.mat_norm.txt",quote = FALSE,sep = "\t")
save(tpm.mat.norm,file = "tpm.mat_norm.Rdata")


# small tips
#1. for human and mouse, we need to change different database in useMart(), we can do this by the following command
??useMart
mart = useMart('ensembl')
listDatasets(mart)


#2.for human and mouse, we need to change different attributes in getBM()
attributes = listAttributes(ensembl)

#3. convert entrezgene_id to hgnc_symbol
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

?org.Hs.eg.db
g2s=toTable(org.Hs.egSYMBOL)
head(g2s)

genelength <-  merge(g_l,g2s,by='gene_id')
genelength <- genelength[,c(3,2)]
colnames(genelength) <- c("hgnc_symbol","length")

#4. remove duplicated
genelength<-genelength[!duplicated(genelength$hgnc_symbol),] 


