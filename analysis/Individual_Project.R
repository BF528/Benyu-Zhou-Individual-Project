##cuflinks
cufflinks<-read.table("/projectnb2/bf528/students/bzhou7/P0_1_tophat/P0_1_cufflinks/genes.fpkm_tracking",header=TRUE,sep='\t')
#filter out fpkm>1 (one transcript per cell  in at least one out of six samples was used to obtain expressed transcripts). 
cufflinks<-subset(cufflinks,FPKM>0.01)
#Histogram of FPKM
hist(log10(cufflinks$FPKM),main = "Histogram of log10(FPKM)", xlab = "log10(FPKM)",breaks=20)
#Read the P0vsAd differential expression result
P0vsAd<-read.csv("/projectnb2/bf528/students/bzhou7/cuffdiff_out/gene_exp.diff",sep='\t')
#Sort the q-value 
temp<-order(P0vsAd$q_value)
#reorder the p4vsp7 dataframe by q vlaue
P0vsAd<-P0vsAd[temp,]
#Top 10 differentially expressed genes
top10<-P0vsAd[1:10,]
#Histogram of log2fold change
hist(P0vsAd$log2.fold_change.,breaks = 20,main = "Histogram of log2 fold change", xlab = "log2 fold change")
#subset sifnificant genes
P0vsAd.sub <- subset(P0vsAd,significant=="yes")
#Histogram of log2fold change for significant genes only
hist(P0vsAd.sub$log2.fold_change.,breaks = 20,main = "Histogram of log2 fold change", xlab = "log2 fold change")
#subset up and down regulated genes
up<-subset(P0vsAd.sub,log2.fold_change.>0)
down<-subset(P0vsAd.sub,log2.fold_change.<0)
#write out files
write(up$gene,file="up.csv",sep='\n')
write(down$gene,file="down.csv",sep='\n')
