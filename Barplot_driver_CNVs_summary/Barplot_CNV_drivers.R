###Barplot summary for CNV-drivers in PDAC:

genes_amp=c("GATA6","MYC", "RIT1", "ZBTB7B", "AKT2","KRAS","ERBB2")
genes_del=c(cnvs_genes_del$Gene[cnvs_genes_del$Percentage>0.15],"PTEN")

arms_amp=c('1q')
arms_del=c('6p','6q','8p','9p','17p','17q','18p','18q')

order=c(genes_amp,arms_amp,genes_del,arms_del)

chr_map=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/hg38_gene_bands_sorted.txt',sep="\t",header=FALSE)
colnames(chr_map)[3:4]=c('arm','Gene')

chr_map$arm=gsub('([pq]).*','\\1',chr_map$arm)
chr_map$arm=paste(chr_map$V1, chr_map$arm,sep="")

chr_map_s=chr_map[(chr_map$arm %in% c(arms_amp, arms_del)) | (chr_map$Gene %in% c(genes_amp, genes_del)), ]

#Ampl_genes:
cnvs_amp=cnv[rownames(cnv) %in% genes_amp,]
cnvs_amp[is.na(cnvs_amp)]=0
cnvs_amp=ifelse(cnvs_amp>0.2,1,0)
cnvs_amp=as.data.frame(cnvs_amp)
cnvs_amp$Sum=rowSums(cnvs_amp)
cnvs_amp$Percentage=cnvs_amp$Sum/n_samples
cnvs_amp$Type="gene_amp"
cnvs_amp$Gene=rownames(cnvs_amp)
cnvs_genes_amp=cnvs_amp %>% select ('Gene','Percentage','Type')


#Del_genes:
cnvs_del=cnv[rownames(cnv) %in% genes_del,]
cnvs_del[is.na(cnvs_del)]=0
cnvs_del=ifelse(-cnvs_del>0.2,1,0)
cnvs_del=as.data.frame(cnvs_del)
cnvs_del$Sum=rowSums(cnvs_del)
cnvs_del$Percentage=cnvs_del$Sum/n_samples
cnvs_del$Type="gene_del"
cnvs_del$Gene=rownames(cnvs_del)
cnvs_genes_del=cnvs_del %>% select ('Gene','Percentage','Type')



for (arm in arms_amp){
	genes=chr_map$Gene[chr_map$arm==arm]
	cnv_s=cnv[rownames(cnv) %in% genes,]
	x=colMeans(cnv_s,na.rm=TRUE)
	x=x[!is.na(x)]
	perc=length(x[x>0.2])/n_samples
	cnvs_arms_amp=cbind(arm,perc,'arm_amp')
}
cnvs_arms_amp=as.data.frame(cnvs_arms_amp)
colnames(cnvs_arms_amp)=c('Gene','Percentage','Type')

cnvs_arms_del=NULL
for (arm in arms_del){
	genes=chr_map$Gene[chr_map$arm==arm]
	cnv_s=cnv[rownames(cnv) %in% genes,]
	x=colMeans(cnv_s,na.rm=TRUE)
	x=x[!is.na(x)]
	perc=length(x[-x>0.2])/n_samples
	st=cbind(arm,perc,'arm_del')
	cnvs_arms_del=rbind(cnvs_arms_del,st)
}
cnvs_arms_del=as.data.frame(cnvs_arms_del)
colnames(cnvs_arms_del)=c('Gene','Percentage','Type')

cnvs_arms=rbind(cnvs_arms_amp, cnvs_arms_del)
cnvs_genes=rbind(cnvs_genes_amp,cnvs_genes_del)

cnvs_arms$Percentage=as.numeric(as.character(unlist(cnvs_arms$Percentage)))

#all_cnvs=rbind(cnvs_arms,cnvs_genes)
all_cnvs=cnvs_genes
all_cnvs$Type_cnv=gsub('.*_(.*)','\\1',all_cnvs$Type)

all_cnvs=all_cnvs %>% mutate(Gene = factor(Gene, levels = rev(order))) %>% arrange(Gene)
all_cnvs$Percentage=all_cnvs$Percentage*100
###Making barplot (c("#253494","#8da0cb")):

p <- ggplot(all_cnvs, aes(fill=Type_cnv, y=Percentage, x=Gene))

p <- p + geom_bar(position="stack", stat="identity",width=0.5)

p <- p + scale_fill_manual(name="State", values=c("red","blue"), na.value=NA)+ scale_y_continuous(expand = c(0, 0), limits = c(0, 60))

p <- p + theme(axis.ticks = element_blank())+coord_flip() + theme_bw() + theme_classic()

pdf(paste("Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-17/CNVs/cnvs_summary_Drivers_added_2020-08-17.pdf",sep=""), width=5, height=4,useDingbats=FALSE)

print(p)

dev.off()    
