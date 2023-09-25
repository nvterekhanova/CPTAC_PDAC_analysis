#library(copynumber)
library(maftools)
library(bnstruct)
library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(plyr)

library(ggplot2)
library(ggrepel)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)


file=read.table('Desktop/Projects/CPTAC3_PDAC/Data_Freeze_1.0_08262020/SCNA_log2_segment_level.cct',header=TRUE,sep="\t")


meta=read_delim('Desktop/Projects/CPTAC3_PDAC/Data_Freeze_v0.2_07282020/meta_table_clinic.tsv',delim="\t")
meta=as.data.frame(meta)
samples=meta$case_id[!(meta$histology_diagnosis=="Adenosquamous carcinoma") &meta$tumor_included_for_the_study=="yes"]

####Using only High Quality samples, upd 2020-08-17:
samples_hq=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-17/Samples_HQ.txt',sep="\t",header=FALSE)
samples=intersect(samples, samples_hq$V1)
n_samples=length(samples)

file=file[file$Sample %in% samples,]
file=file[!(file$Chromosome %in% c("chrX","chrY")),]

write.table(file,'Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/Filtered_samples_SCNA_log2_segment_level_withoutXY.cct',sep="\t",quote=FALSE,row.names=FALSE)



all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- system.file("extdata", "scores.gistic", package = "maftools")

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)

gisticChromPlot(gistic = laml.gistic, markBands = "all")



dir="Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/271838/"
all.lesions=paste(dir,"all_lesions.conf_95.txt",sep="")
amp.genes=paste(dir,"amp_genes.conf_95.txt",sep="")
del.genes=paste(dir,"del_genes.conf_95.txt",sep="")
scores.gis=paste(dir,"scores.gistic",sep="")

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = FALSE)

gisticChromPlot(gistic = laml.gistic, markBands = "all")

pdf(paste("Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/cnvs_GISTIC_2020-08-29.pdf",sep=""), width=9, height=5,useDingbats=FALSE)
gisticChromPlot(gistic = laml.gistic, markBands = "all")
dev.off()   



gisticOncoPlot(gistic = laml.gistic, sortByAnnotation = TRUE, top = 10)


####Making new volcano-plots:

###CNV-summary plot:
library(bnstruct)
library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(plyr)

library(ggplot2)
library(ggrepel)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)

meta=read_delim('Desktop/Projects/CPTAC3_PDAC/Data_Freeze_v0.2_07282020/meta_table_clinic.tsv',delim="\t")
meta=as.data.frame(meta)
samples=meta$case_id[!(meta$histology_diagnosis=="Adenosquamous carcinoma") &meta$tumor_included_for_the_study=="yes"]

####Using only High Quality samples, upd 2020-08-17:
samples_hq=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-17/Samples_HQ.txt',sep="\t",header=FALSE)
samples=intersect(samples, samples_hq$V1)
n_samples=length(samples)
####################################

genes_amp_gistic=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/genes_amp_gistic.txt',sep="\t",header=FALSE)
genes_del_gistic=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/genes_del_gistic.txt',sep="\t",header=FALSE)

genes_amp_gistic=as.character(genes_amp_gistic$V1)
genes_del_gistic=as.character(genes_del_gistic$V1)

cnv=read_delim('Desktop/Projects/CPTAC3_PDAC/Data_Freeze_v0.2_07282020/PDA_WGS_CNA_gene_level_rerun_v1.cct',delim="\t")
cnv=as.data.frame(cnv)
#ignore the warning
rownames(cnv)=cnv[,1]
cnv=cnv[,-1]
cnv=(cnv[,colnames(cnv) %in% samples])

###CNV-impacts on expression
cnvs_amp=cnv
cnvs_amp=t(cnvs_amp)
cnvs_amp=as.data.frame(cnvs_amp)
cnvs_amp$Sample=rownames(cnvs_amp)

cnvs_del=cnv
cnvs_del=t(cnvs_del)
cnvs_del=as.data.frame(cnvs_del)
cnvs_del$Sample=rownames(cnvs_del)

###Protein:
prot=read.table('Desktop/Projects/CPTAC3_PDAC/Data_Freeze_v0.2_07282020/proteomics_gene_level_MD_ratio_tumor.cct',header=TRUE,sep="\t")
rownames(prot)=prot[,1]
prot=prot[,-1]
table=prot

colnames(table) <- gsub('(C3[LN]).([0-9]+)','\\1-\\2',colnames(table))
table=table[,colnames(table) %in% samples]

prot=table
####Looking only at genes from the ampl. regions:
prot_amp=prot[rownames(prot) %in% genes_amp_gistic,]
prot_del=prot[rownames(prot) %in% genes_del_gistic,]
prot_all=prot

all_st=NULL
prot=prot_amp
for (i in 1:nrow(prot)){
	gene=rownames(prot)[i]
#	samp_cnv=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]==1]
	
	samp_cnv=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]>0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	wt=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]<0.4 & cnvs_amp[,colnames(cnvs_amp)==gene]>-0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	pr_cnv=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% samp_cnv])))
	pr_wt=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% wt])))
	pr_cnv=pr_cnv[!is.na(pr_cnv)]
	pr_wt=pr_wt[!is.na(pr_wt)]
	if (length(pr_cnv)>4 & length(pr_wt)>4){
		x=wilcox.test(pr_cnv,pr_wt,alternative='greater')
		f_ch=mean(pr_cnv)-mean(pr_wt)
		st=cbind(gene,f_ch,x$p.value,length(pr_cnv))
		all_st=rbind(all_st,st)
	}
}
all_st=as.data.frame(all_st)
all_st$V3=as.numeric(as.character(unlist(all_st$V3)))
all_st$f_ch=as.numeric(as.character(unlist(all_st$f_ch)))
all_st$FDR=p.adjust(all_st$V3,method='fdr')
all_st=all_st[order(all_st$FDR),]
colnames(all_st)=c('Gene','Log2FoldChange','P_value','N_carriers','FDR')
#all_st$V2=round(all_st$V2,5)
write.table(all_st,"Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/protein_ampl_impact.tsv",sep='\t',quote=FALSE,row.names=FALSE)
g_prot=as.character(all_st$Gene[all_st$FDR<0.05])

all_st=NULL
prot=prot_del
for (i in 1:nrow(prot)){
	gene=rownames(prot)[i]
	
	cnvs_amp=cnvs_del
	samp_cnv=cnvs_amp$Sample[-cnvs_amp[,colnames(cnvs_amp)==gene]>0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	wt=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]<0.4 & cnvs_amp[,colnames(cnvs_amp)==gene]>-0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	pr_cnv=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% samp_cnv])))
	pr_wt=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% wt])))
	pr_cnv=pr_cnv[!is.na(pr_cnv)]
	pr_wt=pr_wt[!is.na(pr_wt)]
	if (length(pr_cnv)>4 & length(pr_wt)>4){
		x=wilcox.test(pr_cnv,pr_wt,alternative='less')
		f_ch=mean(pr_cnv)-mean(pr_wt)
		st=cbind(gene,f_ch,x$p.value,length(pr_cnv))
		all_st=rbind(all_st,st)
	}
}
all_st=as.data.frame(all_st)
all_st$V3=as.numeric(as.character(unlist(all_st$V3)))
all_st$f_ch=as.numeric(as.character(unlist(all_st$f_ch)))
all_st$FDR=p.adjust(all_st$V3,method='fdr')
all_st=all_st[order(all_st$FDR),]
colnames(all_st)=c('Gene','Log2FoldChange','P_value','N_carriers','FDR')
#all_st$V2=round(all_st$V2,5)
write.table(all_st,"Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/protein_del_impact.tsv",sep='\t',quote=FALSE,row.names=FALSE)


###Making Volcanos:
library(ggrepel)
library(grid)
library(png)
library(gridExtra)
library(ggrastr)

del=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/protein_del_impact.tsv',sep='\t',header=TRUE)
amp=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/protein_ampl_impact.tsv',sep='\t',header=TRUE)
del$Type="Del"
amp$Type="Amp"
all_w_test_stat=rbind(del,amp)

all_fdr=all_w_test_stat[all_w_test_stat$FDR<0.05,]
results=all_w_test_stat

f_ch_cut=0.7


all_subs_fdr=results[(results$FDR<0.05),]

all_subs_fdr=all_subs_fdr[all_subs_fdr$Gene %in% g_16,]

results_g=results
results_g=results_g[(results_g$Log2FoldChange>0 & results_g$Type=="Amp") | (results_g$Log2FoldChange<0 & results_g$Type=="Del"),]
res_subs_fdr_gr=all_subs_fdr

results_g=results_g[(results_g$Log2FoldChange>0 & results_g$Type=="Amp") | (results_g$Log2FoldChange<0 & results_g$Type=="Del"),]

p <- ggplot(res_subs_fdr_gr, aes(y=-log10(FDR),x= Log2FoldChange)) + geom_point_rast(alpha=0.5, aes(colour=Type), data = results_g[results_g$FDR<0.05,])

p <- p + geom_point_rast(alpha=0.5, colour="grey", data = results_g[results_g$FDR>=0.05,])

p <- p + theme_minimal()  + geom_vline(xintercept=0, alpha=0.5,size=0.2)+ geom_point(data = res_subs_fdr_gr, aes(col = Type), alpha=0.95)

p <- p + ylab("-log10(FDR)")

p <- p + xlab("log2(Protein Expression Fold Change)") 

p <- p + geom_text_repel(aes(y=-log10(FDR),x= Log2FoldChange,label=ifelse(!is.na(Gene), as.character(Gene),NA)),alpha=1,size=3,segment.size  = 0.2,color = "black")

p <- p + scale_color_manual(values=c("red","blue"))

p <- p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p <- p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank(),legend.position="none")+ xlim(-1.5,1.5)

p <- p + geom_vline(xintercept=f_ch_cut, size=0.4,alpha=0.7,linetype=2) + geom_vline(xintercept=-f_ch_cut, size=0.4,alpha=0.7,linetype=2)+geom_hline(yintercept=-log10(0.05), alpha=0.7,size=0.4,linetype=2) + geom_hline(yintercept=0, alpha=0.5,size=0.2) +ylim(0,6)

pdf(paste("Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/Protein_volcano.pdf",sep=""), width=6, height=5,useDingbats=FALSE)
print(p)
dev.off()


#######################################
###Making the same for gene expression:
#######################################

gexpr=read.table("Desktop/Projects/CPTAC3_PDAC/Data_Freeze_v0.2_07282020/mRNA_RSEM_UQ_log2_Tumor.cct",header=TRUE,sep="\t")
table=gexpr
rownames(table)=table[,1]
table=table[,-1]

colnames(table) <- gsub('(C3[LN]).([0-9]+)','\\1-\\2',colnames(table))
table=table[,colnames(table) %in% samples]

prot=table
prot_amp=prot[rownames(prot) %in% genes_amp_gistic,]
prot_del=prot[rownames(prot) %in% genes_del_gistic,]
prot_all=prot


################
#Amplifications:
################
all_st=NULL
prot=prot_amp
for (i in 1:nrow(prot)){
	gene=rownames(prot)[i]
	
	samp_cnv=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]>0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	wt=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]<0.4 & cnvs_amp[,colnames(cnvs_amp)==gene]>-0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	pr_cnv=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% samp_cnv])))
	pr_wt=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% wt])))
	pr_cnv=pr_cnv[!is.na(pr_cnv)]
	pr_wt=pr_wt[!is.na(pr_wt)]
	if (length(pr_cnv)>4 & length(pr_wt)>4){
		x=wilcox.test(pr_cnv,pr_wt,alternative='greater')
		f_ch=mean(pr_cnv)-mean(pr_wt)
		st=cbind(gene,f_ch,x$p.value,length(pr_cnv))
		all_st=rbind(all_st,st)
	}
}
all_st=as.data.frame(all_st)
all_st$V3=as.numeric(as.character(unlist(all_st$V3)))
all_st$f_ch=as.numeric(as.character(unlist(all_st$f_ch)))
all_st$FDR=p.adjust(all_st$V3,method='fdr')
all_st=all_st[order(all_st$FDR),]
colnames(all_st)=c('Gene','Log2FoldChange','P_value','N_carriers','FDR')
write.table(all_st,"Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/gene_expression_ampl_impact.tsv",sep='\t',quote=FALSE,row.names=FALSE)
g_expr=as.character(all_st$Gene[all_st$FDR<0.05])

###########
#Deletions:
###########
all_st=NULL
prot=prot_del
for (i in 1:nrow(prot)){
	gene=rownames(prot)[i]
	
	cnvs_amp=cnvs_del
	samp_cnv=cnvs_amp$Sample[-cnvs_amp[,colnames(cnvs_amp)==gene]>0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	wt=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]<0.4 & cnvs_amp[,colnames(cnvs_amp)==gene]>-0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	pr_cnv=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% samp_cnv])))
	pr_wt=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% wt])))
	pr_cnv=pr_cnv[!is.na(pr_cnv)]
	pr_wt=pr_wt[!is.na(pr_wt)]
	if (length(pr_cnv)>4 & length(pr_wt)>4){
		x=wilcox.test(pr_cnv,pr_wt,alternative='less')
		f_ch=mean(pr_cnv)-mean(pr_wt)
		st=cbind(gene,f_ch,x$p.value,length(pr_cnv))
		all_st=rbind(all_st,st)
	}
}
all_st=as.data.frame(all_st)
all_st$V3=as.numeric(as.character(unlist(all_st$V3)))
all_st$f_ch=as.numeric(as.character(unlist(all_st$f_ch)))
all_st$FDR=p.adjust(all_st$V3,method='fdr')
all_st=all_st[order(all_st$FDR),]
colnames(all_st)=c('Gene','Log2FoldChange','P_value','N_carriers','FDR')
#all_st$V2=round(all_st$V2,5)
write.table(all_st,"Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/gene_expression_del_impact.tsv",sep='\t',quote=FALSE,row.names=FALSE)

################
###volcano-plot:
################
del=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/gene_expression_del_impact.tsv',sep='\t',header=TRUE)
amp=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/gene_expression_ampl_impact.tsv',sep='\t',header=TRUE)
del$Type="Del"
amp$Type="Amp"
all_w_test_stat=rbind(del,amp)


all_fdr=all_w_test_stat[all_w_test_stat$FDR<0.05,]
results=all_w_test_stat


f_ch_cut=2

all_subs_fdr=results[(results$FDR<0.05),]

all_subs_fdr=all_subs_fdr[all_subs_fdr$Gene %in% g_16,]

results_g=results
results_g=results_g[(results_g$Log2FoldChange>0 & results_g$Type=="Amp") | (results_g$Log2FoldChange<0 & results_g$Type=="Del"),]
res_subs_fdr_gr=all_subs_fdr

results_g=results_g[(results_g$Log2FoldChange>0 & results_g$Type=="Amp") | (results_g$Log2FoldChange<0 & results_g$Type=="Del"),]

p <- ggplot(res_subs_fdr_gr, aes(y=-log10(FDR),x= Log2FoldChange)) + geom_point_rast(alpha=0.5, aes(colour=Type), data = results_g[results_g$FDR<0.05,])

p <- p + geom_point_rast(alpha=0.5, colour="grey", data = results_g[results_g$FDR>=0.05,])

p <- p + theme_minimal()  + geom_vline(xintercept=0, alpha=0.5,size=0.2)+ geom_point(data = res_subs_fdr_gr, aes(col = Type), alpha=0.95)

p <- p + ylab("-log10(FDR)")

p <- p + xlab("log2(Gene Expression Fold Change)") 

p <- p + geom_text_repel(aes(y=-log10(FDR),x= Log2FoldChange,label=ifelse(!is.na(Gene), as.character(Gene),NA)),alpha=1,size=3,segment.size  = 0.2,color = "black")

p <- p + scale_color_manual(values=c("red","blue"))

p <- p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p <- p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank(),legend.position="none")+xlim(-3.7,3.7)

p <- p + geom_vline(xintercept=f_ch_cut, size=0.4,alpha=0.7,linetype=2) + geom_vline(xintercept=-f_ch_cut, size=0.4,alpha=0.7,linetype=2)+geom_hline(yintercept=-log10(0.05), alpha=0.7,size=0.4,linetype=2) + geom_hline(yintercept=0, alpha=0.5,size=0.2) +ylim(0,6)

pdf(paste("Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/Gene_expression_volcano.pdf",sep=""), width=6, height=5,useDingbats=FALSE)
print(p)
dev.off()

#####

g_16=intersect(g_prot,g_expr)

####Intersected between the Gene expression and Protein:
#[1] "SUPT5H"  "MEPCE"   "XPO5"    "GBA2"    "AGFG2"   "AKAP9"   "GTPBP2"  "MRPS18A" "SAMD4B"  "POLR1C" 
#[11] "MED29"   "HINT2"   "TSC22D4" "DSC2"    "TLN1"    "LRCH4" 
SUPT5H
MEPCE
XPO5
GBA2
AGFG2
AKAP9
GTPBP2
MRPS18A
SAMD4B
POLR1C
MED29
HINT2
TSC22D4
DSC2
TLN1
LRCH4

##################
###Phosphoprotein:
##################
phospho=read.table("Desktop/Projects/CPTAC3_PDAC/Data_Freeze_v0.2_07282020/phosphoproteomics_gene_level_MD_ratio_tumor.cct",header=TRUE,sep="\t")
table=phospho
rownames(table)=table[,1]
table=table[,-1]

colnames(table) <- gsub('(C3[LN]).([0-9]+)','\\1-\\2',colnames(table))
table=table[,colnames(table) %in% samples]

prot=table
prot_amp=prot[rownames(prot) %in% genes_amp_gistic,]
prot_del=prot[rownames(prot) %in% genes_del_gistic,]
prot_all=prot

################
#Amplifications:
################
all_st=NULL
prot=prot_amp
for (i in 1:nrow(prot)){
	gene=rownames(prot)[i]
	
	samp_cnv=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]>0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	wt=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]<0.4 & cnvs_amp[,colnames(cnvs_amp)==gene]>-0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	pr_cnv=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% samp_cnv])))
	pr_wt=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% wt])))
	pr_cnv=pr_cnv[!is.na(pr_cnv)]
	pr_wt=pr_wt[!is.na(pr_wt)]
	if (length(pr_cnv)>4 & length(pr_wt)>4){
		x=wilcox.test(pr_cnv,pr_wt,alternative='greater')
		f_ch=mean(pr_cnv)-mean(pr_wt)
		st=cbind(gene,f_ch,x$p.value,length(pr_cnv))
		all_st=rbind(all_st,st)
	}
}
all_st=as.data.frame(all_st)
all_st$V3=as.numeric(as.character(unlist(all_st$V3)))
all_st$f_ch=as.numeric(as.character(unlist(all_st$f_ch)))
all_st$FDR=p.adjust(all_st$V3,method='fdr')
all_st=all_st[order(all_st$FDR),]
colnames(all_st)=c('Gene','Log2FoldChange','P_value','N_carriers','FDR')
write.table(all_st,"Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/phosphoprotein_ampl_impact.tsv",sep='\t',quote=FALSE,row.names=FALSE)

###########
#Deletions:
###########
all_st=NULL
prot=prot_del
for (i in 1:nrow(prot)){
	gene=rownames(prot)[i]
	
	cnvs_amp=cnvs_del
	samp_cnv=cnvs_amp$Sample[-cnvs_amp[,colnames(cnvs_amp)==gene]>0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	wt=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]<0.4 & cnvs_amp[,colnames(cnvs_amp)==gene]>-0.4 & !is.na(cnvs_amp[,colnames(cnvs_amp)==gene])]
	pr_cnv=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% samp_cnv])))
	pr_wt=as.numeric(as.character(unlist(prot[i,colnames(prot) %in% wt])))
	pr_cnv=pr_cnv[!is.na(pr_cnv)]
	pr_wt=pr_wt[!is.na(pr_wt)]
	if (length(pr_cnv)>4 & length(pr_wt)>4){
		x=wilcox.test(pr_cnv,pr_wt,alternative='less')
		f_ch=mean(pr_cnv)-mean(pr_wt)
		st=cbind(gene,f_ch,x$p.value,length(pr_cnv))
		all_st=rbind(all_st,st)
	}
}
all_st=as.data.frame(all_st)
all_st$V3=as.numeric(as.character(unlist(all_st$V3)))
all_st$f_ch=as.numeric(as.character(unlist(all_st$f_ch)))
all_st$FDR=p.adjust(all_st$V3,method='fdr')
all_st=all_st[order(all_st$FDR),]
colnames(all_st)=c('Gene','Log2FoldChange','P_value','N_carriers','FDR')
#all_st$V2=round(all_st$V2,5)
write.table(all_st,"Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/phosphoprotein_del_impact.tsv",sep='\t',quote=FALSE,row.names=FALSE)

################
###volcano-plot:
################
del=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/phosphoprotein_del_impact.tsv',sep='\t',header=TRUE)
amp=read.table('Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/phosphoprotein_ampl_impact.tsv',sep='\t',header=TRUE)
del$Type="Del"
amp$Type="Amp"
all_w_test_stat=rbind(del,amp)


all_fdr=all_w_test_stat[all_w_test_stat$FDR<0.05,]
results=all_w_test_stat

f_ch_cut=0.5

#all_subs_fdr=results[((results$Gene %in% c('GATA6','MYC','AKT2','KRAS','ERBB2','CDKN2A','SMAD4','ARID1A','PTEN')) & results$FDR<0.05) | (abs(results$Log2FoldChange)>f_ch_cut & results$FDR<0.05),]

all_subs_fdr=results[(abs(results$Log2FoldChange)>f_ch_cut & results$FDR<0.05),]

all_subs_fdr=all_subs_fdr[(all_subs_fdr$Log2FoldChange>0 & all_subs_fdr$Type=="Amp") | (all_subs_fdr$Log2FoldChange<0 & all_subs_fdr$Type=="Del"),]


results_g=results
results_g=results_g[(results_g$Log2FoldChange>0 & results_g$Type=="Amp") | (results_g$Log2FoldChange<0 & results_g$Type=="Del"),]
res_subs_fdr_gr=all_subs_fdr

results_g=results_g[(results_g$Log2FoldChange>0 & results_g$Type=="Amp") | (results_g$Log2FoldChange<0 & results_g$Type=="Del"),]

p <- ggplot(res_subs_fdr_gr, aes(y=-log10(FDR),x= Log2FoldChange)) + geom_point_rast(alpha=0.5, aes(colour=Type), data = results_g[results_g$FDR<0.05,])

p <- p + geom_point_rast(alpha=0.5, colour="grey", data = results_g[results_g$FDR>=0.05,])

p <- p + theme_minimal()  + geom_vline(xintercept=0, alpha=0.5,size=0.2)+ geom_point(data = res_subs_fdr_gr, aes(col = Type), alpha=0.95)

p <- p + ylab("-log10(FDR)")

p <- p + xlab("log2(Phosphoprotein Expression Fold Change)") 

p <- p + geom_text_repel(aes(y=-log10(FDR),x= Log2FoldChange,label=ifelse(!is.na(Gene), as.character(Gene),NA)),alpha=1,size=3,segment.size  = 0.2,color = "black")

p <- p + scale_color_manual(values=c("red","blue"))

p <- p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p <- p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank(),legend.position="none")+ xlim(-1.5,1.5)

p <- p + geom_vline(xintercept=f_ch_cut, size=0.4,alpha=0.7,linetype=2) + geom_vline(xintercept=-f_ch_cut, size=0.4,alpha=0.7,linetype=2)

p <- p + geom_hline(yintercept=-log10(0.05), alpha=0.7,size=0.4,linetype=2) + geom_hline(yintercept=0, alpha=0.5,size=0.2)

pdf(paste("Desktop/Projects/CPTAC3_PDAC/Analysis/2020-08-29/CNV_impact/Phosphoprotein_volcano.pdf",sep=""), width=6, height=5,useDingbats=FALSE)
print(p)
dev.off()
