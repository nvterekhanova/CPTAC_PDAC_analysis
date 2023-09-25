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
genes=intersect(rownames(prot),colnames(cnvs_amp))
for (gene in genes){
#	samp_cnv=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]==1]
	
	cnv_tab=cnvs_amp %>% select ('Sample',gene)
	cnv_tab=cnv_tab[!is.na(cnv_tab[,2]),]
	cnv_tab[,2]=ifelse(cnv_tab[,2]>0.4,1,0)
	prot_tab=prot[rownames(prot)==gene,]
	prot_tab=t(prot_tab)
	prot_tab=as.data.frame(prot_tab)
	prot_tab$Sample=rownames(prot_tab)
	colnames(prot_tab)[1]=paste('Protein_',colnames(prot_tab)[1],sep='')
	res=merge(prot_tab, cnv_tab)
	x=cor.test(res[,2],res[,3],method="spearman")	
	
#		f_ch=mean(pr_cnv)-mean(pr_wt)
		st=cbind(gene,x$p.value,x$estimate)
		all_st=rbind(all_st,st)
}
all_st=as.data.frame(all_st)
all_st$V2=as.numeric(as.character(unlist(all_st$V2)))
all_st$FDR=p.adjust(all_st$V2,method='fdr')
all_st=all_st[order(all_st$FDR),]
#colnames(all_st)=c('Gene','Log2FoldChange','P_value','N_carriers','FDR')
#all_st$V2=round(all_st$V2,5)
all_st=all_st[all_st$FDR<0.05,]
all_st=all_st[!is.na(all_st$FDR),]
all_st$V3=as.numeric(as.character(unlist(all_st$V3)))
all_st=all_st[all_st$V3>0,]
prot_st=all_st


#####Gene expression analysis:
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
for (gene in genes){
	
	cnv_tab=cnvs_amp %>% select ('Sample',gene)
	cnv_tab=cnv_tab[!is.na(cnv_tab[,2]),]
	cnv_tab[,2]=ifelse(cnv_tab[,2]>0.4,1,0)
	prot_tab=prot[rownames(prot)==gene,]
	prot_tab=t(prot_tab)
	prot_tab=as.data.frame(prot_tab)
	prot_tab$Sample=rownames(prot_tab)
	colnames(prot_tab)[1]=paste('Protein_',colnames(prot_tab)[1],sep='')
	res=merge(prot_tab, cnv_tab)
	x=cor.test(res[,2],res[,3],method="spearman")	
	#		f_ch=mean(pr_cnv)-mean(pr_wt)
		st=cbind(gene,x$p.value,x$estimate)
		all_st=rbind(all_st,st)
}
all_st=as.data.frame(all_st)
all_st$V2=as.numeric(as.character(unlist(all_st$V2)))
all_st$FDR=p.adjust(all_st$V2,method='fdr')
all_st=all_st[order(all_st$FDR),]
#colnames(all_st)=c('Gene','Log2FoldChange','P_value','N_carriers','FDR')
#all_st$V2=round(all_st$V2,5)
all_st=all_st[all_st$FDR<0.05,]
all_st=all_st[!is.na(all_st$FDR),]
all_st$V3=as.numeric(as.character(unlist(all_st$V3)))
all_st=all_st[all_st$V3>0,]
gexpr_st=all_st

prot_st_sel=prot_st[prot_st$gene %in% gexpr_st$gene,]

###if spearman:
###gexpr:165,
###protein:23


###################################
###Making impact analysis:
###################################
genes_amp=c('XPO5', 'USP14','KRAS','NOTCH2','CDK6','POLB')


order=c(genes_amp)
cnvs_amp=cnv[rownames(cnv) %in% genes_amp,]
cnvs_amp[is.na(cnvs_amp)]=0
cnvs_amp=ifelse(cnvs_amp>0.4,1,0)
cnvs_amp=t(cnvs_amp)
cnvs_amp=as.data.frame(cnvs_amp)
cnvs_amp$Sample=rownames(cnvs_amp)


###Protein:
prot=read.table('Desktop/Projects/CPTAC3_PDAC/Data_Freeze_v0.2_07282020/proteomics_gene_level_MD_ratio_tumor.cct',header=TRUE,sep="\t")
rownames(prot)=prot[,1]
prot=prot[,-1]
table=prot

colnames(table) <- gsub('(C3[LN]).([0-9]+)','\\1-\\2',colnames(table))
table=table[,colnames(table) %in% samples]

prot=table

genes=c(genes_amp)

prot_s=prot[rownames(prot) %in% genes,]
prot_s=t(prot_s)
prot_s=scale(prot_s)
prot_s=t(prot_s)
prot_s=as.data.frame(prot_s)
#for(i in 1:nrow(prot_s)){
#	prot_s[i,]=ecdf_fun(prot_s[i,],prot_s[i,])
#}


prot_s$Gene=rownames(prot_s)
prot_s=melt(prot_s,by=list(prot_s$Gene))
colnames(prot_s)=c('Gene','Sample','Expr')
prot_s$Type="wt"

for (gene in genes_amp){
	samp_cnv=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]==1]
	prot_s$Type=ifelse(prot_s$Sample %in% samp_cnv & prot_s$Gene==gene,"cnv",prot_s$Type)
}

prot_s=prot_s %>% mutate(Type = factor(Type, levels = c('wt','cnv'))) %>% arrange(Type)
prot_s=prot_s %>% mutate(Gene = factor(Gene, levels = order)) %>% arrange(Gene)

protein_s=prot_s
###Making boxplots:

p <- ggplot(data = prot_s, aes(x=Gene, y=Expr))

p <- p + geom_boxplot(position="dodge",outlier.shape = NA,aes(fill=Type),width=0.9)
     
p <- p + scale_fill_manual(values = c("wt"="dark blue","cnv"="gold"))

p <- p + geom_point(aes(x=Gene, y=Expr,fill=Type),position=position_jitterdodge(0.5),shape=21,size=0.5)

p <- p + theme_bw() + theme_nogrid() + labs(title="",x="",y="Expression")+theme(axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))+theme(legend.position = "none")

p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))

p <- p + facet_grid( .~Gene, scales="free",space="free")


all_st=NULL
for (gene in c(genes_amp)){
	wt=prot_s$Expr[prot_s$Gene==gene & prot_s$Type=="wt"]
	cnv=prot_s$Expr[prot_s$Gene==gene & prot_s$Type=="cnv"]
if (length(wt)>1 & length(cnv)>1){
		x=wilcox.test(wt,cnv)
	st=cbind(gene,x$p.value)
	all_st=rbind(all_st,st)
}
}
all_st=as.data.frame(all_st)
all_st$V2=as.numeric(as.character(unlist(all_st$V2)))
all_st$FDR=p.adjust(all_st$V2,method='fdr')
all_st$V2=round(all_st$V2,5)

#P-vals:

    gene      V2         FDR
1   XPO5 0.00026 0.001552851	***
2  USP14 0.00125 0.003742990	***
3   KRAS 0.00304 0.003858813	***
4 NOTCH2 0.00304 0.003858813	***
5   CDK6 0.00360 0.003858813	***
6   POLB 0.00386 0.003858813	***





###Making same for the gene expression
gexpr=read.table("Desktop/Projects/CPTAC3_PDAC/Data_Freeze_v0.2_07282020/mRNA_RSEM_UQ_log2_Tumor.cct",header=TRUE,sep="\t")
table=gexpr
rownames(table)=table[,1]
table=table[,-1]

colnames(table) <- gsub('(C3[LN]).([0-9]+)','\\1-\\2',colnames(table))
table=table[,colnames(table) %in% samples]

prot=table

#genes=c(genes_amp, genes_del)

prot_s=prot[rownames(prot) %in% genes,]
prot_s=t(prot_s)
prot_s=scale(prot_s)
prot_s=t(prot_s)
prot_s=as.data.frame(prot_s)

#for(i in 1:nrow(prot_s)){
#	prot_s[i,]=ecdf_fun(prot_s[i,],prot_s[i,])
#}

prot_s$Gene=rownames(prot_s)
prot_s=melt(prot_s,by=list(prot_s$Gene))
colnames(prot_s)=c('Gene','Sample','Expr')
prot_s$Type="wt"

for (gene in genes_amp){
	samp_cnv=cnvs_amp$Sample[cnvs_amp[,colnames(cnvs_amp)==gene]==1]
	prot_s$Type=ifelse(prot_s$Sample %in% samp_cnv & prot_s$Gene==gene,"cnv",prot_s$Type)
}

prot_s=prot_s %>% mutate(Type = factor(Type, levels = c('wt','cnv'))) %>% arrange(Type)
prot_s=prot_s %>% mutate(Gene = factor(Gene, levels = order)) %>% arrange(Gene)

gexpr_s=prot_s

###Making boxplots:
p <- ggplot(data = prot_s, aes(x=Gene, y=Expr))

p <- p + geom_boxplot(position="dodge",outlier.shape = NA,aes(fill=Type),width=0.9)

p <- p + scale_fill_manual(values = c("wt"="dark blue","cnv"="gold"))

p <- p + geom_point(aes(x=Gene, y=Expr,fill=Type),position=position_jitterdodge(0.5),shape=21,size=0.5)

p <- p + theme_bw() + theme_nogrid() + labs(title="",x="",y="Expression")+theme(axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))+theme(legend.position = "none")

p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))

p <- p + facet_grid( .~Gene, scales="free",space="free")

#pdf(paste("Desktop/Projects/CPTAC3_PDAC/Analysis/CNVs/cnvs_gexpr_2020-08-04.pdf",sep=""), width=7, height=2.5,useDingbats=FALSE)
#print(p)
#dev.off()

all_st=NULL
for (gene in c(genes_amp)){
	wt=prot_s$Expr[prot_s$Gene==gene & prot_s$Type=="wt"]
	cnv=prot_s$Expr[prot_s$Gene==gene & prot_s$Type=="cnv"]
if (length(wt)>1 & length(cnv)>1){
	x=wilcox.test(wt,cnv)
	st=cbind(gene,x$p.value)
	all_st=rbind(all_st,st)
}
}
all_st=as.data.frame(all_st)
all_st$V2=as.numeric(as.character(unlist(all_st$V2)))
all_st$FDR=p.adjust(all_st$V2,method='fdr')
all_st$V2=round(all_st$V2,5)

#####FDR:
    gene      V2         FDR
1   XPO5 0.00018 0.001068329
2  USP14 0.00169 0.002535913
3   KRAS 0.00081 0.001611085
4 NOTCH2 0.00081 0.001611085
5   CDK6 0.00227 0.002729189
6   POLB 0.00362 0.003617935



###Combining all:
gexpr_s$Data="Gexpr"
protein_s$Data="Protein"
#phospho_s$Data="Phospho"

#all_s=rbind(gexpr_s,protein_s,phospho_s)
all_s=rbind(gexpr_s,protein_s)
###Making boxplots:
all_s=all_s %>% mutate(Type = factor(Type, levels = c('wt','cnv'))) %>% arrange(Type)
all_s=all_s %>% mutate(Data = factor(Data, levels = c('Gexpr','Protein','Phospho'))) %>% arrange(Gene)
all_s=all_s %>% mutate(Gene = factor(Gene, levels = order)) %>% arrange(Gene)


###Making violin:
p <- ggplot(data = all_s, aes(x=Gene, y=Expr))

p <- p + geom_violin(aes(fill=Type),trim=FALSE)

p <- p + scale_fill_manual(values = c("wt"="dark blue","cnv"="gold"))

p <- p + geom_point(aes(x=Gene, y=Expr,fill=Type),position=position_jitterdodge(0.3,dodge.width = 0.9),shape=21,size=0.5)

p <- p + theme_bw() + theme_nogrid() + labs(title="",x="",y="Expression")+theme(axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))+theme(legend.position = "none")

p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))

p <- p + facet_grid( Data~Gene, scales="free",space="free") +theme_minimal()


pdf(paste("Desktop/Projects/CPTAC3_PDAC/Analysis/2020-09-22/cnvs_allData.V2.2020-09-22.pdf",sep=""), width=10, height=4,useDingbats=FALSE)
print(p)
dev.off()
