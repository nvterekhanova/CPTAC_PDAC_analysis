###Making barplot (c("#253494","#8da0cb")):
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

bp=read_delim('/Users/nadezhdaterekhanova/Desktop/Projects/Additional/CPTAC3_PDAC/Analysis/2020-09-11/GO_terms_enriched.txt',delim="\t",show_col_types=F)
bp=as.data.frame(bp)
colnames(bp)[1]="Type"
colnames(bp)[4]="Number"
colnames(bp)[7]="FDR"
bp$Type=factor(bp$Type,levels=rev(unique(bp$Type)))

cols <- c("#F9EBEA","#CD6155")
getPalette= colorRampPalette(cols)

p <- ggplot(bp, aes(y=Number, x=Type,fill=-log10(FDR))) + scale_fill_gradientn(name= "-log10(FDR)", colours=getPalette(100), na.value="grey")


p <- p + geom_bar(position="stack", stat="identity",width=0.5,)+

p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, 8))

p <- p + theme(axis.ticks = element_blank())+coord_flip() + theme_bw() + theme_classic()+ggtitle("GO-terms enrichment")+ylab("Overlap")

#options(repr.plot.width=9, repr.plot.height=3)
print(p) 

pdf(paste("Desktop/Projects/CPTAC3_PDAC/Analysis/2020-09-11/GO_terms_2020-09-11.pdf",sep=""), width=9, height=3,useDingbats=FALSE)
print(p)
dev.off()    

	