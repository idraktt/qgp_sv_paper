library(tidyr)
library(ggplot2)


svtypes.file <- "data/SV_colors.txt"
pops.file <- "data/QGP_population_colors.txt"
QGP_annotated <- read.delim("data/QGP_hq.tsv.processed.tsv")
table(QGP_annotated$SV.type)
QGP_annotated<-QGP_annotated[!duplicated(QGP_annotated$AnnotSV.ID), ]
table(QGP_annotated$SV.type)
qgp_subpop_freq <-read.delim("data/QGP_hq_allele_freq_subpop.tsv")
qgp_subpop_freq <- qgp_subpop_freq[order(qgp_subpop_freq$ID, -abs(qgp_subpop_freq$AF) ), ] 
nrow(qgp_subpop_freq)
qgp_subpop_freq<-qgp_subpop_freq[!duplicated(qgp_subpop_freq$ID), ]
nrow(qgp_subpop_freq)


qgp_white_list <- read.table("data/qgp_bam_white.list", quote="\"", comment.char="")
autosomes <- read.delim("data/autosomes.bed")
qgp_mapping <- read.delim("data/mapping.tsv")
qgp_svtk_annotation <- read.delim("data/QGP_hq.svtk.bed")
cols.to.drop <- c("X.chrom","start","end","svtype","END","MAPQ",
                  "RE","IMPRECISE","PRECISE","AVGLEN","SVMETHOD",
                  "SVTYPE","SUPP","STRANDS","AF","NSAMP","QGP_ADM_AF",
                  "QGP_ADM_AC","GP_ADM_AF","QGP_ADM_AN",
                  "QGP_ADM_HOMO","QGP_AFR_AC","QGP_AFR_AF",
                  "QGP_AFR_AN","QGP_AFR_HOMO","QGP_GAR_AC",
                  "QGP_GAR_AF","QGP_GAR_AN","QGP_GAR_HOMO",
                  "QGP_PAR_AC","QGP_PAR_AF","QGP_PAR_AN",
                  "QGP_PAR_HOMO","QGP_SAS_AC","QGP_SAS_AF",
                  "QGP_SAS_AN","QGP_SAS_HOMO","QGP_WEP_AC",
                  "QGP_WEP_AF","QGP_WEP_AN","QGP_WEP_HOMO")
qgp_svtk_annotation <- qgp_svtk_annotation[, -which(colnames(qgp_svtk_annotation) %in% cols.to.drop)]


###Sets populations & colors
if(!is.null(pops.file)){
  pops <- read.table(pops.file,sep="\t",header=T,comment.char="",check.names=F)
}

Population_structure_labels <- read.delim("data/Population_structure_labels", header=FALSE)
qgp_population_structure_clusters <- read.delim("data/qgp_population_structure_clusters.txt", header=FALSE)
population<-merge(qgp_population_structure_clusters, Population_structure_labels, by.x="V2", by.y="V1")
population<-merge(population, pops, by.x="V2.y", by.y="name")

###Sets sv types & colors
if(!is.null(svtypes.file)){
  svtypes <- read.table(svtypes.file,sep="\t",header=F,comment.char="",check.names=F)
  svtypes <- as.data.frame(apply(svtypes,2,as.character))
  colnames(svtypes) <- c("svtype","color")
}else{
  require(RColorBrewer,quietly=T)
  svtypes.v <- unique(dat$SVTYPE)
  svtypes.c <- brewer.pal(length(svtypes.v),"Dark2")
  svtypes <- data.frame("svtype"=svtypes.v,
                        "color"=svtypes.c)
}



dat<-merge(QGP_annotated, qgp_subpop_freq, by.x="ID", by.y="ID")
white_list_samples<-merge(qgp_mapping,qgp_white_list,by.x = "BAM_LOC",by.y = "V1")
qgp_pop_mapping_final<-merge(white_list_samples,population,by.x = "NAME", by.y = "V1")

write.table(qgp_pop_mapping_final, file='data/QGP_pop_mapping_final.tsv', quote=FALSE, sep='\t')



dat$AN<-dat$NSAMP
dat$AC<-dat$QGP_ADM_AC + dat$QGP_AFR_AC + dat$QGP_WEP_AC + dat$QGP_GAR_AC + dat$QGP_PAR_AC + dat$QGP_SAS_AC
dat$SVLEN<-dat$SV.end - dat$SV.start
dat$SVTYPE<-dat$SV.type
dat<-dat[which(dat$AN>0),]

dat$N_BI_GENOS<-6141
dat$QGP_ADM_N_BI_GENOS<-1186
dat$QGP_AFR_N_BI_GENOS<-179
dat$QGP_GAR_N_BI_GENOS<-2312
dat$QGP_PAR_N_BI_GENOS<-1052
dat$QGP_SAS_N_BI_GENOS<-38
dat$QGP_WEP_N_BI_GENOS<-1374

dat$N_HOMALT<-dat$QGP_ADM_HOMO + dat$QGP_AFR_HOMO + dat$QGP_GAR_HOMO + dat$QGP_PAR_HOMO + dat$QGP_SAS_HOMO + dat$QGP_WEP_HOMO
dat$QGP_ADM_N_HOMALT<-dat$QGP_ADM_HOMO
dat$QGP_AFR_N_HOMALT<-dat$QGP_AFR_HOMO
dat$QGP_GAR_N_HOMALT<-dat$QGP_GAR_HOMO
dat$QGP_PAR_N_HOMALT<-dat$QGP_PAR_HOMO
dat$QGP_SAS_N_HOMALT<-dat$QGP_SAS_HOMO
dat$QGP_WEP_N_HOMALT<-dat$QGP_WEP_HOMO
dat$QGP_ADM_N_HET<-(dat$QGP_ADM_AC - (dat$QGP_ADM_N_HOMALT * 2))
dat$QGP_AFR_N_HET<-(dat$QGP_AFR_AC - (dat$QGP_AFR_N_HOMALT * 2))
dat$QGP_GAR_N_HET<-(dat$QGP_GAR_AC - (dat$QGP_GAR_N_HOMALT * 2))
dat$QGP_PAR_N_HET<-(dat$QGP_PAR_AC - (dat$QGP_PAR_N_HOMALT * 2))
dat$QGP_SAS_N_HET<-(dat$QGP_SAS_AC - (dat$QGP_SAS_N_HOMALT * 2))
dat$QGP_WEP_N_HET<-(dat$QGP_WEP_AC - (dat$QGP_WEP_N_HOMALT * 2))
dat$QGP_ADM_N_HOMREF<-dat$QGP_ADM_N_BI_GENOS - dat$QGP_ADM_N_HOMALT - dat$QGP_ADM_N_HET
dat$QGP_AFR_N_HOMREF<-dat$QGP_AFR_N_BI_GENOS - dat$QGP_AFR_N_HOMALT - dat$QGP_AFR_N_HET
dat$QGP_GAR_N_HOMREF<-dat$QGP_GAR_N_BI_GENOS - dat$QGP_GAR_N_HOMALT - dat$QGP_GAR_N_HET
dat$QGP_PAR_N_HOMREF<-dat$QGP_PAR_N_BI_GENOS - dat$QGP_PAR_N_HOMALT - dat$QGP_PAR_N_HET
dat$QGP_SAS_N_HOMREF<-dat$QGP_SAS_N_BI_GENOS - dat$QGP_SAS_N_HOMALT - dat$QGP_SAS_N_HET
dat$QGP_WEP_N_HOMREF<-dat$QGP_WEP_N_BI_GENOS - dat$QGP_WEP_N_HOMALT - dat$QGP_WEP_N_HET

dat$N_HET<-(dat$AC - (dat$N_HOMALT * 2))
dat$N_HOM_REF<-dat$N_BI_GENOS - dat$N_HOMALT - dat$N_HET

write.table(dat, file='data/QGP.tsv', quote=FALSE, sep='\t')

n.genos.idx <- which(colnames(dat)==paste(pops$pop, "_N_BI_GENOS", sep=""))
n.homref.idx <- which(colnames(dat)==paste(pops$pop, "_N_HOMREF", sep=""))
n.het.idx <- which(colnames(dat)==paste(pops$poppop, "_N_HET", sep=""))
n.homalt.idx <- which(colnames(dat)==paste(pops$poppop, "_N_HOMALT", sep=""))


autosomes.in <- "data/autosomes.bed"
centromeres.in <- "data/centromeres.bed"
nmask.in <- "data/nmasks.bed"

sites.counts.all <- read.delim("data/QGP_hq.total_counts.txt")
sites.counts.all.wide <- spread(sites.counts.all, svtype, count)
sites.counts.homozygous <- read.delim("data/QGP_hq_total_counts.homo.txt")
sites.counts.homozygous.wide <- spread(sites.counts.homozygous, svtype, count)
sites.counts.rare <-read.delim("data/QGP_hq_total_counts.rare.txt")
sites.counts.rare.wide <- spread(sites.counts.rare, svtype, count)
sites.counts.singleton <-read.delim("data/QGP_hq_total_counts.singleton.txt")
sites.counts.singleton.wide <- spread(sites.counts.singleton, svtype, count)

sites.counts.all.wide.pop<-merge(sites.counts.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "SAMPLE_NAME")
sites.counts.homozygous.wide.pop<-merge(sites.counts.homozygous.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "SAMPLE_NAME")
sites.counts.rare.wide.pop<-merge(sites.counts.rare.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "SAMPLE_NAME")
sites.counts.singleton.wide.pop<-merge(sites.counts.singleton.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "SAMPLE_NAME")

all.sites<-setNames(data.frame(matrix(ncol = 9, nrow = 6140)), c("sample","pop",
                                                                 "all.DEL","all.MCNV_DEL",
                                                                 "all.DUP","all.MCNV_DUP",
                                                                 "all.INS","all.INV","all.CPX"))

all.sites$sample<-sites.counts.all.wide.pop$sample
all.sites$pop<-sites.counts.all.wide.pop$pop
all.sites$all.DEL<-sites.counts.all.wide.pop$DEL
all.sites$all.MCNV_DEL<-0
all.sites$all.DUP<-sites.counts.all.wide.pop$DUP
all.sites$all.MCNV_DUP<-0
all.sites$all.INS<-0
all.sites$all.INV<-sites.counts.all.wide.pop$INV
all.sites$all.CPX<-0

QGP_hq.persample_counts.all <- read.delim("data/QGP_hq.persample_counts.all.tsv", header=TRUE)

sites.counts.homozygous.wide.pop<-merge(QGP_hq.persample_counts.all,qgp_pop_mapping_final,by.x = "SAMPLE",by.y = "SAMPLE_NAME")

homo.sites<-setNames(data.frame(matrix(ncol = 9, nrow = 6140)), c("sample","pop",
                                                                 "all.DEL","all.MCNV_DEL",
                                                                 "all.DUP","all.MCNV_DUP",
                                                                 "all.INS","all.INV","all.CPX"))

homo.sites$sample<-sites.counts.homozygous.wide.pop$SAMPLE
homo.sites$pop<-sites.counts.homozygous.wide.pop$pop
homo.sites$all.DEL<-sites.counts.homozygous.wide.pop$HOM_DEL
homo.sites$all.MCNV_DEL<-0
homo.sites$all.DUP<-sites.counts.homozygous.wide.pop$HOM_DUP
homo.sites$all.MCNV_DUP<-0
homo.sites$all.INS<-0
homo.sites$all.INV<-sites.counts.homozygous.wide.pop$HOM_INV
homo.sites$all.CPX<-0

rare.sites<-setNames(data.frame(matrix(ncol = 9, nrow = 6140)), c("sample","pop",
                                                                  "all.DEL","all.MCNV_DEL",
                                                                  "all.DUP","all.MCNV_DUP",
                                                                  "all.INS","all.INV","all.CPX"))

rare.sites$sample<-sites.counts.rare.wide.pop$sample
rare.sites$pop<-sites.counts.rare.wide.pop$pop
rare.sites$all.DEL<-sites.counts.rare.wide.pop$DEL
rare.sites$all.MCNV_DEL<-0
rare.sites$all.DUP<-sites.counts.rare.wide.pop$DUP
rare.sites$all.MCNV_DUP<-0
rare.sites$all.INS<-0
rare.sites$all.INV<-sites.counts.rare.wide.pop$INV
rare.sites$all.CPX<-0

singleton.sites<-setNames(data.frame(matrix(ncol = 9, nrow = 5081)), c("sample","pop",
                                                                  "all.DEL","all.MCNV_DEL",
                                                                  "all.DUP","all.MCNV_DUP",
                                                                  "all.INS","all.INV","all.CPX"))

singleton.sites$sample<-sites.counts.singleton.wide.pop$sample
singleton.sites$pop<-sites.counts.singleton.wide.pop$pop
singleton.sites$all.DEL<-sites.counts.singleton.wide.pop$DEL
singleton.sites$all.MCNV_DEL<-0
singleton.sites$all.DUP<-sites.counts.singleton.wide.pop$DUP
singleton.sites$all.MCNV_DUP<-0
singleton.sites$all.INS<-0
singleton.sites$all.INV<-sites.counts.singleton.wide.pop$INV
singleton.sites$all.CPX<-0


all.sites.wdat <- getWaterfallData(dat=all.sites)
homo.sites.wdat <- getWaterfallData(dat=homo.sites)
rare.sites.wdat <- getWaterfallData(dat=rare.sites)
singleton.sites.wdat <- getWaterfallData(dat=singleton.sites)


output_folder<-"/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/reports/sample_final_annotated/"

genic.counts.wide.all.wide<-read.delim(paste(output_folder,"genic_counts_wide.tsv"))
exonic_counts_wide.all.wide<-read.delim(paste(output_folder,"exonic_counts_wide.tsv"))
genic_counts_novel_wide.all.wide<-read.delim(paste(output_folder,"genic_counts_novel_wide.tsv"))
exonic_counts_novel_wide.all.wide<-read.delim(paste(output_folder,"exonic_counts_novel_wide.tsv"))
total_genes_wide.all.wide<-read.delim(paste(output_folder,"total_genes_wide.tsv"))
total_exonic_genes_wide.all.wide<-read.delim(paste(output_folder,"total_exonic_genes_counts.tsv"))
omim_genes_wide.all.wide<-read.delim(paste(output_folder,"omim_genes_counts.tsv"))
acmg_counts_wide.all.wide<-read.delim(paste(output_folder,"acmg_counts.tsv"))

genic_counts_wide.all.wide.pop<-merge(genic.counts.wide.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "NAME")
exonic_counts_wide.all.wide.pop<-merge(exonic_counts_wide.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "NAME")
genic_counts_novel_wide.all.pop<-merge(genic_counts_novel_wide.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "NAME")
exonic_counts_novel_wide.all.wide.pop<-merge(exonic_counts_novel_wide.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "NAME")
total_genes_wide.all.wide.pop<-merge(total_genes_wide.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "NAME")
total_exonic_genes_wide.all.wide.pop<-merge(total_exonic_genes_wide.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "NAME")
omim_genes_wide.all.wide.pop<-merge(omim_genes_wide.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "NAME")
acmg_counts_wide.all.wide.pop<-merge(acmg_counts_wide.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "NAME")

genic.sites<-setNames(data.frame(matrix(ncol = 9, nrow = 6141)), c("sample","pop",
                                                                  "all.DEL","all.MCNV_DEL",
                                                                  "all.DUP","all.MCNV_DUP",
                                                                  "all.INS","all.INV","all.CPX"))

genic.sites$sample<-genic_counts_wide.all.wide.pop$SAMPLE
genic.sites$pop<-genic_counts_wide.all.wide.pop$pop
genic.sites$all.DEL<-genic_counts_wide.all.wide.pop$DEL
genic.sites$all.MCNV_DEL<-0
genic.sites$all.DUP<-genic_counts_wide.all.wide.pop$DUP
genic.sites$all.MCNV_DUP<-0
genic.sites$all.INS<-0
genic.sites$all.INV<-genic_counts_wide.all.wide.pop$INV
genic.sites$all.CPX<-0

exonic.sites<-setNames(data.frame(matrix(ncol = 9, nrow = 6141)), c("sample","pop",
                                                                   "all.DEL","all.MCNV_DEL",
                                                                   "all.DUP","all.MCNV_DUP",
                                                                   "all.INS","all.INV","all.CPX"))

exonic.sites$sample<-exonic_counts_wide.all.wide.pop$SAMPLE
exonic.sites$pop<-exonic_counts_wide.all.wide.pop$pop
exonic.sites$all.DEL<-exonic_counts_wide.all.wide.pop$DEL
exonic.sites$all.MCNV_DEL<-0
exonic.sites$all.DUP<-exonic_counts_wide.all.wide.pop$DUP
exonic.sites$all.MCNV_DUP<-0
exonic.sites$all.INS<-0
exonic.sites$all.INV<-exonic_counts_wide.all.wide.pop$INV
exonic.sites$all.CPX<-0

genic.sites.novel<-setNames(data.frame(matrix(ncol = 9, nrow = 6141)), c("sample","pop",
                                                                    "all.DEL","all.MCNV_DEL",
                                                                    "all.DUP","all.MCNV_DUP",
                                                                    "all.INS","all.INV","all.CPX"))

genic.sites.novel$sample<-genic_counts_novel_wide.all.pop$SAMPLE
genic.sites.novel$pop<-genic_counts_novel_wide.all.pop$pop
genic.sites.novel$all.DEL<-genic_counts_novel_wide.all.pop$DEL
genic.sites.novel$all.MCNV_DEL<-0
genic.sites.novel$all.DUP<-genic_counts_novel_wide.all.pop$DUP
genic.sites.novel$all.MCNV_DUP<-0
genic.sites.novel$all.INS<-0
genic.sites.novel$all.INV<-genic_counts_novel_wide.all.pop$INV
genic.sites.novel$all.CPX<-0

exonic.sites.novel<-setNames(data.frame(matrix(ncol = 9, nrow = 6141)), c("sample","pop",
                                                                   "all.DEL","all.MCNV_DEL",
                                                                   "all.DUP","all.MCNV_DUP",
                                                                   "all.INS","all.INV","all.CPX"))

exonic.sites.novel$sample<-exonic_counts_novel_wide.all.wide.pop$SAMPLE
exonic.sites.novel$pop<-exonic_counts_novel_wide.all.wide.pop$pop
exonic.sites.novel$all.DEL<-exonic_counts_novel_wide.all.wide.pop$DEL
exonic.sites.novel$all.MCNV_DEL<-0
exonic.sites.novel$all.DUP<-exonic_counts_novel_wide.all.wide.pop$DUP
exonic.sites.novel$all.MCNV_DUP<-0
exonic.sites.novel$all.INS<-0
exonic.sites.novel$all.INV<-exonic_counts_novel_wide.all.wide.pop$INV
exonic.sites.novel$all.CPX<-0

total.genes<-setNames(data.frame(matrix(ncol = 9, nrow = 6141)), c("sample","pop",
                                                                          "all.DEL","all.MCNV_DEL",
                                                                          "all.DUP","all.MCNV_DUP",
                                                                          "all.INS","all.INV","all.CPX"))

total.genes$sample<-total_genes_wide.all.wide.pop$SAMPLE
total.genes$pop<-total_genes_wide.all.wide.pop$pop
total.genes$all.DEL<-total_genes_wide.all.wide.pop$DEL
total.genes$all.MCNV_DEL<-0
total.genes$all.DUP<-total_genes_wide.all.wide.pop$DUP
total.genes$all.MCNV_DUP<-0
total.genes$all.INS<-0
total.genes$all.INV<-total_genes_wide.all.wide.pop$INV
total.genes$all.CPX<-0

total.exonic.genes<-setNames(data.frame(matrix(ncol = 9, nrow = 6141)), c("sample","pop",
                                                                   "all.DEL","all.MCNV_DEL",
                                                                   "all.DUP","all.MCNV_DUP",
                                                                   "all.INS","all.INV","all.CPX"))

total.exonic.genes$sample<-total_exonic_genes_wide.all.wide.pop$SAMPLE
total.exonic.genes$pop<-total_exonic_genes_wide.all.wide.pop$pop
total.exonic.genes$all.DEL<-total_exonic_genes_wide.all.wide.pop$DEL
total.exonic.genes$all.MCNV_DEL<-0
total.exonic.genes$all.DUP<-total_exonic_genes_wide.all.wide.pop$DUP
total.exonic.genes$all.MCNV_DUP<-0
total.exonic.genes$all.INS<-0
total.exonic.genes$all.INV<-total_exonic_genes_wide.all.wide.pop$INV
total.exonic.genes$all.CPX<-0

omim.genes<-setNames(data.frame(matrix(ncol = 9, nrow = 6141)), c("sample","pop",
                                                                          "all.DEL","all.MCNV_DEL",
                                                                          "all.DUP","all.MCNV_DUP",
                                                                          "all.INS","all.INV","all.CPX"))

omim.genes$sample<-omim_genes_wide.all.wide.pop$SAMPLE
omim.genes$pop<-omim_genes_wide.all.wide.pop$pop
omim.genes$all.DEL<-omim_genes_wide.all.wide.pop$DEL
omim.genes$all.MCNV_DEL<-0
omim.genes$all.DUP<-omim_genes_wide.all.wide.pop$DUP
omim.genes$all.MCNV_DUP<-0
omim.genes$all.INS<-0
omim.genes$all.INV<-omim_genes_wide.all.wide.pop$INV
omim.genes$all.CPX<-0

acmg.genes<-setNames(data.frame(matrix(ncol = 9, nrow = 6141)), c("sample","pop",
                                                                  "all.DEL","all.MCNV_DEL",
                                                                  "all.DUP","all.MCNV_DUP",
                                                                  "all.INS","all.INV","all.CPX"))

acmg.genes$sample<-acmg_counts_wide.all.wide.pop$SAMPLE
acmg.genes$pop<-acmg_counts_wide.all.wide.pop$pop
acmg.genes$all.DEL<-acmg_counts_wide.all.wide.pop$DEL
acmg.genes$all.MCNV_DEL<-0
acmg.genes$all.DUP<-acmg_counts_wide.all.wide.pop$DUP
acmg.genes$all.MCNV_DUP<-0
acmg.genes$all.INS<-0
acmg.genes$all.INV<-acmg_counts_wide.all.wide.pop$INV
acmg.genes$all.CPX<-0

genic.sites.wdat <- getWaterfallData(dat=genic.sites)
exonic.sites.wdat <- getWaterfallData(dat=exonic.sites)
genic.sites.novel.wdat <- getWaterfallData(dat=genic.sites.novel)
exonic.sites.novel.wdat <- getWaterfallData(dat=exonic.sites.novel)
total.genes.wdat <- getWaterfallData(dat=total.genes)
total.exonic.genes.wdat <- getWaterfallData(dat=total.exonic.genes)
omim.genes.wdat <- getWaterfallData(dat=omim.genes)
acmg.genes.wdat <- getWaterfallData(dat=acmg.genes)


###Process input data
cat("NOW LOADING DATA\n")
autosomes.raw <- loadPileup(autosomes.in,centromeres.in,smooth=F)
autosomes <- loadPileup(autosomes.in,centromeres.in)

###Set analysis parameters & read helper files
meta.svtypes <- c("ALL","DEL","DUP","INV","INS","TRA")
Nmasks <- read.table(nmask.in,header=F)
names(Nmasks) <- c("chr","start","end")

###Gather meta-chromosome averages per percentile
meta.means <- lapply(meta.svtypes,function(svtype){
  m <- metaAverage(dat=autosomes,SVTYPE=svtype,n.bins=500)
})
names(meta.means) <- meta.svtypes

dat <- merge(dat,qgp_svtk_annotation,by.x = "ID",by.y = "name")

# Basic barplot
df <- data.frame(Study=c("GNOMAD-SV", "QGP", "1000G","GONL","GTEX"),
                 Counts=c(14237, 6141, 2504, 750,147))
p<-ggplot(data=df, aes(x=Study, y=Counts, fill=Study)) +
  geom_text(aes(label=Counts), vjust=0.6, color="white", size=6.5)+
  theme_minimal() + 
  geom_bar(stat="identity") + 
  scale_x_discrete(limits=c("GTEX", "GONL", "1000G","QGP","GNOMAD-SV")) +
  labs(title="Counts by study", 
       x="Study", y = "Sample count") + 
  geom_text(aes(label=Counts),hjust = -0.5, vjust=0, color="black", size=4.5)
# Horizontal bar plot
p + coord_flip()
###Plot simple bars of counts -- done
cat("NOW STARTING COUNT OF SV BY TYPE\n")
pdf(paste(OUTDIR, "/", prefix, ".site_counts_by_type.pdf", sep=""), 
    height=2.75, width=2.25)
plot.totalCounts(dat=dat, svtypes=svtypes, thousandG=F)
dev.off()

names(dat)[names(dat) == "SVLEN.x"] <- "SVLEN"


###Plot size distribution
cat("NOW STARTING SIZE DISTRIBUTIONS\n")
pdf(paste(OUTDIR, "/", prefix, ".size_distribution_by_type.pdf", sep=""), 
    height=2, width=2.9)
plot.sizes(dat=dat, svtypes=svtypes)
dev.off()

###Plot frequency distribution by SVTYPE
cat("NOW STARTING FREQUENCY DISTRIBUTIONS\n")
pdf(paste(OUTDIR, "/", prefix, ".site_frequency_distributions_bySVtype.pdf", sep=""), 
    height=2, width=2.3)
plot.freqByType(dat=dat)
dev.off()



cat("NOW STARTING HARDY-WEINBERG EQUILIBRIUM CALCULATIONS\n")
png(paste(OUTDIR, "/", prefix, ".HWE_all_samples.png", sep=""), 
    height=1000, width=1000, res=400)
plot.HWE(dat=dat[which(dat$SV.type == "DEL"),], pop=NULL, title="All Samples")
dev.off()
if(!is.null(pops.file)){
  sapply(1:nrow(pops), function(i){
    pop <- pops$pop[i]
    if(pop=="AFR"){
      pop.name <- "African"
    }else{
      pop.name <- pops$name[i]
    }
    print(pop)
    png(paste(OUTDIR, "/", prefix, ".HardyWeinberg_ternary_plot.", pop, "_samples.png", sep=""), 
        height=1000, width=1000, res=400)
    plot.HWE(dat=dat[which(dat$SV.type == "DEL"),], pop=pop, title=pop.name)
    dev.off()
  })
}


###Plot cross-population allele frequency correlations -- Done
cat("NOW STARTING CROSS-POPULATION ALLELE FREQUENCY CORRELATIONS\n")
sapply(1:nrow(pops), function(a){
  sapply(2:nrow(pops), function(b){
    if(b>a){
      popA <- pops$pop[a]
      popB <- pops$pop[b]
      if(popA != "OTH" & popB != "OTH"){
        png(paste(OUTDIR, "/", prefix, ".frequency_correlation.", popA, "_vs_", popB, ".png", sep=""), 
            height=1000, width=1000, res=400)
        plot.crossPopCorr(dat, pops, popA=popA, popB=popB)
        dev.off()
      }
    }
  })
})

###Generate genome-wide panel of idiograms - ALL SV
png(paste(OUTDIR,"/",prefix,".SV_density_idiograms.png",sep=""),
    height=800,width=2600,res=300,bg="white")
generateCovPlotsPerChrom(mat=lapply(1:22,function(k){return(as.data.frame(autosomes[which(autosomes$chr==k),c(1:3,which(colnames(autosomes)=="ALL"))]))}),
                         colors="gray10",
                         contigs.top=1:5,
                         contigs.middle=6:12,
                         contigs.bottom=13:22,
                         fill=T)
dev.off()
png(paste(OUTDIR,"/",prefix,".SV_density_idiograms.shorter.png",sep=""),
    height=600,width=(2/3)*2600,res=300,bg="white")
generateCovPlotsPerChrom(mat=lapply(1:22,function(k){return(as.data.frame(autosomes[which(autosomes$chr==k),c(1:3,which(colnames(autosomes)=="ALL"))]))}),
                         colors="gray10",
                         contigs.top=1:5,
                         contigs.middle=6:12,
                         contigs.bottom=13:22,
                         fill=T)
dev.off()

###Generate "meta-chromosome" plots per SVTYPE
#Large plot: all svtypes
png(paste(OUTDIR,"/",prefix,".SV_density_metaIdiograms.All_SV.png",sep=""),
    height=750,width=2.5*(2600/7),res=300,bg="white")
par(mar=c(1.5,2,0.2,1))
plot.metaDist(meta.dat=meta.means[[which(names(meta.means)=="ALL")]],
              color="gray10",norm=T,
              xlabel="Meta-chromosome (Mean of Autosomes)")
dev.off() 

#Smaller panels: individual svtypes
png(paste(OUTDIR, "/", prefix, ".SV_density_metaIdiograms.bySVTYPE.png", sep=""), 
    height=650, width=5*(2600/7), res=300, bg="white")
par(mfrow=c(2, length(meta.svtypes[-1])+1/2), 
    mar=c(1.5, 2, 1.5, 1))
lapply(meta.svtypes[-1], function(svtype){
  svt.col <- svtypes$color[which(svtypes$svtype==svtype)]
  plot.metaDist(meta.dat=meta.means[[which(names(meta.means)==svtype)]], 
                color=svt.col, norm=T)
})
dev.off()

###Generate genome-wide panel of idiograms - ALL SV
png(paste(OUTDIR,"/",prefix,".SV_density_idiograms_bySVTYPE.png",sep=""),
    height=800,width=2600,res=300,bg="white")
generateCovPlotsPerChrom(mat=lapply(1:22,function(k){return(as.data.frame(autosomes[which(autosomes$chr==k),c(1:3,which(colnames(autosomes) %in% c("DEL","DUP","INV")))]))}),
                         colors=c(svtypes$color[which(svtypes$svtype=="DEL")],
                                  svtypes$color[which(svtypes$svtype=="DUP")],
                                  svtypes$color[which(svtypes$svtype=="INV")],
                                  svtypes$color[which(svtypes$svtype=="INS")],
                                  svtypes$color[which(svtypes$svtype=="TRA")]),
                         contigs.top=1:5,
                         contigs.middle=6:12,
                         contigs.bottom=13:22,
                         fill=F,norm=T)
dev.off()

###Plot meta means by ter/int/cen
pdf(paste(OUTDIR,"/",prefix,".SV_density_meta.byChromContext.dotplots.pdf",sep=""),
    height=1.9,width=3.25)
plot.metaByContext(autosomes,meta.svtypes,
                   colors=c("gray10",
                            svtypes$color[which(svtypes$svtype=="DEL")],
                            svtypes$color[which(svtypes$svtype=="DUP")],
                            svtypes$color[which(svtypes$svtype=="INV")],
                            svtypes$color[which(svtypes$svtype=="INS")],
                            svtypes$color[which(svtypes$svtype=="TRA")]))
dev.off()

#Waterfall
#Generate waterfall plots
###Waterfall plots of variant counts per sample
site.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")],
                 svtypes$color[which(svtypes$svtype=="DEL")],
                 svtypes$color[which(svtypes$svtype=="DUP")],
                 svtypes$color[which(svtypes$svtype=="DUP")],
                 svtypes$color[which(svtypes$svtype=="INS")],
                 svtypes$color[which(svtypes$svtype=="INV")],
                 svtypes$color[which(svtypes$svtype=="CPX")])
site.labels <- c("DEL","MCNV (Loss)","DUP","MCNV (Gain)",
                 "INS","INV","CPX")

png(paste(OUTDIR,"/",prefix,".sites_perSample.waterfall.png",sep=""),
    width=9*300,height=8*300,res=400)
layout(matrix(1:6,nrow=6),heights=c(3,2,1.5,1.5))
par(mar=c(0.5,3.5,1.75,5.5),bty="n")
plotWaterfall(wdat=all.sites.wdat,colors=site.colors,
              ylabel="All",titles=T,cat.labels=site.labels)
par(mar=c(0.5,3.5,0.5,5.5),bty="n")
plotWaterfall(wdat=homo.sites.wdat,colors=site.colors,
              ylabel="Homozygous",y.at.nmax=4,cat.labels=site.labels)
plotWaterfall(wdat=rare.sites.wdat,colors=site.colors,
              ylabel="Rare",y.at.nmax=3,cat.labels=site.labels)
plotWaterfall(wdat=singleton.sites.wdat,colors=site.colors,
              ylabel="Singleton",y.at.nmax=3,cat.labels=site.labels)
dev.off()

#Waterfall
#Generate waterfall Gene plots
###Waterfall plots of variant counts per sample
site.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")],
                 svtypes$color[which(svtypes$svtype=="DEL")],
                 svtypes$color[which(svtypes$svtype=="DUP")],
                 svtypes$color[which(svtypes$svtype=="DUP")],
                 svtypes$color[which(svtypes$svtype=="INS")],
                 svtypes$color[which(svtypes$svtype=="INV")],
                 svtypes$color[which(svtypes$svtype=="CPX")])
site.labels <- c("DEL","MCNV (Loss)","DUP","MCNV (Gain)",
                 "INS","INV","CPX")

png(paste(OUTDIR,"/",prefix,".sites_perSample.waterfall.png",sep=""),
    width=9*300,height=8*300,res=400)
layout(matrix(1:8,nrow=8),heights=c(3,2,2,2,2,2,2,2))
par(mar=c(0.5,3.5,1.75,5.5),bty="n")
plotWaterfall(wdat=total.genes.wdat,colors=site.colors,
              ylabel="Genes",titles=T,y.at.nmax=4,cat.labels=site.labels)
par(mar=c(0.5,3.5,0.5,5.5),bty="n")
plotWaterfall(wdat=genic.sites.wdat,colors=site.colors,
              ylabel="Genic SVs",y.at.nmax=4,cat.labels=site.labels)
plotWaterfall(wdat=genic.sites.novel.wdat,colors=site.colors,
              ylabel="Novel Genic",y.at.nmax=4,cat.labels=site.labels)
plotWaterfall(wdat=exonic.sites.wdat,colors=site.colors,
              ylabel="Exonic SVs",y.at.nmax=4,cat.labels=site.labels)
plotWaterfall(wdat=total.exonic.genes.wdat,colors=site.colors,
              ylabel="Exonic Genes",y.at.nmax=4,cat.labels=site.labels)
plotWaterfall(wdat=exonic.sites.novel.wdat,colors=site.colors,
              ylabel="Novel Exonic",y.at.nmax=4,cat.labels=site.labels)
plotWaterfall(wdat=omim.genes.wdat,colors=site.colors,
              ylabel="OMIM genes",y.at.nmax=4,cat.labels=site.labels)
plotWaterfall(wdat=acmg.genes.wdat,colors=site.colors,
              ylabel="ACMG genes",y.at.nmax=4,cat.labels=site.labels)
dev.off()


#Generate smaller plot of just counts per sample for main figure
png(paste(OUTDIR,"/",prefix,".sites_perSample.waterfall.small.png",sep=""),
    width=10.5*300,height=1.5*300,res=400)
par(mar=c(0.5,3.5,1.1,5.5),bty="n")
plotWaterfall(wdat=all.sites.wdat,colors=site.colors,
              ylabel="SVs per Genome",titles=T,cat.labels=site.labels,y.at.nmax=10,cex.ylabs=0.6)
dev.off()

#Generate sina plot of counts per class by ancestry
png(paste(OUTDIR,"/",prefix,".sites_perSample.sina.small.png",sep=""),
    width=8*300,height=1.1*300,res=400)
wrapper.gridSinaSvtypeByPop(wdat=all.sites.wdat)
dev.off()

#Generate sina plot of counts per class by ancestry
png(paste(OUTDIR,"/",prefix,".sites_perSample.sina.small.homo.png",sep=""),
    width=8*300,height=1.1*300,res=400)
wrapper.gridSinaSvtypeByPop(wdat=homo.sites.wdat)
dev.off()

##PCA PLOT
pca0.5 <-read.table("/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/pca/0.5.eigenvec")
colnames(pca0.5)[1:12]=c("INDEX","ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
pca0.5 <- merge(pca0.5,qgp_pop_mapping_final,by.x = "ID", by.y = "SAMPLE_NAME")

ggplot(pca0.5,aes(x=PC1,y=PC2,color=pop))+geom_point()+ scale_color_manual(values=c("#BAB0AC", "#F28E2B", "#4E79A7","#76B7B2", "#EDC948", "#E15759"),
                                                                           name="Populations",
                                                                           labels=c("AdmixedPop", "East Africans", "Common Arabs","Peninsular Arabs","South Asians","Persians")) + theme_minimal()

ggplot(pca0.5,aes(x=PC1,y=PC3,color=pop))+geom_point()+ scale_color_manual(values=c("#BAB0AC", "#F28E2B", "#4E79A7","#76B7B2", "#EDC948", "#E15759"),
                                                                           name="Populations",
                                                                           labels=c("AdmixedPop", "East Africans", "Common Arabs","Peninsular Arabs","South Asians","Persians")) + theme_minimal()

ggplot(pca0.5,aes(x=PC2,y=PC3,color=pop))+geom_point()+ scale_color_manual(values=c("#BAB0AC", "#F28E2B", "#4E79A7","#76B7B2", "#EDC948", "#E15759"),
                                                                           name="Populations",
                                                                           labels=c("AdmixedPop", "East Africans", "Common Arabs","Peninsular Arabs","South Asians","Persians")) + theme_minimal()

#Plot three-panel of SV with and without predicted coding effects
pdf(paste(OUTDIR, "/", prefix, ".site_frequency_distributions_bySVtype.by_genomic_context.three_panels.pdf", sep=""), 
    height=1.75, width=5.5)
par(mfrow=c(1, 3))

plot.freqByType(dat=dat[which(!is.na(dat$LOF)
                              | !is.na(dat$COPY_GAIN)
                              | !is.na(dat$DUP_LOF)), ], 
                ymin=0.20, axlabel.cex=0.8)
plot.freqByType(dat=dat[which(is.na(dat$LOF)
                              & is.na(dat$COPY_GAIN)
                              & is.na(dat$DUP_LOF)
                              & (!is.na(dat$INTRONIC))
                              | !is.na(dat$DUP_PARTIAL)), ], 
                ymin=0.20, axlabel.cex=0.8)
plot.freqByType(dat=dat[which(is.na(dat$INTRONIC)), ], 
                ymin=0.20, axlabel.cex=0.8)
dev.off()

#Plot frequency distribution per SVTYPE by coding context
pdf(paste(OUTDIR, "/", prefix, ".site_frequency_distributions_bySVtype.by_genomic_context.five_panels.pdf", sep=""), 
    height=2*1.75, width=5.5)
par(mfrow=c(2, 3))
sapply(c("DEL", "DUP", "INV"), function(svtype){
  plot.freqByContext(dat=dat, svtype=svtype, ymin=0.25, axlabel.cex=0.8, parmar=c(2.5, 4.25, 2, 0.5))
})
dev.off()

###Get fraction of SV class with functional consequences by frequency bin
all.fracs <- count.fracLOFSV(dat[which(!(dat$chrom %in% c("X", "Y"))), ])
common.fracs <- count.fracLOFSV(dat[which(dat$AF>=0.01 & !(dat$chrom %in% c("X", "Y"))), ])
rare.fracs <- count.fracLOFSV(dat[which(dat$AF<0.01 & dat$AC>1 & !(dat$chrom %in% c("X", "Y"))), ])
singleton.fracs <- count.fracLOFSV(dat[which(dat$AC==1 & !(dat$chrom %in% c("X", "Y"))), ])

#Alessia figures
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(qqman)
library(readxl)
library(ggpubr)

graphic.settings <- theme_classic(base_size = 14) + theme(axis.ticks = element_line(size = 0.3)) + 
  theme(legend.position="none") +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.subtitle=element_text(size=9), plot.title=element_text(size=18))


homo.sites.alessia<-homo.sites
homo.sites.alessia$HOM_TOT<-homo.sites.alessia$all.DEL + homo.sites.alessia$all.DUP + homo.sites.alessia$all.INV
homo.sites.alessia<-merge(homo.sites.alessia, pops, by.x="pop", by.y="pop")
all.sites.alessia<-all.sites
all.sites.alessia$HET_TOT<-all.sites.alessia$all.DEL + all.sites.alessia$all.DUP + all.sites.alessia$all.INV
all.sites.alessia<-merge(all.sites.alessia, pops, by.x="pop", by.y="pop")

rare.sites.alessia<-rare.sites
rare.sites.alessia$RARE_TOTAL<-rare.sites.alessia$all.DEL + rare.sites.alessia$all.DUP + rare.sites.alessia$all.INV
rare.sites.alessia<-merge(rare.sites.alessia, pops, by.x="pop", by.y="pop")

singleton.sites.alessia<-singleton.sites
singleton.sites.alessia$SINGLETON_TOTAL<-singleton.sites.alessia$all.DEL + singleton.sites.alessia$all.DUP + singleton.sites.alessia$all.INV
singleton.sites.alessia<-merge(singleton.sites.alessia, pops, by.x="pop", by.y="pop")


my_comparisons <- list( c("QGP_PAR", "QGP_ADM"),  c("QGP_PAR", "QGP_AFR"),  c("QGP_PAR", "QGP_GAR"), c("QGP_PAR", "QGP_SAS"),  c("QGP_PAR", "QGP_WEP"))
p1 <- ggplot( homo.sites.alessia, aes(x=pop, y=HOM_TOT, fill=color)) + geom_boxplot(outlier.shape = NA) + labs(y="Homozygous burden", x="", title="", subtitle="") + 
  graphic.settings +
  scale_fill_manual(values=c("#4E79A7","#76B7B2","#941594","#E15759","#EDC948","#F28E2B")) + 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label="p.format")

pdf(paste(OUTDIR,"/",prefix,".2h.pdf",sep=""),
    width=7.6,height=4.3)
p1
dev.off()

my_comparisons <- list( c("QGP_AFR", "QGP_ADM"), c("QGP_AFR", "QGP_PAR"), c("QGP_AFR", "QGP_GAR"), c("QGP_AFR", "QGP_SAS"), c("QGP_AFR", "QGP_WEP"))
p2 <- ggplot(all.sites.alessia, aes(x=pop, y=HET_TOT, fill=color)) + geom_boxplot(outlier.shape = NA) + labs(y="Heterozygous burden", x="", title="", subtitle="") + 
  graphic.settings + 
  scale_fill_manual(values=c("#4E79A7","#76B7B2","#941594","#E15759","#EDC948","#F28E2B")) + 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label="p.format")

pdf(paste(OUTDIR,"/",prefix,".2i.pdf",sep=""),
    width=7.6,height=4.3)
p2
dev.off()

library(dplyr)
library(ggpubr)
# ---- helpers ----
fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 2.2e-16, "< 2.2e-16", formatC(p, format="e", digits=2)))
}

# ---- Figure 2h: HOM vs PAR ----
cmp_hom <- c("QGP_ADM","QGP_AFR","QGP_GAR","QGP_SAS","QGP_WEP")
res_hom <- compare_means(
  HOM_TOT ~ pop, data = homo.sites.alessia,
  method = "wilcox.test", ref.group = "QGP_PAR", exact = FALSE
) %>%
  filter(group2 %in% cmp_hom) %>%
  mutate(p_exact = fmt_p(p)) %>%
  arrange(factor(group2, levels = cmp_hom)) %>%
  select(group1, group2, p_exact)

# ---- Figure 2i: HET vs AFR ----
cmp_het <- c("QGP_ADM","QGP_PAR","QGP_GAR","QGP_SAS","QGP_WEP")
res_het <- compare_means(
  HET_TOT ~ pop, data = all.sites.alessia,
  method = "wilcox.test", ref.group = "QGP_AFR", exact = FALSE
) %>%
  filter(group2 %in% cmp_het) %>%
  mutate(p_exact = fmt_p(p)) %>%
  arrange(factor(group2, levels = cmp_het)) %>%
  select(group1, group2, p_exact)

# ---- Print legend-ready text to console ----
cat("Fig. 2h (Homozygous deletions): two-sided Wilcoxon rank-sum vs Peninsular Arabs (PAR)\n")
apply(res_hom, 1, function(r) cat(sprintf("PAR vs %s: p = %s\n", r[['group2']], r[['p_exact']])))

cat("\nFig. 2i (Heterozygous deletions): two-sided Wilcoxon rank-sum vs African Arabs (AFR)\n")
apply(res_het, 1, function(r) cat(sprintf("AFR vs %s: p = %s\n", r[['group2']], r[['p_exact']])))

# ---- Also save tidy CSVs you can paste from ----
write.csv(res_hom, file = file.path(OUTDIR, paste0(prefix, ".2h_pvalues.csv")), row.names = FALSE)
write.csv(res_het, file = file.path(OUTDIR, paste0(prefix, ".2i_pvalues.csv")), row.names = FALSE)




rare<-ggplot(rare.sites.alessia, aes(x=pop, y=RARE_TOTAL, fill=color)) + geom_boxplot(outlier.shape = NA) + labs(y="Rare burden", x="", title="", subtitle="") + 
  graphic.settings + 
  scale_fill_manual(values=c("#4E79A7","#76B7B2","#941594","#E15759","#EDC948","#F28E2B")) + 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label="p.format") 
  
singleton<-ggplot(singleton.sites.alessia, aes(x=pop, y=SINGLETON_TOTAL, fill=color)) + geom_boxplot(outlier.shape = NA) + labs(y="Singleton burden", x="", title="", subtitle="") + 
  graphic.settings + 
  scale_fill_manual(values=c("#4E79A7","#76B7B2","#941594","#E15759","#EDC948","#F28E2B")) + 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label="p.format",label.y = c(50,65,80,95,110,125))



print(grid.arrange(p1, p2,singleton,rare, ncol=2))


total_size <-read.table("data/total_size.tsv")
names(total_size) <- c("sample", "TOT_HET_DEL", "TOT_HOM_DEL", "TOT_HET_DUP", "TOT_HOM_DUP", "TOT_HET_INV", "TOT_HOM_INV")  

total_size <- merge(total_size, qgp_pop_mapping_final,by.x = "sample" , by.y = "NAME")

total_size$HOM_TOTAL<-total_size$TOT_HOM_DEL+total_size$TOT_HOM_DUP+total_size$TOT_HOM_INV
total_size$HET_TOTAL<-total_size$TOT_HET_DEL+total_size$TOT_HET_DUP+total_size$TOT_HET_INV

p3 <- ggplot(total_size, aes(x=V2.y, y=HOM_TOTAL, fill=color)) + geom_boxplot(outlier.shape = NA) + labs(y="Distribution of total size of homozygous SVs per sample", x="", title="", subtitle="") + 
  graphic.settings +
  scale_fill_manual(values=c("#4E79A7","#76B7B2","#941594","#E15759","#EDC948","#F28E2B")) + 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label="p.format")

my_comparisons <- list( c("EastAfricans", "AdmixedPop"), c("EastAfricans", "CommonArabs"), c("EastAfricans", "PeninsularArabs"), c("EastAfricans", "South Asians"), c("EastAfricans", "Persians"))
p4 <- ggplot(total_size, aes(x=V2.y, y=HET_TOTAL, fill=color)) + geom_boxplot(outlier.shape = NA) + labs(y="Distribution of total size of heterozygous SVs per sample", x="", title="", subtitle="") + 
  graphic.settings + 
  scale_fill_manual(values=c("#4E79A7","#76B7B2","#941594","#E15759","#EDC948","#F28E2B")) + 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label="p.format")
print(grid.arrange(p3, p4, ncol=2))





dev.off()


#941594
#F28E2B
#4E79A7
#76B7B2
#EDC948
#E15759



