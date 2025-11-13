library(tidyr)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(tidyverse)

setwd("C:/Dropbox/SV Project - Paper Submission/Figure_generation")

QGP_hq.tsv.processed <- read.delim("data/QGP_hq.tsv.processed.tsv")
QGP_hq.tsv.processed_no_tra <- QGP_hq.tsv.processed[QGP_hq.tsv.processed$SV.type != "TRA", ]
QGP_hq.tsv.processed_no_tra <- QGP_hq.tsv.processed_no_tra[!duplicated(QGP_hq.tsv.processed_no_tra$AnnotSV.ID), ]
table(QGP_hq.tsv.processed_no_tra$SV.type)

QGP_hq.tsv.processed_tra <- QGP_hq.tsv.processed[QGP_hq.tsv.processed$SV.type == "TRA", ]
QGP_hq.tsv.processed_tra$AnnotSV.ID<-paste(QGP_hq.tsv.processed_tra$ID,"_",QGP_hq.tsv.processed_tra$ALT,sep = "")
QGP_hq.tsv.processed_tra <- QGP_hq.tsv.processed_tra[!duplicated(QGP_hq.tsv.processed_tra$AnnotSV.ID), ]


qgp_tra_allele_freq_subpop <- read.delim("data/QGP_hq_allele_freq_subpop.tra.tsv")
qgp_tra_allele_freq_subpop$ID<-paste(qgp_tra_allele_freq_subpop$ID,"_",qgp_tra_allele_freq_subpop$ALT,sep = "")
qgp_tra_allele_freq_subpop = subset(qgp_tra_allele_freq_subpop, select = -c(ALT) )
qgp_tra_allele_freq_subpop <- qgp_tra_allele_freq_subpop[order(qgp_tra_allele_freq_subpop$ID, -abs(qgp_tra_allele_freq_subpop$AF) ), ] 
nrow(qgp_tra_allele_freq_subpop)
qgp_tra_allele_freq_subpop<-qgp_tra_allele_freq_subpop[!duplicated(qgp_tra_allele_freq_subpop$ID), ]
nrow(qgp_tra_allele_freq_subpop)


qgp_subpop_freq <-read.delim("data/QGP_hq_allele_freq_subpop.tsv")
qgp_subpop_freq <- qgp_subpop_freq[order(qgp_subpop_freq$ID, -abs(qgp_subpop_freq$AF) ), ] 
nrow(qgp_subpop_freq)
qgp_subpop_freq<-qgp_subpop_freq[!duplicated(qgp_subpop_freq$ID), ]
nrow(qgp_subpop_freq)
qgp_ins_allele_freq_subpop <- read.delim("data/QGP_hq_allele_freq_subpop.ins.tsv")
qgp_ins_allele_freq_subpop <- qgp_ins_allele_freq_subpop[order(qgp_ins_allele_freq_subpop$ID, -abs(qgp_ins_allele_freq_subpop$AF) ), ] 
nrow(qgp_ins_allele_freq_subpop)
qgp_ins_allele_freq_subpop<-qgp_ins_allele_freq_subpop[!duplicated(qgp_ins_allele_freq_subpop$ID), ]
nrow(qgp_ins_allele_freq_subpop)
AF_merged<-rbind(qgp_subpop_freq,qgp_tra_allele_freq_subpop,qgp_ins_allele_freq_subpop)
QGP_hq.tsv.processed<-rbind(QGP_hq.tsv.processed_tra,QGP_hq.tsv.processed_no_tra)
QGP_hq.tsv.processed<-merge(QGP_hq.tsv.processed, AF_merged, by.x="AnnotSV.ID", by.y="ID")

QGP_hq.tsv.processed$AN<-QGP_hq.tsv.processed$NSAMP
QGP_hq.tsv.processed$AC<-QGP_hq.tsv.processed$QGP_ADM_AC + QGP_hq.tsv.processed$QGP_AFR_AC + QGP_hq.tsv.processed$QGP_WEP_AC + QGP_hq.tsv.processed$QGP_GAR_AC + QGP_hq.tsv.processed$QGP_PAR_AC + QGP_hq.tsv.processed$QGP_SAS_AC

dat.all<-QGP_hq.tsv.processed %>% select(AnnotSV.ID, SV.chrom, SV.start, SV.end, SV.length, SV.type,AN,AF,AC)
dat.all$SVLEN<-dat.all$SV.length
dat.all[which(dat.all$SV.type == "DEL"),]$SVLEN<-dat.all[which(dat.all$SV.type == "DEL"),]$SV.end-dat.all[which(dat.all$SV.type == "DEL"),]$SV.start

dat.all$SVTYPE<-dat.all$SV.type


table<-QGP_hq.tsv.processed_no_tra %>%
  filter(SV.type=="DEL") %>%
  group_by(SV.type,DISEASE_GENES,Exonic,highest_PLI > 0.95, AF < 0.01) %>%
  tally()

pops.file <- "data/QGP_population_colors.txt"
if(!is.null(pops.file)){
  pops <- read.table(pops.file,sep="\t",header=T,comment.char="",check.names=F)
}

svtypes.file <- "data/SV_colors.txt"
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

qgp_mapping <- read.delim("data/mapping.tsv")
Population_structure_labels <- read.delim("data//Population_structure_labels", header=FALSE)
qgp_population_structure_clusters <- read.delim("data//qgp_population_structure_clusters.txt", header=FALSE)
population<-merge(qgp_population_structure_clusters, Population_structure_labels, by.x="V2", by.y="V1")
population<-merge(population, pops, by.x="V2.y", by.y="name")
qgp_white_list <- read.table("data//qgp_bam_white.list", quote="\"", comment.char="")

white_list_samples<-merge(qgp_mapping,qgp_white_list,by.x = "BAM_LOC",by.y = "V1")

qgp_pop_mapping_final<-merge(white_list_samples,population,by.x = "NAME", by.y = "V1")

OUTDIR <- "output"
prefix <- "QGP_SV"




###Plot simple bars of counts -- done
cat("NOW STARTING COUNT OF SV BY TYPE\n")
pdf(paste(OUTDIR, "/", prefix, ".site_counts_by_type_all.pdf", sep=""), 
    height=2.75, width=2.25)
plot.totalCounts(dat=dat.all, svtypes=svtypes, thousandG=F)
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
                 svtypes$color[which(svtypes$svtype=="TRA")])
site.labels <- c("DEL","MCNV (Loss)","DUP","MCNV (Gain)",
                 "INS","INV","TRA")

sites.counts.all <- read.delim("data/QGP_hq.total_counts.txt")
sites.counts.ins <- read.delim("data/QGP_hq.total_counts.ins.txt")
sites.counts.tra <- read.delim("data/QGP_hq.total_counts.tra.txt")
sites.counts.tra<-sites.counts.tra[which(sites.counts.tra$svtype == "TRA"),]


sites.counts.all.wide <- spread(sites.counts.all, svtype, count)
sites.counts.ins.wide <- spread(sites.counts.ins, svtype, count)
sites.counts.tra.wide <- spread(sites.counts.tra, svtype, count)
sites.counts.all.wide<-merge(sites.counts.all.wide,sites.counts.ins.wide,by = "sample")
sites.counts.all.wide<-merge(sites.counts.all.wide,sites.counts.tra.wide,by = "sample")

sites.counts.all.wide.pop<-merge(sites.counts.all.wide,qgp_pop_mapping_final,by.x = "sample",by.y = "SAMPLE_NAME")
all.sites<-setNames(data.frame(matrix(ncol = 9, nrow = 6116)), c("sample","pop",
                                                                 "all.DEL","all.MCNV_DEL",
                                                                 "all.DUP","all.MCNV_DUP",
                                                                 "all.INS","all.INV","all.TRA"))

all.sites$sample<-sites.counts.all.wide.pop$sample
all.sites$pop<-sites.counts.all.wide.pop$pop
all.sites$all.DEL<-sites.counts.all.wide.pop$DEL
all.sites$all.MCNV_DEL<-0
all.sites$all.DUP<-sites.counts.all.wide.pop$DUP
all.sites$all.MCNV_DUP<-0
all.sites$all.INS<-sites.counts.all.wide.pop$INS
all.sites$all.INV<-sites.counts.all.wide.pop$INV
all.sites$all.TRA<-sites.counts.all.wide.pop$TRA

all.sites.wdat <- getWaterfallData(dat=all.sites)

all.sites$total<-all.sites$all.DEL + all.sites$all.DUP + all.sites$all.INS + all.sites$all.TRA + all.sites$all.INV


#Generate smaller plot of just counts per sample for main figure
png(paste(OUTDIR,"/",prefix,".sites_perSample.waterfall.small.all.png",sep=""),
    width=10.5*300,height=1.5*300,res=400)
par(mar=c(0.5,3.5,1.1,5.5),bty="n")
plotWaterfall(wdat=all.sites.wdat,colors=site.colors,
              ylabel="SVs per Genome",titles=T,cat.labels=site.labels,y.at.nmax=10,cex.ylabs=0.6)
dev.off()

pdf(paste(OUTDIR,"/",prefix,".sites_perSample.waterfall.small.all.pdf",sep=""),
    width=10.5,height=1.5)
par(mar=c(0.5,3.5,1.1,5.5),bty="n")
plotWaterfall(wdat=all.sites.wdat,colors=site.colors,
              ylabel="SVs per Genome",titles=T,cat.labels=site.labels,y.at.nmax=10,cex.ylabs=0.6)
dev.off()



#Plot distribution of allele counts by SVTYPE
plot.freqByTypeLC <- function(dat, ymin=NULL, axlabel.cex=1){
  #Gather data
  exclude <- grep("MULTIALLELIC", dat$FILTER, fixed=T)
  if(length(exclude)>0){
    dat <- dat[-exclude, ]
  }
  dat <- dat[which(dat$SVTYPE!="MCNV" & !(dat$SV.chrom %in% c("X", "Y"))), ]
  AC.cutoffs <- c(1:9, seq(10, 90, 10), seq(100, 900, 100), seq(1000, 25000, 1000))
  AF.dat <- as.data.frame(sapply(c("DUP", "DEL", "INV","INS","TRA"), function(svtype){
    sapply(AC.cutoffs, function(AC){
      length(which(dat$AC[which(dat$SVTYPE==svtype)]<=AC))/length(which(dat$SVTYPE==svtype))
    })
  }))
  # AF.dat$OTH <- as.numeric(sapply(AC.cutoffs, function(AC){
  #   return(length(which(dat$AC[which(!(dat$SVTYPE %in% svtypes$svtype[which(svtypes$svtype!="OTH")]))]<=AC))/length(which(!(dat$SVTYPE %in% svtypes$svtype[which(svtypes$svtype!="OTH")]))))
  # }))
  #Plot
  if(is.null(ymin)){
    ymin <- log10(floor(100*min(AF.dat, na.rm=T))/100)
  }else{
    ymin <- log10(ymin)
  }
  AF.dat <- as.data.frame(apply(AF.dat, 2, log10))
  xrange <- c(-0.25, max(log10(AC.cutoffs)))
  common.threshold <- min(as.numeric(dat$AC[which(as.numeric(dat$AF)>=0.01 & as.numeric(dat$AN)==max(dat$AN, na.rm=T))]), na.rm=T)
  common.threshold <-120
  par(mar=c(2, 3.25, 0.5, 0.5), bty="n")
  plot(x=xrange, y=c(ymin, log10(1)), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="")
  # rect(xleft=-0.1, xright=0.1, ybottom=par("usr")[3], ytop=par("usr")[4], 
  #      border=NA, bty="n", col="gray90")
  # abline(v=log10(common.threshold), col="gray80")
  segments(x0=log10(common.threshold), x1=log10(common.threshold), 
           y0=par("usr")[3], y1=log10(1), col="gray80")
  text(x=log10(common.threshold)+(0.025*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])), 
       labels="Rare\n(AF<1%)", pos=2, cex=0.7)
  text(x=log10(common.threshold)-(0.025*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])), 
       labels="Common\n(AF>1%)", pos=4, cex=0.7)
  sapply(c("DUP", "DEL", "INV","INS","TRA"), function(svtype){
    points(x=log10(AC.cutoffs), y=AF.dat[, which(colnames(AF.dat)==svtype)], 
           col=svtypes$color[which(svtypes$svtype==svtype)], type="l", lwd=1.5)
  })
  sapply(c("DUP", "DEL", "INV","INS","TRA"), function(svtype){
    # points(x=log10(AC.cutoffs)[1], y=AF.dat[, which(colnames(AF.dat)==svtype)][1], 
    #        col="black", pch="-", cex=1.4)
    # points(x=log10(AC.cutoffs)[1], y=AF.dat[, which(colnames(AF.dat)==svtype)][1], 
    #        col=svtypes$color[which(svtypes$svtype==svtype)], pch="-", cex=1.2)
    rect(xleft=log10(AC.cutoffs)[1]-(0.03*(par("usr")[2]-par("usr")[1])), 
         xright=log10(AC.cutoffs)[1]+(0.03*(par("usr")[2]-par("usr")[1])), 
         ybottom=AF.dat[, which(colnames(AF.dat)==svtype)][1]-(0.01*(par("usr")[4]-par("usr")[3])), 
         ytop=AF.dat[, which(colnames(AF.dat)==svtype)][1]+(0.01*(par("usr")[4]-par("usr")[3])), 
         col=svtypes$color[which(svtypes$svtype==svtype)])
  })
  logscale.all <- log10(as.numeric(unlist(sapply(c(0:9), function(i){(1:9)*(10^i)}))))
  logscale.major <- 0:9
  axis(1, at=logscale.all, labels=NA, tck=-0.015, lwd=0.7)
  axis(1, at=logscale.major, labels=NA, tck=-0.03, lwd=1.1)
  sapply(1:5, function(i){
    axis(1, at=i-1, tick=F, line=-0.9, cex.axis=0.8, 
         labels=c("1", "10", "100", "1k", "10k")[i])
  })
  mtext(1, text="Allele Count", line=1, cex=axlabel.cex)
  logscale.pct.all <- log10((1:100)/100)
  logscale.pct.major <- log10(seq(10, 100, 10)/100)
  axis(2, at=logscale.pct.all, labels=NA, tck=-0.015, lwd=0.9)
  axis(2, at=logscale.pct.major, labels=NA, tck=-0.03, lwd=1.1)
  axis(2, at=logscale.pct.major, tick=F, line=-0.5, cex.axis=0.8, las=2, 
       labels=paste(seq(10, 100, 10), "%", sep=""))
  mtext(2, text="Fraction of SVs", line=2.25, cex=axlabel.cex)
}

###Plot frequency distribution by SVTYPE
cat("NOW STARTING FREQUENCY DISTRIBUTIONS\n")
pdf(paste(OUTDIR, "/", prefix, ".site_frequency_distributions_bySVtype_all.pdf", sep=""), 
    height=2.8, width=3.5)
plot.freqByTypeLC(dat=dat.all)
dev.off()

###Plot size distribution
cat("NOW STARTING SIZE DISTRIBUTIONS\n")
pdf(paste(OUTDIR, "/", prefix, ".size_distribution_by_type_all.pdf", sep=""), 
    height=2, width=2.9)
plot.sizes(dat=dat.all[which(dat.all$SV.type != "TRA"),], svtypes=svtypes)
dev.off()

table(dat.all$SV.type)

write.table(QGP_hq.tsv.processed, file='/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/publication_plots/data/QGP.tsv', quote=FALSE, sep='\t')
write.table(QGP_hq.tsv.processed, file='/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_Deliverables/QGP.tsv', quote=FALSE, sep='\t')



