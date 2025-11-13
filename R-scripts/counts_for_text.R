library(tidyr)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)

QGP_hq.tsv.processed <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq.tsv.processed.tsv")
qgp_subpop_freq <-read.delim("~/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq_allele_freq_subpop.tsv")
qgp_tra_allele_freq_subpop <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq_allele_freq_subpop.tra.tsv")
qgp_tra_allele_freq_subpop = subset(qgp_tra_allele_freq_subpop, select = -c(ALT) )

qgp_ins_allele_freq_subpop <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq_allele_freq_subpop.ins.tsv")
AF_merged<-rbind(qgp_subpop_freq,qgp_tra_allele_freq_subpop,qgp_ins_allele_freq_subpop)
QGP_hq.tsv.processed<-merge(QGP_hq.tsv.processed, AF_merged, by.x="ID", by.y="ID")


QGP_hq.tsv.processed$N_HOMALT<-QGP_hq.tsv.processed$QGP_ADM_HOMO + QGP_hq.tsv.processed$QGP_AFR_HOMO + QGP_hq.tsv.processed$QGP_GAR_HOMO + QGP_hq.tsv.processed$QGP_PAR_HOMO + QGP_hq.tsv.processed$QGP_SAS_HOMO + QGP_hq.tsv.processed$QGP_WEP_HOMO

67276+13377+5396+36156+20628

nrow(filter(QGP_hq.tsv.processed, AF < 0.01, NSAMP > 1, Exonic == "true",SV.type=="DEL",IMH_AF < 0.01,X1000g_max_AF < 0.01,GD_POPMAX_AF<0.01,DGV_LOSS_Frequency <0.01,DISEASE_GENES=="true"))

nrow(filter(QGP_hq.tsv.processed, Exonic == "true",SV.type=="DEL",NSAMP > 1,AF < 0.01,IMH_AF < 0.01,X1000g_max_AF < 0.01,GD_POPMAX_AF<0.01,DGV_LOSS_Frequency <0.01))

singleton_count<-nrow(filter(QGP_hq.tsv.processed, NSAMP == 1))/nrow(QGP_hq.tsv.processed)

QGP_hq.tsv.processed$LENGTH<-QGP_hq.tsv.processed$SV.end-QGP_hq.tsv.processed$SV.start

QGP_hq.tsv.processed %>%

  summarize(MeanLength = mean(LENGTH, na.rm = TRUE))

result <- QGP_hq.tsv.processed %>%
  summarise(standard_deviation = sd(LENGTH, na.rm = TRUE))

result <- QGP_hq.tsv.processed %>%
  summarise(interquartile_range = IQR(LENGTH, na.rm = TRUE))

QGP_hq.tsv.processed.deldupinv %>%
  summarize(MeanLength = mean(LENGTH, na.rm = TRUE))

hist(QGP_hq.tsv.processed$LENGTH,xlim = c(0,10000))

QGP_hq.tsv.processed %>%
  group_by(SV.type) %>%
  summarize(MeanLength = max(LENGTH, na.rm = TRUE))

#Biological impact
QGP_hq.tsv.processed.deldupinv<-QGP_hq.tsv.processed[which(QGP_hq.tsv.processed$SV.type=="DEL" | QGP_hq.tsv.processed$SV.type=="DUP" |  QGP_hq.tsv.processed$SV.type=="INV"),]
QGP_hq.tsv.processed.deldupinv$LENGTH<-QGP_hq.tsv.processed.deldupinv$SV.end-QGP_hq.tsv.processed.deldupinv$SV.start

QGP_hq.tsv.processed %>%
  filter(SV.type=="DEL",NSAMP==1,N_HOMALT > 0,highest_PLI > 0.95) %>%
  group_by(Genic,Exonic) %>%
  tally()

sd(all.sites$total)

average(all.sites$total)

median(all.sites$total)

table<-QGP_hq.tsv.processed.deldupinv %>%
  filter(SV.type=="DEL") %>%
  group_by(SV.type,DISEASE_GENES,Exonic,highest_PLI > 0.95, AF < 0.01) %>%
  tally()



QGP_hq.tsv.processed.deldupinv %>%
  filter(grepl('ENC', ENCODEexperiments)) %>%
  group_by(SV.type) %>%
  tally()

QGP_hq.tsv.processed.deldupinv %>%
  filter(grepl('MIR', Gene.name)) %>%
  group_by(SV.type) %>%
  tally()

QGP_hq.tsv.processed.deldupinv %>%
  filter(grepl('txStart-txEnd', exon_merged)) %>%
  group_by(SV.type) %>%
  tally()

QGP_hq.tsv.processed.deldupinv %>%
  group_by(SV.type,Genic) %>% tally()
  summarize(MeanLength = max(LENGTH, na.rm = TRUE))

results<-QGP_hq.tsv.processed %>%
  group_by(SV.type,Exonic,Novel,NSAMP=1) %>%
  tally()

results<-QGP_hq.tsv.processed %>%
  group_by(SV.type,Exonic,Novel,AF>0.01) %>%
  tally()

results<-QGP_hq.tsv.processed %>%
  group_by(Exonic,Novel,AF > 0.01) %>%
  tally()

results<-QGP_hq.tsv.processed %>%
  filter(Novel=="true",SV.type!="TRA") %>%
  tally()

table(QGP_hq.tsv.processed.deldupinv$Genic)

QGP_hq.tsv.processed.deldupinv %>%
  group_by(SV.type,Exonic,N_HOMALT > 0) %>%
  tally()

counts<-QGP_hq.tsv.processed.deldupinv %>%
  group_by(SV.type,Exonic,Novel,DISEASE_GENES,N_HOMALT > 0)%>%
tally()

hom<-filter(QGP_hq.tsv.processed.deldupinv,SV.type == "DEL",Exonic == "true",DISEASE_GENES == "true",N_HOMALT > 0)
novel_hom<-filter(QGP_hq.tsv.processed.deldupinv,SV.type == "DEL",Exonic == "true",Novel == "true",DISEASE_GENES == "true",N_HOMALT > 0)

#Biological impact
disease_causing<-filter(QGP_hq.tsv.processed, AF < 0.01,
                        highest_PLI > 0.8, 
                        NSAMP > 1, 
                        Exonic == "true",
                        SV.type!="INV",
                        SV.type!="INS",
                        SV.type!="TRA",
                        IMH_AF < 0.01,
                        X1000g_max_AF < 0.01,
                        GD_POPMAX_AF<0.01,
                        DGV_LOSS_Frequency <0.01,
                        DISEASE_GENES=="true")

disease_causing$SIZE<-disease_causing$SV.end - disease_causing$SV.start

disease_causing_subset<-
  select(
    disease_causing,
    SV.chrom,
    SV.start,
    SV.end,
    SV.type,
    SIZE,
    Gene.name,
    AnnotSV.ranking,
    phenotypes_merged,
    dbVar_event,
    dbVar_variant,
    dbVar_status,
    AF,
    NSAMP,
    N_HOMALT,
    exon_merged,
    highest_PLI,
    Novel
  )

sites.counts.all <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq.total_counts.txt")
sites.counts.ins <- read.delim("/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq.total_counts.ins.txt")
sites.counts.tra <- read.delim("/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq.total_counts.tra.txt")
sites.counts.tra<-sites.counts.tra[which(sites.counts.tra$svtype == "TRA"),]


sites.counts.all.wide <- spread(sites.counts.all, svtype, count)
sites.counts.ins.wide <- spread(sites.counts.ins, svtype, count)
sites.counts.tra.wide <- spread(sites.counts.tra, svtype, count)
sites.counts.all.wide<-merge(sites.counts.all.wide,sites.counts.ins.wide,by = "sample")
sites.counts.all.wide<-merge(sites.counts.all.wide,sites.counts.tra.wide,by = "sample")

sites.counts.all.wide$total<-sites.counts.all.wide$DEL+
  sites.counts.all.wide$DUP+
  sites.counts.all.wide$INV+
  sites.counts.all.wide$INS+
  sites.counts.all.wide$TRA

min(sites.counts.all.wide$DEL)
max(sites.counts.all.wide$DEL)

apply(sites.counts.all.wide,2,min)
apply(sites.counts.all.wide,2,max)

sites.counts.all.wide
median(sites.counts.all.wide$total)
sd(sites.counts.all.wide$total)



