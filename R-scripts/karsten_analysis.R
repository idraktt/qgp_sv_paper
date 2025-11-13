library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)

QGP_hq <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq.tsv")
Soma_position <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QGP-RES-PUB-003/QGP_shared/Aziz_Elbay/Soma_position.txt")
QGP_hq_allele_freq_subpop <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq_allele_freq_subpop.tsv")

all_sv_annotated<-merge(QGP_hq, QGP_hq_allele_freq_subpop, by.x="ID", by.y="ID")

all_sv_annotated$N_HOMALT<-all_sv_annotated$QGP_ADM_HOMO + all_sv_annotated$QGP_AFR_HOMO + all_sv_annotated$QGP_GAR_HOMO + all_sv_annotated$QGP_PAR_HOMO + all_sv_annotated$QGP_SAS_HOMO + all_sv_annotated$QGP_WEP_HOMO

all_sv_annotated_sep<-separate(data = all_sv_annotated, col = location, into = c("loc1", "loc2"), sep = "-",remove = FALSE)

all_sv_annotated_sep_hom<-filter(all_sv_annotated_sep,N_HOMALT > 0,SV.type == "DEL",AnnotSV.type ==  "split")

all_sv_annotated_sep_hom_in_gene_list_exonic<-filter(all_sv_annotated_sep_hom,grepl("ex|tx",location))

hom_intron<-filter(all_sv_annotated_sep_hom,grepl("intron",loc1),grepl("intron",loc2))

hom_ex_dif_introns<-filter(hom_intron,loc1!=loc2)

exonic_sv<-rbind(all_sv_annotated_sep_hom_in_gene_list_exonic,hom_ex_dif_introns)

exonic_in_soma<-filter(exonic_sv,exonic_sv$Gene.name %in% Soma_position$NAME)



all_sv_annotated_sep_hom_in_gene_list<-filter(all_sv_annotated_sep_hom,all_sv_annotated_sep_hom$Gene.name %in% Soma_position$NAME)

all_sv_annotated_sep_hom_in_gene_list_exonic<-filter(all_sv_annotated_sep_hom_in_gene_list,grepl("ex|tx",location))

hom_intron<-filter(all_sv_annotated_sep_hom_in_gene_list,grepl("intron",loc1),grepl("intron",loc2))

hom_ex_dif_introns<-filter(hom_intron,loc1!=loc2)

karsten_svs<-rbind(all_sv_annotated_sep_hom_in_gene_list_exonic,hom_ex_dif_introns)

write.table(karsten_svs,"/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Karsten_Proteomics/analysis_18_11_2021/annotated_sv_filtered.tsv",sep="\t",row.names = FALSE,quote = FALSE)
write.table(karsten_svs$ID,"/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Karsten_Proteomics/analysis_18_11_2021/ids.list",sep="\t",row.names = FALSE,quote = FALSE)

length(unique(all_sv_annotated$ID))

all_sv_annotated_sep_hom<-filter(all_sv_annotated_sep,SV.type == "DEL",N_HOMALT > 0,AnnotSV.type ==  "split")

length(unique(exonic_sv$ID))


