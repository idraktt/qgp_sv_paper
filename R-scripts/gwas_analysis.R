library(ggplot2)
library(hrbrthemes)
library(readxl)
library(dplyr)
library(tidyr)

##NEw results - Data generation
finalmerged_withPCs_FILTERED_BoneFerroni <- read.delim("C:/Users/ealiyev/Dropbox/SV Project - Data/4 - Analysis/GWAS/finalmerged_INVERSE_NORMALZIED.withPCs_BoneFerroni.tsv")

snps <- read.delim("C:/Users/ealiyev/Dropbox/SV Project - Data/4 - Analysis/SNP-Tagging/snps_tagged.tsv")


finalmerged_withPCs_FILTERED_BoneFerroni$ID<-lapply(strsplit(finalmerged_withPCs_FILTERED_BoneFerroni$MARKER_ID, "_"), tail, n=4)

finalmerged_withPCs_FILTERED_BoneFerroni$ID<-as.character(lapply(finalmerged_withPCs_FILTERED_BoneFerroni$ID, paste, collapse = "_"))

QGP_hq.tsv.processed <- read.delim("C:/Users/ealiyev/Desktop/QGP_annotated_master.tsv")

QGP_merged_with_association<-merge(x = finalmerged_withPCs_FILTERED_BoneFerroni, 
                                   y = QGP_hq.tsv.processed, 
                                   by.x="ID", by.y="AnnotSV.ID",all.x = TRUE)

pheno_list<-as.data.frame(unique(QGP_merged_with_association$PHENO))

karsten_najeeb_mapping <- read.delim("C:/Users/ealiyev/Dropbox/SV Project - Data/4 - Analysis/GWAS/karsten_najeeb_mapping.tsv")
write.table(pheno_list,file="C:/Users/ealiyev/Dropbox/SV Project - Data/4 - Analysis/GWAS/pheno.tsv",quote = FALSE,row.names = FALSE,sep = "\t" )

QGP_merged_with_association_karsten <- QGP_merged_with_association %>% 
  filter(PHENO %in% karsten_najeeb_mapping$Najeeb_pheno)  %>%
  filter(PHENO != "bicarbonate") %>%
  filter(PHENO != "c_peptide") 

count_by_sv_type <- QGP_merged_with_association_karsten %>%
  group_by(SV.type) %>%
  summarize(unique_ID_count = n_distinct(ID))

p_value<-0.05/(length(unique(QGP_merged_with_association$ID))*42)

QGP_merged_with_association_karsten$SV.length<-QGP_merged_with_association_karsten$SV.end-QGP_merged_with_association_karsten$SV.start

p_value<-0.05/(length(unique(QGP_merged_with_association_karsten[which(QGP_merged_with_association_karsten$SV.length < 2000000),]$ID))*42)

p_value<-0.05/(length(unique(QGP_merged_with_association_karsten[which(QGP_merged_with_association_karsten$SV.length < 2000000 & QGP_merged_with_association_karsten$AF > 0.01),]$ID))*42)


count_by_sv_type <- QGP_merged_with_association_karsten[which(QGP_merged_with_association_karsten$SV.length < 2000000),] %>%
  group_by(SV.type) %>%
  summarize(unique_ID_count = n_distinct(ID))

count_by_sv_type <- QGP_merged_with_association_karsten[which(QGP_merged_with_association_karsten$SV.length < 2000000),] %>%
  group_by(SV.type) %>%
  summarize(unique_ID_count = n_distinct(ID))


QGP_merged_with_association_significant<-QGP_merged_with_association_karsten %>% 
  filter(PVALUE < p_value)

QGP_merged_with_association_significant %>%
  group_by(PHENO) %>%
  tally()

QGP_merged_with_association_significant_bonferroni<-QGP_merged_with_association_karsten %>% 
  filter(bonf_PVALUE < p_value)



QGP_merged_with_association_significant_length_filtered<-QGP_merged_with_association_karsten %>% 
  filter(PVALUE < p_value & SV.length < 2000000)

QGP_merged_with_association_significant_bonferroni_length_filtered<-QGP_merged_with_association_karsten %>% 
  filter(bonf_PVALUE < 0.05 & SV.length < 2000000)

length(unique(QGP_merged_with_association_significant_length_filtered$PHENO))

final_table<-select(QGP_merged_with_association_significant_length_filtered,
                    ID,SV.chrom,SV.start,SV.end,SV.length,SV.type,PHENO,PVALUE,Gene.name,exon_merged,phenotypes_merged,
                    GENOCNT,MAF,STAT,BETA,SEBETA,AF,
                    ADM_AF,AFR_AF,GAR_AF,PAR_AF,SAS_AF,WEP_AF,
                    GLOBAL_MAX_DEL,GLOBAL_MAX_DUP,GLOBAL_MAX_INV,GLOBAL_MAX_INS)

sorted_snp <- snps %>%
  arrange(ID, desc(R2))

top_snps <- sorted_snp %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID,SNP_B,R2)

final_df <- final_table %>%
  left_join(top_snps, by = "ID")

write.table(final_df,file="C:/Users/ealiyev/Dropbox/SV Project - Data/4 - Analysis/GWAS/pub_table.tsv",quote = FALSE,row.names = FALSE,sep = "\t" )




  


