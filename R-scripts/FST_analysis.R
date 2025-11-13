library(dplyr)
library(tidyr)
library(CMplot)

install.packages("CMplot")

mergedFSTcalculation <- read.delim("C:/Users/ealiyev/Dropbox/SV Project - Data/4 - Analysis/FST/mergedFSTcalculation.tsv")

mergedFSTcalculation <- mergedFSTcalculation %>%
  distinct(ID, .keep_all = TRUE)

df <- mergedFSTcalculation %>%
  separate(ID, into = c("CHROM", "START","END","SVTYPE","POPS"), sep = "_", extra = "merge")

df$ID<-paste(df$CHROM,df$START,df$END,df$SVTYPE,sep = "_")

df <- df %>% select(-c(X.CHROM, POS))

df_filtered<-df %>%
  filter(HUDSON_FST > 0.2)

QGP_hq.tsv.processed <- read.delim("C:/Users/ealiyev/Dropbox/SV Project - Data/QGP_annotated_master.tsv")

df_filtered_annotated<-merge(x = df_filtered,
                                   y = QGP_hq.tsv.processed,
                                   by.x="ID", by.y="AnnotSV.ID",all.x = TRUE)

df_filtered_plof<-df_filtered_annotated %>%
  filter( PREDICTED_LOF_count> 0)

write.table(df_filtered_annotated,file="C:/Users/ealiyev/Dropbox/SV Project - Data/4 - Analysis/FST/fst_higher0.2_annotated.tsv",quote = FALSE,row.names = FALSE,sep = "\t" )


# Reshape the dataframe
df_new <- df %>%
  spread(POPS, HUDSON_FST) %>% group_by(ID) %>%
  summarise(across(starts_with("ADM.AFR"):starts_with("SAS.WEP"), sum, na.rm = TRUE, .names = "sum_{.col}")) %>%
  ungroup()

df_cm <- df_new %>%
  separate(ID, into = c("CHROM", "START","END","SVTYPE"), sep = "_", remove = FALSE)

df_cm<-df_cm %>%  select(ID,CHROM,START,sum_GAR.PAR,sum_PAR.WEP,sum_PAR.SAS,sum_AFR.PAR)

CMplot(df_cm, plot.type="c", LOG10=FALSE, cir.chr.h=1.3,chr.den.col="black",r=0.4,cir.axis=TRUE,
       outward=FALSE,cir.axis.col="black", trait.legend.ncol=4, trait.legend.pos="left",file.output = FALSE)


QGP_SV.annotated<-left_join(QGP_SV.annotated, df_cm, by = c("ID" = "ID"))

df_merged <- left_join(df, QGP_SV.annotated, by = c("ID" = "ID"))

test<-df_merged %>%
  filter(DISEASE_GENES == "TRUE" & HUDSON_FST > 0.2 )%>%
  tally()

library(igraph)
library(dplyr)
setwd("C:/Dropbox/SVAssociation/fstCalculation")
file_filtered = read.csv("FST_filtered.fst.summary",  sep = "\t", header = TRUE)
colnames(file_filtered) = c("POP1","POP2","FST")
file_filtered = file_filtered %>% filter(POP1 != "QGP_EUR" ) %>% filter(POP2 != "QGP_EUR") %>%
  filter(POP1 != "QGP_ADM" ) %>% filter(POP2 != "QGP_ADM")
g_fil <- graph.data.frame(file_filtered, directed=FALSE)
# add value as a weight attribute
m_fil <- get.adjacency(g_fil, attr="FST", sparse=FALSE)
hc <- hclust(dist(m_fil),method="complete" )
png("Figure3.png", res = 600, height = 5, width = 5, units = "in")
ph <- pheatmap::pheatmap(m_fil,
                         display_numbers = TRUE,  # Add numbers to the heatmap cells
                         number_format = "%.2f",  # Format numbers to 2 decimal places
                         fontsize_number = 8)
# m_fil[lower.tri(m_fil)] <- NA
# pheatmap::pheatmap(m_fil,cluster_rows=F, cluster_cols=T)
dev.off()

pdf(paste( "Figure3", ".pdf", sep=""), 
    height=5, width=5)
print(ph)
dev.off()

library(grid)



file_filtered = read.csv("plink2.fst.summary",  sep = "\t", header = TRUE)
colnames(file_filtered) = c("POP1","POP2","FST")
file_filtered = file_filtered %>% filter(POP1 != "QGP_EUR" ) %>% filter(POP2 != "QGP_EUR") %>%
  filter(POP1 != "QGP_ADM" ) %>% filter(POP2 != "QGP_ADM")
g_fil <- graph.data.frame(file_filtered, directed=FALSE)
# add value as a weight attribute
m_fil <- get.adjacency(g_fil, attr="FST", sparse=FALSE)
hc <- hclust(dist(m_fil),method="complete" )
png("Figure1.png", res = 600, height = 6, width = 8, units = "in")
ph <- pheatmap::pheatmap(m_fil,
                         cluster_rows = FALSE,    # Remove row dendrogram
                         cluster_cols = FALSE,    # Remove column dendrogram
                         display_numbers = TRUE,  # Add numbers to the heatmap cells
                         number_format = "%.2f",  # Format numbers to 2 decimal places
                         fontsize_number = 10)
# m_fil[lower.tri(m_fil)] <- NA
# pheatmap::pheatmap(m_fil,cluster_rows=F, cluster_cols=T)
dev.off()



