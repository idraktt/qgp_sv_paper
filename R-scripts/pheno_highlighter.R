library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(tidyr)
library(readxl)
library(qqman)




1.825E-11
2.283E-06
1.024E-07
4.806E-08

##CD36


phenotypes_df <- data.frame(
  name = c(
    "White.Blood.Cell",
    "Red.Blood.Cell",
    "Mean.Cell.Hemoglobin.Concentration",
    "Mean.Cell.Hemoglobin",
    "Platelet",
    "Mean.Platelet.Volume",
    "Hemoglobin",
    "Hematocrit",
    "Neutrophil.Auto.."
  ),
  label = c(
    "White Blood Cell",
    "Red Blood Cell",
    "Mean Cell Hemoglobin Concentration",
    "Mean Cell Hemoglobin",
    "Platelet Count",
    "Mean Platelet Volume",
    "Hemoglobin",
    "Hematocrit",
    "Neutrophils"
  ),
  stringsAsFactors = FALSE
)
for (i in 1:nrow(phenotypes_df)) {
  df <- read.delim("data/CD36.tsv", header=FALSE)
  colnames(df)<-c("SAMPLE","ID","GT")
  merged_ranked_pheno <- read.delim("data/merged_ranked_pheno.tsv")
  merged_ranked_pheno[merged_ranked_pheno == -9] <- NA
  
  df<-merge(df,merged_ranked_pheno,by.x = "SAMPLE",by.y = "SAMPLE_NAME")
  df$GT <- factor(df$GT , levels=c("0/0", "0/1", "1/1"))
  
  my_comparisons <- list( c("0/0", "0/1"), c("0/0", "1/1"), c("0/1", "1/1") )
  pheno<-phenotypes_df$name[i]
  display_label <- phenotypes_df$label[i]
  
  kruskal.test(df[[pheno]] ~ df$GT, data = df)
  
  theme_set(theme_bw(base_size = 12))
  stat_box_data <- function(y, upper_limit = max(df[[pheno]]) * 1.15) {
    return( 
      data.frame(
        y = 0.95 * upper_limit,
        label = paste('count =', length(y), '\n',
                      'mean =', round(mean(y), 1), '\n')
      )
    )
  }
  max_value<-max(merged_ranked_pheno[[pheno]],na.rm = TRUE)
  
  p<-ggboxplot(df, x = "GT", y = paste(pheno),title = paste(""),conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
    stat_summary(
      fun.data = stat_box_data, 
      geom = "text" )+ stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = max_value*1.5) +     # global
    theme(
      axis.title.x = element_text(size = 11, face = "bold"),  # X-axis title
      axis.title.y = element_text(size = 11, face = "bold")  # Y-axis title
      #axis.text.x  = element_text(size = 12, face = "bold"),  # X-axis text (e.g., 0/0, 0/1)
      #axis.text.y  = element_text(size = 12, face = "bold")   # Y-axis tick labels
    ) +
    labs(y = display_label)
  print(p)
  
  ggsave(paste("Figure_4/CD36","_",pheno,".pdf",sep = ""), plot = p, width = 4, height = 5)
}
##GBP3

df <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/GBP3.tsv", header=FALSE)
colnames(df)<-c("SAMPLE","ID","GT")
merged_ranked_pheno <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/Phenotypes_6218/merged_ranked_pheno.tsv")
merged_ranked_pheno[merged_ranked_pheno == -9] <- NA

df<-merge(df,merged_ranked_pheno,by.x = "SAMPLE",by.y = "SAMPLE_NAME")
df$GT <- factor(df$GT , levels=c("0/0", "0/1", "1/1"))

my_comparisons <- list( c("0/0", "0/1"), c("0/0", "1/1"), c("0/1", "1/1") )

max_value<-max(merged_ranked_pheno[[pheno]],na.rm = TRUE)

p<-ggboxplot(df, x = "GT", y = paste(pheno),title = paste(""),conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text" )+ stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = max_value*1.5)     # Add global p-value
print(p)

##SLC2A9

df <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/SLC2A9.tsv", header=FALSE)
colnames(df)<-c("SAMPLE","ID","GT")
merged_ranked_pheno <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/Phenotypes_6218/merged_ranked_pheno.tsv")
merged_ranked_pheno[merged_ranked_pheno == -9] <- NA

df<-merge(df,merged_ranked_pheno,by.x = "SAMPLE",by.y = "SAMPLE_NAME")
df$GT <- factor(df$GT , levels=c("0/0", "0/1", "1/1"))

my_comparisons <- list( c("0/0", "0/1"), c("0/0", "1/1"), c("0/1", "1/1") )

max_value<-max(merged_ranked_pheno[[pheno]],na.rm = TRUE)

p<-ggboxplot(df, x = "GT", y = paste(pheno),title = paste(""),conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text" )+ stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = max_value*1.5)     # Add global p-value
print(p)

##TMEM175
phenotypes<-colnames(merged_ranked_pheno)

pheno<-data.frame(phenotypes[12:240])
names(pheno) <- c("pheno")  

pheno<-pheno %>%
  filter(across(everything(), ~ !grepl('Dummy', .))) %>%
  filter(across(everything(), ~ !grepl('Rank', .)))

df <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/CYTH2.tsv", header=FALSE)
colnames(df)<-c("SAMPLE","ID","GT")
merged_ranked_pheno <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/Phenotypes_6218/merged_ranked_pheno.tsv")
merged_ranked_pheno[merged_ranked_pheno == -9] <- NA

df<-merge(df,merged_ranked_pheno,by.x = "SAMPLE",by.y = "SAMPLE_NAME")
df$GT <- factor(df$GT , levels=c("0/0", "0/1", "1/1"))

my_comparisons <- list( c("0/0", "0/1"), c("0/0", "1/1"), c("0/1", "1/1") )

pdf("/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/QGP_hq.CYTH2.pdf")

for (row in 1:nrow(pheno)) {
  
  cur_pheno<-pheno[row, "pheno"]
  max_value<-max(merged_ranked_pheno[[cur_pheno]],na.rm = TRUE)
  

  p<-ggboxplot(df, x = "GT", y = paste(cur_pheno),title = paste(""),conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
    stat_summary(
      fun.data = stat_box_data, 
      geom = "text" )+ stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = max_value*1.5)     # Add global p-value
  print(p)
}
dev.off()


##Najeeb file
results <- read.delim("/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/QGP_hq.najeeb_result_no_prothrombin.tsv")
merged_ranked_pheno <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/Phenotypes_6218/merged_ranked_pheno.tsv")
merged_ranked_pheno[merged_ranked_pheno == -9] <- NA
plot_list = list()
list<-system2(command = "/gpfs/software/genomics/bcftools/1.9/bin/bcftools",args    = c("query", 
                                                                                        "-f","'[%SAMPLE\\t%CHROM/_%POS/_%END/_%INFO/SVTYPE\\t%GT\\n]'",
                                                                                        "/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/QGP_hq.managenese.vcf"), stdout = TRUE)
df <- data.frame(matrix(unlist(list), nrow=length(list), byrow=T))
df<-separate(data = df, col = matrix.unlist.list...nrow...length.list...byrow...T., into = c("SAMPLE", "ID","GT"), sep = "\\t")
distribution<-merge(df,merged_ranked_pheno,by.x = "SAMPLE",by.y = "SAMPLE_NAME.y")

distribution$ID <- gsub("/","",distribution$ID)

distribution$GT <- factor(distribution$GT , levels=c("0/0", "0/1", "1/1"))

pdf("/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/QGP_hq.najeeb_result_plots.pdf")

for (row in 1:nrow(results)) {
  id <- results[row, "ID"]
  pheno<-results[row, "PHENO_CORRECTED"]
  gene_name<-results[row, "Gene.name"]
  p_value<-results[row, "PVALUE"]
  p<-ggboxplot(distribution[which(distribution$ID == id),], x = "GT", y = paste(pheno),title = paste("ID:",id,"p-value:",p_value,"\n","Gene:",gene_name),conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
    stat_summary(
      fun.data = stat_box_data, 
      geom = "text", 
      hjust = 0.5,
      vjust = 0.9
    )
  print(p)
  
}
dev.off()


##Najeeb results end
top400 <- read_excel("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/Publication/PLINK/top400.xlsx", 
                     sheet = "nteresting")
merged_ranked_pheno <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/Phenotypes_6218/merged_ranked_pheno.tsv")
plot_list = list()
pdf("/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/PLINK_output/plots.pdf")

for (row in 1:nrow(top400)) {
  id <- top400[row, "ID"]
  pheno<-top400[row, "PHENO\\n"]
  gene_name<-top400[row, "Gene name"]
  p_value<-top400[row, "P...98"]
  pheno<-gsub("plink.","",pheno)
  pheno<-gsub(".qassoc","",pheno)
  list<-system2(command = "/gpfs/software/genomics/bcftools/1.9/bin/bcftools",args    = c("query", 
                                                                                          "-f","'[%SAMPLE\\t%GT\\n]'",
                                                                                          "/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/QGP_hq.mersad_filtered.vcf"), stdout = TRUE)
  df <- data.frame(matrix(unlist(list), nrow=length(list), byrow=T))
  df<-separate(data = df, col = matrix.unlist.list...nrow...length.list...byrow...T., into = c("SAMPLE", "GT"), sep = "\\t")
  distribution<-merge(df,merged_ranked_pheno,by.x = "SAMPLE",by.y = "SAMPLE_NAME.y")
  distribution$GT <- factor(distribution$GT , levels=c("0/0", "0/1", "1/1"))
  p<-ggboxplot(distribution, x = "GT", y = paste(pheno),title = paste("ID:",id,"p-value:",p_value,"Gene:",gene_name),conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
    stat_summary(
      fun.data = stat_box_data, 
      geom = "text", 
      hjust = 0.5,
      vjust = 0.9
    )
  print(p)
  
}
dev.off()


qq(plink.White.Blood.Cell$P)

nrow(plink.Activated.Partial.Thromboplastin.Time)


system2("module",args = c("load","bcftools/1.9"))

list<-system2(command = "/gpfs/software/genomics/bcftools/1.9/bin/bcftools",args    = c("query", 
                                                                                        "-f","'[%SAMPLE\\t%GT\\t%ID\\n]'",
                                                                                        "-i",'\'ID="DEL00105232SUR"\'',
                                                                                        "/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/QGP_SV/reports/population_final/QGP_hq.sorted.no_supp_vec.sorted.vcf.gz"), stdout = TRUE)
df <- data.frame(matrix(unlist(list), nrow=length(list), byrow=T))

separate(data = df, col = matrix.unlist.list...nrow...length.list...byrow...T., into = c("SAMPLE", "GT"), sep = "\\t")


#Plot distribution
DEL00105232SUR<-DEL00105232SUR %>% 
  rename(
    GT = V2
  )
distribution<-merge(DEL00105232SUR,merged_ranked_pheno,by.x = "V1",by.y = "SAMPLE_NAME.y")



merged_ranked_pheno <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/Phenotypes_6218/merged_ranked_pheno.tsv")
merged_ranked_pheno <- merged_ranked_pheno[which(merged_ranked_pheno$V2.y != "1European"),]
merged_ranked_pheno[which(merged_ranked_pheno$Iron < 0),]$Iron<-NA
merged_ranked_pheno[which(merged_ranked_pheno$Platelet < 0),]$Platelet<-NA
merged_ranked_pheno[which(merged_ranked_pheno$BMI < 0),]$BMI<-NA
merged_ranked_pheno[which(merged_ranked_pheno$BMI > 60),]$BMI<-NA


iron_df <- merged_ranked_pheno %>% 
  filter(QBB_ID=="SI000020000523" | QBB_ID=="SI000020000524")

samples_highlight<-c

platelet_df <-merged_ranked_pheno %>% 
  filter(QBB_ID=="SI000020100136")


merged_ranked_pheno %>%
  ggplot( aes(x=V2.y, y=BMI, fill=V2.y)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.4) +
  geom_jitter(data=platelet_df, aes(x=V2.y,y=BMI), color='red',size=3) +
  theme_gray() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Distribution of BMI levels in QGP subpopulations") +
  xlab("")

p<-ggboxplot(distribution, x = "GT", y = "Folate",title = "manganese",conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )
print(p)


##SPIRE

df <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/SPIRE2.tsv", header=FALSE)
colnames(df)<-c("SAMPLE","ID","GT")
merged_ranked_pheno <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/Phenotypes_6218/merged_ranked_pheno.tsv")
merged_ranked_pheno[merged_ranked_pheno == -9] <- NA

df<-merge(df,merged_ranked_pheno,by.x = "SAMPLE",by.y = "SAMPLE_NAME")
df$GT <- factor(df$GT , levels=c("0/0", "0/1", "1/1"))

my_comparisons <- list( c("0/0", "0/1"), c("0/0", "1/1"), c("0/1", "1/1") )

max_value<-max(merged_ranked_pheno[[pheno]],na.rm = TRUE)

p<-ggboxplot(df, x = "GT", y = paste(pheno),title = paste(""),conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text" )+ stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = max_value*1.5)     # Add global p-value
print(p)

##MAGI2

df <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/Publication/Mersad_LD_analysis/MAGI2.tsv", header=FALSE)
colnames(df)<-c("SAMPLE","ID","GT")
merged_ranked_pheno <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/Phenotypes_6218/merged_ranked_pheno.tsv")
merged_ranked_pheno[merged_ranked_pheno == -9] <- NA

df<-merge(df,merged_ranked_pheno,by.x = "SAMPLE",by.y = "SAMPLE_NAME")
df$GT <- factor(df$GT , levels=c("0/0", "0/1", "1/1"))

my_comparisons <- list( c("0/0", "0/1"), c("0/0", "1/1"), c("0/1", "1/1") )

max_value<-max(merged_ranked_pheno[[pheno]],na.rm = TRUE)

p<-ggboxplot(df, x = "GT", y = paste(pheno),title = paste(""),conf.int = TRUE, fill = "GT", width = 0.8, palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text" )+ stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = max_value*1.5)     # Add global p-value
print(p)



