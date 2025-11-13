library("dplyr")
library("stringr")
library("janitor")
library("ggplot2")
library("patchwork")
library("uwot")

setwd("C:/Dropbox/SVAssociation/snpDataPCA/")

eigen_file=list.files(pattern = "\\.eigenvec$")
isgr_pheno = read.csv("../mergedQGP1KData/igsr-1000 genomes phase 3 release.tsv.tsv",header = TRUE, sep="\t" ) %>%
  dplyr::select(Population.code,Superpopulation.name)

integrated_tgp = read.csv("../mergedQGP1KData/integrated_call_samples_v3.20200731.ALL.ped",header = TRUE, sep="\t") %>% 
  dplyr::select(Individual.ID,Population) %>% mutate(STUDY="TGP") %>% rename("SAMPLE_NAME"="Individual.ID") %>%
  rename("pop"="Population") %>% left_join(isgr_pheno, c("pop"="Population.code")) %>% select(-pop) %>%
  rename("pop"="Superpopulation.name")

subpop_qgp = read.csv("../mergedQGP1KData/QGP_pop_mapping_final.tsv", header = TRUE, sep="\t") %>% 
  select(NAME,pop) %>% rename("SAMPLE_NAME"="NAME") %>% mutate(STUDY="QGP") 

subpop_qgp_name = read.csv("../mergedQGP1KData/QGP_pop_mapping_final.tsv", header = TRUE, sep="\t")
sv_qgp = read.csv("../mergedQGP1KData/QGP_hq.AF_SORTED_SORT_DEL_hwe.eigenvec", header = TRUE, sep = "\t")
subpop_qgp_name %>% inner_join(sv_qgp, by =c("SAMPLE_NAME"="IID"))

manifest_subpop= rbind(integrated_tgp,subpop_qgp) %>% as.data.frame()


col_vector <- unname(palette)

#col_vector = c("#86E05B","#C99887","#E06DB1","#78BCD6","#B54BDD","#DF7E4C","#D8D87B","#D0DECE","#76DEBB","#D2B1DE","#8277D1")
#col_vector=c("#88CADA","#D5DBBE","#7EE05E","#77DEB8","#7B8CD4","#A248DB","#A248DB","#DADE6E","#D4B5D6","#D9727A","#CB9C63","#DA72CA")
#col_vector=c("#A943DD","#E1D457","#79B6D3","#D6A9D5","#8379D2","#78DEC2","#D57F6C","#CCD396","#E06AC3","#D2D9D1","#81E163")
#col_vector=c("#788BD2","#A146DB","#DA698C","#89CADA","#CF8E64","#80E05D","#D5DABC","#77DDB7","#D776D6","#D4B5D6","#DADA6E")
#col_vector=c("#DE6FA7","#B251DB","#E0D959","#D68967","#D0D99E","#7FE063","#6A8D8D","#AF95DD","#DAD3D4","#87C6E7","#7DE1C5")

col_vector=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

QGP_eigenvec = read.csv("QGP_1KG_merged.eigenvec", header = FALSE, sep =" ")
QGP_eigenvec = QGP_eigenvec %>% inner_join(manifest_subpop, by =c("V1"="SAMPLE_NAME")) 
p1 = ggplot(QGP_eigenvec, aes(x=V3, y=V4, col=pop, shape=STUDY )) + geom_point(alpha = 10/10) +
  scale_color_manual(values=col_vector) +  theme_bw() +
  labs(x = "PC1",
       y = "PC2 ") +
  theme_bw(base_size = 12) 

ggsave(filename = "PCA_QGP_1KG.png", plot = p1, width = 7, height = 6, dpi = 300, units = "in")

pdf(paste( "PCA_QGP_1KG", ".pdf", sep=""), 
    height=6, width=7)
print(p1)
dev.off()

