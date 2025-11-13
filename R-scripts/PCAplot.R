library("dplyr")
library("stringr")
library("janitor")
library("ggplot2")
library("patchwork")
library("uwot")

setwd("C:/Dropbox/SVAssociation/mergedQGP1KData/")

eigen_file=list.files(pattern = "\\.eigenvec$")
isgr_pheno = read.csv("igsr-1000 genomes phase 3 release.tsv.tsv",header = TRUE, sep="\t" ) %>%
            dplyr::select(Population.code,Superpopulation.name)

integrated_tgp = read.csv("integrated_call_samples_v3.20200731.ALL.ped",header = TRUE, sep="\t") %>% 
              dplyr::select(Individual.ID,Population) %>% mutate(STUDY="TGP") %>% rename("SAMPLE_NAME"="Individual.ID") %>%
            rename("pop"="Population") %>% left_join(isgr_pheno, c("pop"="Population.code")) %>% select(-pop) %>%
            rename("pop"="Superpopulation.name")

subpop_qgp = read.csv("QGP_pop_mapping_final.tsv", header = TRUE, sep="\t") %>% 
                select(SAMPLE_NAME,pop) %>% mutate(STUDY="QGP") 

manifest_subpop= rbind(integrated_tgp,subpop_qgp) %>% as.data.frame()


col_vector <- unname(palette)



col_vector=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

eigen_file=list.files(pattern = "\\.eigenvec$")

for (filename in eigen_file) {
  finalName = sub('\\.eigenvec$', '', filename) 
  print (finalName)
QGP_eigenvec = read.csv(filename, header = TRUE, sep="\t")
QGP_eigenvec = QGP_eigenvec %>% inner_join(manifest_subpop, by =c("IID"="SAMPLE_NAME")) 
QGP_eigenvec_DUP =read.csv(filename, header = TRUE, sep="\t")

# palette <- distinctColorPalette(11)
# col_vector <- unname(palette)
QGP_eigenvec$PC1<-QGP_eigenvec$PC1*-1
p1 = ggplot(QGP_eigenvec, aes(x=PC1, y=PC2, col=pop, shape=STUDY )) + geom_point(alpha = 10/10) +
  scale_color_manual(values=col_vector) +  theme_bw()
ggsave(paste0("PCA","/",finalName,"_PC1-2",".png"), width=12, height=8, units="in", p1)

pdf(paste0("PCA","/",finalName,"_PC1-2",".pdf"), 
    height=6, width=7)
print(p1)
dev.off()

p2 = ggplot(QGP_eigenvec, aes(x=PC1, y=PC3, col=pop, shape=STUDY )) + geom_point(alpha = 10/10) +
  scale_color_manual(values=col_vector) +  theme_bw()
ggsave(paste0("PCA","/",finalName,"_PC1-3",".png"), width=12, height=8, units="in", p2)

pdf(paste0("PCA","/",finalName,"_PC1-3",".pdf"), 
    height=8, width=10)
print(p2)
dev.off()

p3 = ggplot(QGP_eigenvec, aes(x=PC2, y=PC3, col=pop, shape=STUDY )) + geom_point(alpha = 10/10) +
  scale_color_manual(values=col_vector) +  theme_bw()
ggsave(paste0("PCA","/",finalName,"_PC2-3",".png"), width=12, height=8, units="in", p3)

pdf(paste0("PCA","/",finalName,"_PC2-3",".pdf"), 
    height=8, width=10)
print(p3)
dev.off()

# QGP_eigenvec_DUP=  read.csv("mergedQGP1KGPShared_hwe.eigenvec", header = TRUE, sep="\t")
# QGP_eigenvec_DUP = QGP_eigenvec_DUP %>% inner_join(manifest_subpop, by =c("IID"="SAMPLE_NAME")) 
# p1 = ggplot(QGP_eigenvec_DUP, aes(x=PC1, y=PC2, col=pop, shape=STUDY )) + geom_point(alpha = 10/10) + 
# scale_color_manual(values=col_vector) + theme_bw()

QGP_eigenvec_DUP$X.FID = NULL
rownames(QGP_eigenvec_DUP) = QGP_eigenvec_DUP$IID
QGP_eigenvec_DUP$IID = NULL
mnist_umap= umap(QGP_eigenvec_DUP, pca = 10, n_neighbors = 16, min_dist = 0.001, init = "spca")
mnist_umap = mnist_umap %>% as.data.frame() 
mnist_umap= tibble::rownames_to_column(mnist_umap, "IID")

mnist_umap= mnist_umap %>% inner_join(manifest_subpop, by =c("IID"="SAMPLE_NAME"))
p4 = ggplot(mnist_umap, aes(x=V1, y=V2,col=pop, shape=STUDY)) + geom_point(alpha = .8) +
  scale_shape_manual(values=c( 16,4)) + 
  scale_size_manual(values=c(3,4)) +
  theme_bw() + scale_color_manual(values=col_vector)
   #scale_color_brewer(palette = "Paired") + theme_bw()
ggsave(paste0("UMAP","/",finalName,"_UMAP",".png"), width=12, height=8, units="in", p4)

if (file.exists(paste0(filename,".var"))){
  
loadings = read.table(paste0(filename,".var"), header = FALSE)
if (loadings[1,5] == "PC1"){
  loadings = loadings[-1,]
}
png(paste0("VAR","/",filename,".var.png"), res = 600, height = 8, width =14 , units = "in")
par(mfrow=c(5,1), mar = c( 1, 4, 1, 2 ))

for( i in 1:5 ) {
      print (i)
      plot( 1:nrow(loadings), abs(as.numeric(loadings[,i+4])), ylab=paste("PC ", i, sep="" ), ylim = c( 1, 6.5 ) )
}
dev.off()
}
}

