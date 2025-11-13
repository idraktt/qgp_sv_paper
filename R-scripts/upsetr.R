library(UpSetR)
library(readr)
library(ggplot2)


upsetr_converted <- read.delim("/gpfs/ngsdata/QGP_projects/QF-QBB-RES-ACC-0032/QGP_SV/Publication/publication_plots/data/upsetr_converted.tsv", header=FALSE)

upsetr_plot<-as.data.frame(upsetr_converted)

myplot <- function(mydata,svtype) {
  plot<-ggplot(data = mydata,aes(x=factor(svtype))) + geom_bar(stat="count",width = 0.7,fill="steelblue") + theme_minimal()
}

colnames(upsetr_plot)[colnames(upsetr_plot)=="V1"] <-"ID"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V2"] <-"SVTYPE"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V3"] <-"Breakdancer"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V4"] <-"Breakseq2"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V5"] <-"CNVnator"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V6"] <-"Delly"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V7"] <-"ERDS"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V8"] <-"Genomestrip"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V9"] <-"Manta"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V10"] <-"Speedseq"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V11"] <-"Svaba"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V12"] <-"WHAM"
colnames(upsetr_plot)[colnames(upsetr_plot)=="V13"] <-"Unknown"

upsetr_plot = subset(upsetr_plot, select = -c(Unknown) )


upsetr_plot_deletions<-upsetr_plot[which(upsetr_plot$SVTYPE=='DEL'),]
upsetr_plot_duplications<-upsetr_plot[which(upsetr_plot$SVTYPE=='DUP'),]
upsetr_plot_inversions<-upsetr_plot[which(upsetr_plot$SVTYPE=='INV'),]

upsetr_plot_duplications_cnvnator_erds<-upsetr_plot_duplications %>% 
  filter(
    Speedseq == 1 & Delly == 1,
    Breakdancer == 0,
    Manta == 0,
    Genomestrip == 0,
    CNVnator == 0,
    ERDS == 0,
    Svaba == 0,
    WHAM == 0
  )

upsetr_plot_duplications = subset(upsetr_plot_duplications, select = -c(Breakdancer,Breakseq2,Genomestrip) )
upsetr_plot_inversions = subset(upsetr_plot_inversions, select = -c(Breakdancer,Breakseq2,Genomestrip,CNVnator,ERDS) )


upset(upsetr_plot_deletions,nsets = 10,show.numbers = FALSE,number.angles = 30,point.size = 1.5,line.size = 1,order.by = "freq",mainbar.y.label = "Caller Intersections",
    sets.bar.color = c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"),nintersects = 50,text.scale = 2,sets.x.label = "Deletions called Per caller")
upset(upsetr_plot_duplications,nsets = 7,show.numbers = FALSE,number.angles = 30,point.size = 1.5,line.size = 1,order.by = "freq",mainbar.y.label = "Caller Intersections",
      sets.bar.color = c("#E15759","#76B7B2","#59A14F","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"),nintersects = 30,text.scale = 2,sets.x.label = "Duplications called Per caller")
upset(upsetr_plot_inversions,nsets = 5,show.numbers = FALSE,number.angles = 20,point.size = 1.5,line.size = 1,order.by = "freq",mainbar.y.label = "Caller Intersections",
      sets.bar.color = c("#59A14F","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"),nintersects = 20,text.scale = 2,sets.x.label = "Inversions called Per caller")


c("#59A14F","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")

upset(upsetr_plot_deletions,nsets = 10,point.size = 2.5,line.size = 1,order.by = "degree",
      sets.bar.color = sample(colours(),10))

upset(upsetr_plot_duplications,nsets = 10,point.size = 2.5,line.size = 1,order.by = "freq",
      sets.bar.color = sample(colours(),10))
upset(upsetr_plot[which(upsetr_plot$SVTYPE=='DUP'),],nsets = 10,point.size = 2.5,line.size = 1,order.by = "degree",
      sets.bar.color = sample(colours(),10),nintersects = 50)

upset(upsetr_plot_inversions,nsets = 10,point.size = 2.5,line.size = 1,order.by = "freq",
      sets.bar.color = sample(colours(),10))

upset(upsetr_plot[which(upsetr_plot$SVTYPE=='INV'),],nsets = 10,point.size = 2.5,line.size = 1,order.by = "degree",
      sets.bar.color = sample(colours(),10),nintersects = 50)

