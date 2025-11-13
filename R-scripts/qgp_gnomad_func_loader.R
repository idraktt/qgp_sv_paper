#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Perform gene-level analyses for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & clean vcf2bed input
read.vcf2bed <- function(vcf2bed.in, aps.in=NULL){
  #Read data
  dat <- read.table(vcf2bed.in, comment.char="", header=T, sep="\t")
  colnames(dat)[1] <- "chrom"
  #Restrict to sites with at least one observed alternative allele
  dat <- dat[union(which(dat$AC>0), grep("MULTIALLELIC", dat$FILTER, fixed=T)), ]
  #Drop columns not being used (to save memory)
  cols.to.drop <- c("CHR2", "CPX_INTERVALS", 
                    "END", "SOURCE", "STRANDS", "UNRESOLVED_TYPE", 
                    "LINCRNA__LOF", "LINCRNA__DUP_LOF", "LINCRNA__COPY_GAIN", 
                    "LINCRNA__DUP_PARTIAL", "LINCRNA__MSV_EXON_OVR", 
                    "LINCRNA__INTRONIC", "LINCRNA__INV_SPAN", "LINCRNA__UTR")
  dat <- dat[, -which(colnames(dat) %in% cols.to.drop)]
  #Convert numeric columns
  numeric.columns <- sort(unique(c(grep("FREQ", colnames(dat), fixed=T), 
                                   grep("AN", colnames(dat), fixed=T), 
                                   grep("AC", colnames(dat), fixed=T), 
                                   grep("AF", colnames(dat), fixed=T))))
  numeric.columns <- setdiff(numeric.columns, grep("SPAN", colnames(dat), fixed=T))
  dat[, numeric.columns] <- apply(dat[, numeric.columns], 2, as.numeric)
  #Read & add APS, if optioned
  if(!is.null(aps.in)){
    aps <- read.table(aps.in, header=T, sep="\t", comment.char="")
    colnames(aps)[1] <- "VID"
    dat <- merge(dat, aps, by.x="name", by.y="VID", all.x=T, all.y=F, sort=F)
    dat$APS[which(dat$chrom %in% c("X", "Y"))] <- NA
  }
  return(dat)
}
#Process SNV gene data
read.snvdata <- function(SNVdata.in, gene.metadata.in){
  #Read & clean
  snv.data <- read.table(SNVdata.in, header=T, comment.char="")
  metadata <- read.table(gene.metadata.in, header=T)
  merged <- merge(x=snv.data, y=metadata, by="gene", sort=F)
  merged <- merged[which(!(merged$chrom %in% c("chrX", "chrY"))), ]
  #Assign oe deciles
  merged$mis_oe_dec <- ceiling(10*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_dec <- ceiling(10*rank(merged$ptv_oe)/(nrow(merged)+1))
  #Assign oe to 40 bins
  merged$mis_oe_binrank <- ceiling(40*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_binrank <- ceiling(40*rank(merged$ptv_oe)/(nrow(merged)+1))
  #Assign oe percentiles
  merged$mis_oe_cent <- ceiling(100*rank(merged$mis_oe)/(nrow(merged)+1))
  merged$ptv_oe_cent <- ceiling(100*rank(merged$ptv_oe)/(nrow(merged)+1))
  #Return formatted data
  return(merged)
}
#Read external gene list and restrict to autosomal genes considered in gnomAD
importGenelist <- function(genelist.dir, filename, gene.data){
  genes <- read.table(paste(genelist.dir, "/", filename, sep=""), header=F)
  genes <- as.character(genes[, 1])
  genes <- genes[which(genes %in% gene.data$gene)]
  return(genes)
}

#Count fraction of SV with functional effects by SVTYPE
count.fracLOFSV <- function(dat){
  subdat <- dat[-unique(c(grep("MULTIALLELIC", dat$FILTER, fixed=T), 
                          grep("UNRESOLVED", dat$FILTER, fixed=T))), ]
  subdat <- dat[which(!(dat$chrom %in% c("X", "Y"))), ]
  #All SV - pLoF
  all.lof <- length(which(!is.na(subdat$PROTEIN_CODING__LOF)))/nrow(subdat)
  #Deletions - pLoF
  del.lof <- length(which(subdat$SVTYPE=="DEL" & !is.na(subdat$PROTEIN_CODING__LOF)))/length(which(subdat$SVTYPE=="DEL"))
  #Insertions - pLoF
  ins.lof <- length(which(subdat$SVTYPE=="INS" & !is.na(subdat$PROTEIN_CODING__LOF)))/length(which(subdat$SVTYPE=="INS"))
  #Inversions - pLoF
  inv.lof <- length(which(subdat$SVTYPE=="INV" & !is.na(subdat$PROTEIN_CODING__LOF)))/length(which(subdat$SVTYPE=="INV"))
  #Complex - pLoF
  cpx.lof <- length(which(subdat$SVTYPE=="CPX" & !is.na(subdat$PROTEIN_CODING__LOF)))/length(which(subdat$SVTYPE=="CPX"))
  #DUP/CPX - CG
  dup.cg <- length(which(subdat$SVTYPE %in% c("DUP", "CPX") & !is.na(subdat$PROTEIN_CODING__COPY_GAIN)))/length(which(subdat$SVTYPE %in% c("DUP", "CPX")))
  #Duplications - pLoF
  dup.lof <- length(which(subdat$SVTYPE=="DUP" & !is.na(subdat$PROTEIN_CODING__DUP_LOF)))/length(which(subdat$SVTYPE=="DUP"))
  # #Duplications - partial gene
  # dup.partial <- length(which(subdat$SVTYPE=="DUP" & !is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)))/length(which(subdat$SVTYPE=="DUP"))
  #Inversions + complex - span
  invcpx.span <- length(which(subdat$SVTYPE %in% c("INV", "CPX") & !is.na(subdat$PROTEIN_CODING__INV_SPAN)))/length(which(subdat$SVTYPE %in% c("INV", "CPX")))
  #Intronic, no coding effect
  all.intronic <- length(which(is.na(subdat$PROTEIN_CODING__LOF) 
                               & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                               & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                               & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                               & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                               & !is.na(subdat$PROTEIN_CODING__INTRONIC)))/nrow(subdat)
  #Promoter, no coding effect
  all.promoter <- length(which(is.na(subdat$PROTEIN_CODING__LOF) 
                               & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                               & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                               & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                               & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                               & !is.na(subdat$PROTEIN_CODING__PROMOTER)))/nrow(subdat)
  #Format output
  v.out <- c("all.lof"=all.lof, "del.lof"=del.lof, "ins.lof"=ins.lof, "inv.lof"=inv.lof, "cpx.lof"=cpx.lof, 
             "dup.cg"=dup.cg, "dup.lof"=dup.lof, "invcpx.span"=invcpx.span, "all.intronic"=all.intronic, "all.promoter"=all.promoter)
  return(v.out)
}
#Barplots of SV fractions by functional category
barplot.fracByEffect <- function(frac.mat, base.colors, category.labels){
  #Prep plot area
  ymax <- 1.05*max(frac.mat, na.rm=T)
  par(mar=c(2.5, 3, 1, 0.5), bty="n")
  plot(x=c(0, nrow(frac.mat)), y=c(0, ymax), type="n", 
       xlab="", xaxt="n", ylab="", yaxt="n", yaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4], 
       col="white", border=NA, bty="n")
  if(length(axTicks(2))>4){
    y.at <- axTicks(2)[seq(1, length(axTicks(2)), 2)]
    y.at <- c(y.at, y.at[length(y.at)]+(y.at[2]-y.at[1]))
  }else{
    y.at <- axTicks(2)
  }
  axis(2, at=y.at, labels=NA, tck=-0.06)
  axis(2, at=y.at, tick=F, cex.axis=0.8, las=2, line=-0.6, 
       labels=paste(round(100*y.at, 0), "%", sep=""))
  mtext(2, line=1.9, text="Pct. of SVs")
  
  #Add bars
  sapply(1:nrow(frac.mat), function(i){
    rect(xleft=i-0.85, xright=i-0.75, ybottom=0, ytop=frac.mat[i, 1], col="black")
    rect(xleft=i-c(0.75, 0.55, 0.35), xright=i-c(0.55, 0.35, 0.15), 
         ybottom=0, ytop=frac.mat[i, 2:4], 
         col=c(adjustcolor(base.colors[i], alpha=0.15), 
               adjustcolor(base.colors[i], alpha=0.5), 
               base.colors[i]))
    axis(1, at=i-0.5, tick=F, labels=category.labels[i], cex.axis=0.8, padj=1, line=-1.5)
    text(x=i-0.8, y=frac.mat[i, 1]-(0.05*(par("usr")[4]-par("usr")[3])), pos=3, xpd=T, 
         labels=paste(round(100*frac.mat[i, 1], 1), "%", sep=""), cex=0.85)
  })
}

#Get proportion of singletons and 95% CI for a set of variants
calc.fracSingletons.singleClass <- function(dat, boot.n=100, conf=0.95, aps=F){
  dat <- dat[which(!(dat$chrom %in% c("X", "Y"))), ]
  ACs <- dat$AC
  v.aps <- dat$APS
  ACs <- ACs[which(!is.na(ACs) & ACs>0)]
  v.aps <- v.aps[which(!is.na(v.aps))]
  if(aps==T){
    helper.getFracSingle <- function(v.aps, indices){sum(v.aps[indices])/length(v.aps[indices])}
    point.est <- helper.getFracSingle(v.aps, indices=1:length(v.aps))
    calc.ci <- function(v.aps, n, conf){
      set.seed(0)
      boot.obj <- boot(data=v.aps, statistic=helper.getFracSingle, R=n)
      ci <- boot.ci(boot.obj, conf=conf, type="basic")$basic[4:5]
      return(ci)
    }
    ci <- calc.ci(v.aps, n=boot.n, conf=conf)
  }else{
    helper.getFracSingle <- function(ACs, indices){length(which(ACs[indices]==1))/length(ACs[indices])}
    point.est <- helper.getFracSingle(ACs, indices=1:length(ACs))
    calc.ci <- function(ACs, n, conf){
      set.seed(0)
      boot.obj <- boot(data=ACs, statistic=helper.getFracSingle, R=n)
      ci <- boot.ci(boot.obj, conf=conf, type="basic")$basic[4:5]
      return(ci)
    }
    ci <- calc.ci(ACs, n=boot.n, conf=conf)
  }
  return(c(point.est, ci))
}
#Count fraction of singletons per functional class
count.fracSingletons <- function(dat, aps=F){
  subdat <- dat[-unique(c(grep("MULTIALLELIC", dat$FILTER, fixed=T), 
                          grep("UNRESOLVED", dat$FILTER, fixed=T))), ]
  subdat <- dat[which(!(dat$chrom %in% c("X", "Y"))), ]
  lof.all.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF))
  lof.del.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF) & subdat$SVTYPE=="DEL")
  lof.invcpx.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF) & subdat$SVTYPE %in% c("INV", "CPX"))
  lof.ins.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF) & subdat$SVTYPE=="INS")
  # lof.dup.idx <- which(!is.na(subdat$PROTEIN_CODING__LOF) & subdat$SVTYPE=="DUP")
  cg.idx <- which(!is.na(subdat$PROTEIN_CODING__COPY_GAIN))
  dup.lof.idx <- which(!is.na(subdat$PROTEIN_CODING__DUP_LOF))
  inv.span.idx <- which(!is.na(subdat$PROTEIN_CODING__INV_SPAN))
  intronic.idx <- which(is.na(subdat$PROTEIN_CODING__LOF) 
                        & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                        & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                        & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                        & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                        & !is.na(subdat$PROTEIN_CODING__INTRONIC))
  promoter.idx <- which(is.na(subdat$PROTEIN_CODING__LOF) 
                        & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                        & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                        & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                        & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                        & is.na(subdat$PROTEIN_CODING__INTRONIC)
                        & !is.na(subdat$PROTEIN_CODING__PROMOTER))
  intergenic.idx <- which(is.na(subdat$PROTEIN_CODING__LOF) 
                          & is.na(subdat$PROTEIN_CODING__DUP_LOF)
                          & is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                          & is.na(subdat$PROTEIN_CODING__DUP_PARTIAL)
                          & is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR)
                          & is.na(subdat$PROTEIN_CODING__INTRONIC)
                          & is.na(subdat$PROTEIN_CODING__PROMOTER)
                          & subdat$PROTEIN_CODING__INTERGENIC=="True")
  all.sv.idx <- 1:nrow(subdat)
  fracSingles <- do.call("rbind", lapply(list(lof.all.idx, lof.del.idx, lof.ins.idx, lof.invcpx.idx, 
                                              cg.idx, dup.lof.idx, inv.span.idx, 
                                              intronic.idx, promoter.idx, intergenic.idx, 
                                              all.sv.idx), function(idx){
                                                c(length(idx), 
                                                  calc.fracSingletons.singleClass(dat=subdat[idx, ], aps=aps))
                                              }))
  if(aps==T){
    colnames(fracSingles) <- c("N_SV", "APS", "lower95CI", "upper95CI")
  }else{
    colnames(fracSingles) <- c("N_SV", "fracSingletons", "lower95CI", "upper95CI")
  }
  #Format output
  fracSingles <- cbind("Annotation"=c("LOF.ANY", "LOF.DEL", "LOF.INS", "LOF.INV_CPX", 
                                      "CG", "DUP_LOF", "INV_SPAN", 
                                      "INTRONIC", "PROMOTER", "INTERGENIC", "ALL.SV"), 
                       fracSingles)
  return(fracSingles)
}
#Count fraction of singletons per SVTYPE, split by those with and without coding effects
count.fracSingletons.coding_vs_noncoding <- function(dat, aps=F){
  subdat <- dat[-unique(c(grep("MULTIALLELIC", dat$FILTER, fixed=T), 
                          grep("UNRESOLVED", dat$FILTER, fixed=T))), ]
  subdat <- dat[which(!(dat$chrom %in% c("X", "Y"))), ]
  sv.with.coding.effects <- which(!is.na(subdat$PROTEIN_CODING__LOF)
                                  | !is.na(subdat$PROTEIN_CODING__DUP_LOF)
                                  | !is.na(subdat$PROTEIN_CODING__COPY_GAIN)
                                  | !is.na(subdat$PROTEIN_CODING__MSV_EXON_OVR))
  intergenic.sv <- setdiff(which(subdat$PROTEIN_CODING__INTERGENIC=="True"), 
                           sv.with.coding.effects)
  all.sv.idx <- 1:nrow(subdat)
  del.lof <- intersect(which(subdat$SVTYPE=="DEL"), sv.with.coding.effects)
  del.noncoding <- intersect(which(subdat$SVTYPE=="DEL"), intergenic.sv)
  dup.ied <- intersect(which(subdat$SVTYPE=="DUP" & !is.na(subdat$PROTEIN_CODING__DUP_LOF)), 
                       sv.with.coding.effects)
  dup.cg <- intersect(which(subdat$SVTYPE=="DUP" & !is.na(subdat$PROTEIN_CODING__COPY_GAIN)), 
                      sv.with.coding.effects)
  dup.noncoding <- intersect(which(subdat$SVTYPE=="DUP"), intergenic.sv)
  ins.lof <- intersect(which(subdat$SVTYPE=="INS"), sv.with.coding.effects)
  ins.noncoding <- intersect(which(subdat$SVTYPE=="INS"), intergenic.sv)
  inv.lof <- intersect(which(subdat$SVTYPE=="INV"), sv.with.coding.effects)
  inv.noncoding <- intersect(which(subdat$SVTYPE=="INV"), intergenic.sv)
  cpx.lof <- intersect(which(subdat$SVTYPE=="CPX" & !is.na(subdat$PROTEIN_CODING__LOF)), 
                       sv.with.coding.effects)
  cpx.cg <- intersect(which(subdat$SVTYPE=="CPX" & !is.na(subdat$PROTEIN_CODING__COPY_GAIN)), 
                      sv.with.coding.effects)
  cpx.noncoding <- intersect(which(subdat$SVTYPE=="CPX"), intergenic.sv)
  fracSingles <- do.call("rbind", lapply(list(all.sv.idx, del.lof, del.noncoding, 
                                              dup.ied, dup.cg, dup.noncoding, 
                                              ins.lof, ins.noncoding, inv.lof, inv.noncoding, 
                                              cpx.lof, cpx.cg, cpx.noncoding), function(idx){
                                                c(length(idx), 
                                                  calc.fracSingletons.singleClass(dat=subdat[idx, ], aps=aps))
                                              }))
  if(aps==T){
    colnames(fracSingles) <- c("N_SV", "APS", "lower95CI", "upper95CI") 
  }else{
    colnames(fracSingles) <- c("N_SV", "fracSingletons", "lower95CI", "upper95CI")
  }
  #Format output
  fracSingles <- cbind("Annotation"=c("ALL.SV", "DEL.LOF", "DEL.NON", 
                                      "DUP.IED", "DUP.CG", "DUP.NON", 
                                      "INS.LOF", "INS.NON", "INV.LOF", "INV.NON", 
                                      "CPX.LOF", "CPX.CG", "CPX.NON"), 
                       fracSingles)
  return(fracSingles)
}
#Dotplot of fraction of singletons per functional annotation by class
dotplot.fracSingletons <- function(fracSingletons, horiz=F, add.pt.labs=F, add.N=T, aps=F){
  plot.dat <- as.data.frame(apply(fracSingletons[, -1], 2, as.numeric))
  colnames(plot.dat)[2] <- "fracSingletons"
  rownames(plot.dat) <- fracSingletons[, 1]
  class.colors <- c(rep(svtypes$color[which(svtypes$svtype=="DEL")], 4), 
                    svtypes$color[which(svtypes$svtype=="DUP")], 
                    svtypes$color[which(svtypes$svtype=="MCNV")], 
                    svtypes$color[which(svtypes$svtype=="INV")], 
                    rep("grey35", 2), "grey70", "black")
  class.labels <- c("pLoF (All)", "pLoF (DEL)", "pLoF (INS)", "pLoF (INV & CPX)", 
                    "CG", "IED", "Whole-Gene INV", 
                    "Intronic", "Promoter", "Intergenic", "All Autosomal SVs")
  new.order <- rev(c(which(rownames(plot.dat)=="ALL.SV"), 
                     which(rownames(plot.dat)=="LOF.ANY"), 
                     intersect(order(-plot.dat$fracSingletons), 
                               which(!(rownames(plot.dat) %in% c("ALL.SV", "LOF.ANY"))))))
  if(horiz==T){
    new.order <- rev(new.order)
  }
  plot.dat <- plot.dat[new.order, ]
  class.colors <- class.colors[new.order]
  class.labels <- class.labels[new.order]
  if(horiz==T){
    par(mar=c(3.5, 2.5, 0.25, 0.25), bty="n")
    plot(y=range(plot.dat[, -1]), 
         x=c(0.4, nrow(plot.dat)-0.1), 
         type="n", xaxt="n", xlab="", yaxt="n", ylab="")
    if(aps==T){
      abline(h=0, lty=2, col="gray70")
    }else{
      abline(h=plot.dat$fracSingletons[which(rownames(plot.dat)=="ALL.SV")], lty=2)
    }
    segments(y0=plot.dat[, 3], y1=plot.dat[, 4], 
             x0=(1:nrow(plot.dat))-0.5, 
             x1=(1:nrow(plot.dat))-0.5, 
             lwd=2, lend="round", 
             col=class.colors)
    points(y=plot.dat[, 2], x=c(1:nrow(plot.dat))-0.5, pch=19, 
           col=class.colors)
    if(add.pt.labs==T){
      text(x=c(2:nrow(plot.dat))-0.65, y=plot.dat[-1, 2], 
           col=class.colors[-1], cex=0.5, pos=4, xpd=T, 
           labels=paste(format(round(100*plot.dat[-1, 2], 1), nsmall=1), "%", sep=""))
      text(x=1-0.65, y=plot.dat[1, 2]+0.01, 
           col=class.colors[1], cex=0.5, pos=4, xpd=T, 
           labels=paste(format(round(100*plot.dat[1, 2], 1), nsmall=1), "%", sep=""))
    }
    par(xpd=T)
    sapply(1:nrow(plot.dat), function(i){
      if(add.N==T){
        text(x=i-0.15, y=par("usr")[3]-(0.03*(par("usr")[4]-par("usr")[3])), 
             labels=class.labels[i], srt=40, pos=2, cex=0.7)
        text(x=i+0.175, y=par("usr")[3]-(0.05*(par("usr")[4]-par("usr")[3])), srt=40, col="gray40", cex=0.5, pos=2, 
             labels=paste("N=", prettyNum(plot.dat[i, 1], big.mark=","), " SVs", sep=""))
      }else{
        text(x=i-0.08, y=par("usr")[3]-(0.03*(par("usr")[4]-par("usr")[3])), 
             labels=class.labels[i], srt=40, pos=2, cex=0.7)
      }
    })
    par(xpd=F)
    if(aps==T){
      y.at <- seq(-1, 1, 0.05)
      axis(2, at=y.at, labels=NA, tck=-0.035)
      sapply(y.at, function(x){
        axis(2, at=x, line=-0.7, cex.axis=0.65, tick=F, las=2)
      })
      mtext(2, text="APS", line=1.5, cex=0.8)
    }else{
      axis(2, at=seq(0.4, 0.7, 0.1), labels=NA, tck=-0.035)
      sapply(seq(0.4, 0.7, 0.1), function(x){
        axis(2, at=x, line=-0.7, cex.axis=0.65, tick=F, las=2, 
             labels=paste(round(x*100, 0), "%", sep=""))
      })
      mtext(2, text="Singleton Proportion", line=1.5, cex=0.8)
    }
  }else{
    par(mar=c(0.5, 6, 2.25, 0.5), bty="n")
    plot(x=range(plot.dat[, -1]), 
         y=c(0.25, nrow(plot.dat)-0.25), 
         type="n", xaxt="n", xlab="", yaxt="n", ylab="")
    if(aps==T){
      abline(v=0, lty=2, col="gray70")
    }else{
      abline(v=plot.dat$fracSingletons[which(rownames(plot.dat)=="ALL.SV")], lty=2)
    }
    segments(x0=plot.dat[, 3], x1=plot.dat[, 4], 
             y0=(1:nrow(plot.dat))-0.5, 
             y1=(1:nrow(plot.dat))-0.5, 
             lwd=2, lend="round", 
             col=class.colors)
    points(x=plot.dat[, 2], y=c(1:nrow(plot.dat))-0.5, pch=19, 
           col=class.colors)
    axis(2, at=(1:nrow(plot.dat))-0.3, tick=F, line=-0.8, las=2, cex.axis=0.8, 
         labels=class.labels)
    sapply(1:nrow(plot.dat), function(i){
      axis(2, at=i-0.7, tick=F, line=-0.8, las=2, cex.axis=0.65, 
           labels=paste("N=", prettyNum(plot.dat[i, 1], big.mark=","), " SVs", sep=""), 
           col.axis=class.colors[i])
    })
    if(aps==T){
      x.at <- seq(-1, 1, 0.05)
      axis(3, at=x.at, labels=NA, tck=-0.035)
      sapply(x.at, function(x){
        axis(3, at=x, line=-0.7, cex.axis=0.65, tick=F)
      })
      mtext(3, text="APS", line=1.15, cex=0.9)
    }else{
      axis(3, at=seq(0.4, 0.7, 0.1), labels=NA, tck=-0.035)
      sapply(seq(0.4, 0.7, 0.1), function(x){
        axis(3, at=x, line=-0.7, cex.axis=0.65, tick=F, 
             labels=paste(round(x*100, 0), "%", sep=""))
      })
      mtext(3, text="Singleton Proportion", line=1.15, cex=0.9)
    }
  }
}
#Dotplot of fraction of singletons per functional annotation by class
dotplot.fracSingletons.coding_vs_non <- function(fracSingletons, add.N=T, aps=F, ylims=NULL){
  rownames <- fracSingletons[, 1]
  plot.dat <- as.data.frame(apply(fracSingletons[, -1], 2, as.numeric))
  colnames(plot.dat)[2] <- "fracSingletons"
  rownames(plot.dat) <- fracSingletons[, 1]
  class.colors <- c("black", 
                    svtypes$color[which(svtypes$svtype=="DEL")], 
                    adjustcolor(svtypes$color[which(svtypes$svtype=="DEL")], alpha=0.3), 
                    rep(svtypes$color[which(svtypes$svtype=="DUP")], 2), 
                    adjustcolor(svtypes$color[which(svtypes$svtype=="DUP")], alpha=0.3), 
                    svtypes$color[which(svtypes$svtype=="INS")], 
                    adjustcolor(svtypes$color[which(svtypes$svtype=="INS")], alpha=0.3), 
                    svtypes$color[which(svtypes$svtype=="INV")], 
                    adjustcolor(svtypes$color[which(svtypes$svtype=="INV")], alpha=0.3), 
                    rep(svtypes$color[which(svtypes$svtype=="CPX")], 2), 
                    adjustcolor(svtypes$color[which(svtypes$svtype=="CPX")], alpha=0.3))
  label.colors <- c("black", 
                    rep(svtypes$color[which(svtypes$svtype=="DEL")], 2), 
                    rep(svtypes$color[which(svtypes$svtype=="DUP")], 3), 
                    rep(svtypes$color[which(svtypes$svtype=="INS")], 2), 
                    rep(svtypes$color[which(svtypes$svtype=="INV")], 2), 
                    rep(svtypes$color[which(svtypes$svtype=="CPX")], 3))
  class.labels <- c("All Autosomal SVs", "DEL (pLoF)", "DEL (Intergenic)", 
                    "DUP (IED)", "DUP (CG)", "DUP (Intergenic)", 
                    "INS (pLoF)", "INS (Intergenic)", "INV (pLoF)", "INV (Intergenic)", 
                    "CPX (pLoF)", "CPX (CG)", "CPX (Intergenic)")
  par(mar=c(3.5, 2.25, 0.25, 0.25), bty="n")
  if(is.null(ylims)){
    ylims <- range(plot.dat[, -1])
  }
  plot(y=ylims, x=c(0.25, nrow(plot.dat)-0.4), 
       type="n", xaxt="n", xlab="", yaxt="n", ylab="")
  if(aps==T){
    abline(h=0, lty=2, col="gray70", lwd=1.5)
  }else{
    abline(h=plot.dat$fracSingletons[which(rownames(plot.dat)=="ALL.SV")], lty=2, lwd=1.5)
  }
  segments(y0=plot.dat[, 3], y1=plot.dat[, 4], 
           x0=(1:nrow(plot.dat))-0.5, 
           x1=(1:nrow(plot.dat))-0.5, 
           lwd=2, lend="round", 
           col="white")
  segments(y0=plot.dat[, 3], y1=plot.dat[, 4], 
           x0=(1:nrow(plot.dat))-0.5, 
           x1=(1:nrow(plot.dat))-0.5, 
           lwd=2, lend="round", 
           col=class.colors)
  points(y=plot.dat[, 2], x=c(1:nrow(plot.dat))-0.5, pch=19, 
         col="white")
  points(y=plot.dat[, 2], x=c(1:nrow(plot.dat))-0.5, pch=19, 
         col=class.colors)
  # text(x=c(1:nrow(plot.dat))-0.625, y=plot.dat[, 2], 
  #      col=label.colors, cex=0.5, pos=4, xpd=T, 
  #      labels=paste(format(round(100*plot.dat[, 2], 1), nsmall=1), "%", sep=""))
  par(xpd=T)
  sapply(1:nrow(plot.dat), function(i){
    if(add.N==T){
      text(x=i-0.15, y=par("usr")[3]-(0.03*(par("usr")[4]-par("usr")[3])), 
           labels=class.labels[i], srt=40, pos=2, cex=0.7)
      text(x=i+0.175, y=par("usr")[3]-(0.05*(par("usr")[4]-par("usr")[3])), srt=40, col="gray40", cex=0.5, pos=2, 
           labels=paste("N=", prettyNum(plot.dat[i, 1], big.mark=","), " SV", sep=""))
    }else{
      text(x=i-0.08, y=par("usr")[3]-(0.03*(par("usr")[4]-par("usr")[3])), 
           labels=class.labels[i], srt=40, pos=2, cex=0.7)
    }
  })
  par(xpd=F)
  if(aps==T){
    y.at <- seq(-1, 1, 0.05)
    axis(2, at=y.at, labels=NA, tck=-0.02)
    sapply(y.at, function(x){
      axis(2, at=x, line=-0.7, cex.axis=0.65, tick=F, las=2)
    })
    mtext(2, text="APS", line=1.5, cex=0.8)
  }else{
    axis(2, at=seq(0.4, 1, 0.1), labels=NA, tck=-0.04)
    sapply(seq(0.4, 1, 0.1), function(x){
      axis(2, at=x, line=-0.7, cex.axis=0.65, tick=F, las=2, 
           labels=paste(round(x*100, 0), "%", sep=""))
    })
    mtext(2, text="Singleton Proportion", line=1.5, cex=0.8)
  }
}


#Organize table of genes from SV data
getSVdat <- function(dat, genes, prefix=NULL){
  #pLoF
  lof.any <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF))], split=",")))
  lof.del <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF)
                                                                        & dat$SVTYPE=="DEL")], split=",")))
  lof.other <- as.character(unlist(strsplit(dat$PROTEIN_CODING__LOF[which(!is.na(dat$PROTEIN_CODING__LOF)
                                                                          & dat$SVTYPE!="DEL")], split=",")))
  #Dup GC
  cg.dup <- as.character(unlist(strsplit(dat$PROTEIN_CODING__COPY_GAIN[which(!is.na(dat$PROTEIN_CODING__COPY_GAIN))], split=",")))
  #Dup pLoF
  plof.dup <- as.character(unlist(strsplit(dat$PROTEIN_CODING__DUP_LOF[which(!is.na(dat$PROTEIN_CODING__DUP_LOF))], split=",")))
  #Inversion span
  inv.span <- as.character(unlist(strsplit(dat$PROTEIN_CODING__INV_SPAN[which(!is.na(dat$PROTEIN_CODING__INV_SPAN))], split=",")))
  #Collect vector per gene
  res <- as.data.frame(t(sapply(genes, function(gene){
    g.lof.any <- length(which(lof.any==gene))
    g.lof.del <- length(which(lof.del==gene))
    g.lof.other <- length(which(lof.other==gene))
    g.cg.dup <- length(which(cg.dup==gene))
    g.plof.dup <- length(which(plof.dup==gene))
    g.inv.span <- length(which(inv.span==gene))
    g.out <- as.integer(c(g.lof.any, g.lof.del, g.lof.other, g.cg.dup, g.plof.dup, g.inv.span))
    g.out[which(is.na(g.out))] <- 0
    g.out <- c(gene, g.out)
    return(g.out)
  })))
  colnames(res) <- c("gene", "lof.any", "lof.del", "lof.other", "cg", "plof", "inv")
  if(!is.null(prefix)){
    colnames(res)[-1] <- paste(prefix, colnames(res)[-1], sep=".")
  }
  rownames(res) <- 1:nrow(res)
  return(res)
}
#Get merged table of SNV stats and SV stats by frequency bin
getSVdat.all <- function(dat, snv.data, include.cpx=T, require.SR=F){
  #Gather data
  genes <- sort(unique(as.character(snv.data$gene)))
  if(include.cpx==F){
    dat <- dat[which(dat$SVTYPE != "CPX"), ]
  }
  if(require.SR==T){
    dat <- dat[grep("SR", dat$EVIDENCE), ]
  }
  all.counts <- getSVdat(dat, genes, prefix="all")
  common.counts <- getSVdat(dat[which(dat$AF>=0.01), ], genes, prefix="common")
  rare.counts <- getSVdat(dat[which(dat$AF<0.01 & dat$AC>1), ], genes, prefix="rare")
  singleton.counts <- getSVdat(dat[which(dat$AC==1), ], genes, prefix="singleton")
  #Merge data
  merged <- merge(x=snv.data, y=all.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=common.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=rare.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=singleton.counts, by="gene", all.x=T, sort=F)
  merged[, -c(which(colnames(merged) %in% c("gene", "chrom")))] <- apply(merged[, -c(which(colnames(merged) %in% c("gene", "chrom")))], 2, as.numeric)
  return(merged)
}

#Stacked barplot of % of genes with at least one SV per category
barplot.genesFracByEffect <- function(gene.data){
  ngenes <- nrow(gene.data)
  #Get cumulative fractions per type
  plot.dat <- sapply(c("lof.any", "lof.del", "lof.other", "plof", "cg", "inv"), function(effect){
    idxs <- grep(paste(".", effect, sep=""), colnames(gene.data), fixed=T)
    sing.hit <- length(which(gene.data[, idxs[4]]>0))/ngenes
    singOrRare.hit <- length(which(gene.data[, idxs[4]]>0 | gene.data[, idxs[3]]>0))/ngenes
    singOrRareOrCommon.hit <- length(which(gene.data[, idxs[4]]>0 | gene.data[, idxs[3]]>0 | gene.data[, idxs[2]]>0))/ngenes
    return(c(sing.hit, singOrRare.hit, singOrRareOrCommon.hit))
  })
  
  #Prepare plot area
  par(bty="n", mar=c(0.1, 5, 2.5, 0.75))
  plot(x=c(0, 1), y=c(0, -6), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="")
  axis(2, at=-(0:5)-0.5, tick=F, las=2, cex.axis=0.8, line=-1.2, 
       labels=c("pLoF (All)", "pLoF (DEL Only)", "pLoF (Not DEL)", 
                "IED", "CG", "Inverted Gene"))
  axis(3, at=seq(0, 1, 0.2), labels=NA, tck=-0.04)
  sapply(seq(0, 1, 0.2), function(l){
    axis(3, at=l, tick=F, line=-0.6, cex.axis=0.8, 
         labels=paste(100*l, "%", sep=""))
  })
  mtext(3, line=1.4, text="Pct. of All Autosomal Genes")
  
  #Plot bars
  bar.base.cols <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                     svtypes$color[which(svtypes$svtype=="DEL")], 
                     svtypes$color[which(svtypes$svtype=="DEL")], 
                     svtypes$color[which(svtypes$svtype=="MCNV")], 
                     svtypes$color[which(svtypes$svtype=="DUP")], 
                     svtypes$color[which(svtypes$svtype=="INV")])
  sapply(1:ncol(plot.dat), function(i){
    rect(xleft=0, xright=1, ybottom=-i+0.8, ytop=-i+0.2, 
         col="white", border="gray85")
    rect(xleft=c(0, plot.dat[, i]), xright=c(plot.dat[, i], 1), 
         ybottom=-i+0.2, ytop=-i+0.8, border=NA, bty="n", 
         col=c(bar.base.cols[i], 
               adjustcolor(bar.base.cols[i], alpha=0.5), 
               adjustcolor(bar.base.cols[i], alpha=0.15), 
               "gray98"))
    rect(xleft=0, xright=max(plot.dat[, i]), 
         ybottom=-i+0.2, ytop=-i+0.8, col=NA)
    text(x=max(plot.dat[, i]), y=-i+0.45, cex=0.7, pos=4, 
         labels=paste(format(round(100*max(plot.dat[, i]), 1), nsmall=1), "%", sep=""))
  })
}

#Mini histogram of # of SV per gene by functional effect
minihist.byFunc <- function(gene.data, func, title=NULL, base.col, x.max=5){
  #Get hist data
  counts <- sapply(0:x.max, function(k){
    length(which(gene.data[, which(colnames(gene.data) == paste("all", func, sep="."))]==k))
  })
  counts <- c(counts, length(which(gene.data[, which(colnames(gene.data) == paste("all", func, sep="."))]>x.max)))
  counts <- counts/1000
  #Plot hist
  par(mar=c(2, 3, 1.5, 0.5), bty="n")
  plot(x=c(0, x.max+2), y=c(0, max(counts)), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=(0:(x.max+1))+0.1, xright=c(0:(x.max+1))+0.9, 
       ybottom=0, ytop=counts, col=base.col)
  sapply((0:x.max), function(k){
    axis(1, at=k+0.5, labels=k, cex.axis=0.8, tick=F, line=-1.3)
  })
  axis(1, at=x.max+1.5, labels=paste(">", x.max, sep=""), cex.axis=0.8, tick=F, line=-1.3)
  
  mtext(1, text="SV per Gene", line=0.75, cex=0.85)
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, labels=paste(axTicks(2), "k", sep=""), 
       line=-0.4, las=2, cex.axis=0.8)
  mtext(2, text="Genes", line=1.75, cex=0.9)
  mtext(3, text=title)
}


##############################
###COUNT PER GENE CORRELATIONS
##############################
#Run permutation test to gather expected correlation between two functional classes
permute.perGeneCounts <- function(gene.data, x.category, y.category, times=1000, max.ax=3){
  #Get actual values
  x <- gene.data[, which(colnames(gene.data)==x.category)]
  y <- gene.data[, which(colnames(gene.data)==y.category)]
  #Get shuffled values
  shuf.res <- lapply(1:times, function(i){
    x.shuf <- sample(x, size=length(x), replace=F)
    y.shuf <- sample(y, size=length(y), replace=F)
    shuf.dat <- sapply(0:(max.ax+1), function(x.ct){
      sapply(0:(max.ax+1), function(y.ct){
        if(x.ct>max.ax & y.ct>max.ax){
          length(which(x.shuf>max.ax &
                         y.shuf>max.ax))
        }else if(x.ct>max.ax & y.ct<=max.ax){
          length(which(x.shuf>max.ax &
                         y.shuf==y.ct))
        }else if(x.ct<=max.ax & y.ct>max.ax){
          length(which(x.shuf==x.ct &
                         y.shuf>max.ax))
        }else{
          length(which(x.shuf==x.ct &
                         y.shuf==y.ct))
        }
      })
    })
  })
  #Get means of shuffled values
  shuf.dat.means <- sapply(0:(max.ax+1), function(x.ct){
    sapply(0:(max.ax+1), function(y.ct){
      mean(sapply(shuf.res, function(df){
        return(df[x.ct+1, y.ct+1])
      }))
    })
  })
  return(shuf.dat.means)
}
#Generic heatmap of counts per gene for two functional classes
plot.perGeneHeatmaps <- function(gene.data, x.category, y.category, max.ax=3, max.col=0.2, 
                                 x.label=NULL, y.label=NULL){
  #Get plot data
  x.idx <- which(colnames(gene.data)==x.category)
  y.idx <- which(colnames(gene.data)==y.category)
  plot.dat <- sapply(0:(max.ax+1), function(x.ct){
    sapply(0:(max.ax+1), function(y.ct){
      if(x.ct>max.ax & y.ct>max.ax){
        length(which(gene.data[, x.idx]>max.ax &
                       gene.data[, y.idx]>max.ax))
      }else if(x.ct>max.ax & y.ct<=max.ax){
        length(which(gene.data[, x.idx]>max.ax &
                       gene.data[, y.idx]==y.ct))
      }else if(x.ct<=max.ax & y.ct>max.ax){
        length(which(gene.data[, x.idx]==x.ct &
                       gene.data[, y.idx]>max.ax))
      }else{
        length(which(gene.data[, x.idx]==x.ct &
                       gene.data[, y.idx]==y.ct))
      }
    })
  })
  plot.dat <- plot.dat/sum(plot.dat)
  plot.dat[which(plot.dat>max.col)] <- max.col
  plot.dat <- floor(1000*plot.dat)
  #Prep plot area
  par(mar=c(3, 3, 0.5, 0.5))
  plot(x=c(0, max.ax+2), y=c(0, max.ax+2), type="n", 
       xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")
  mtext(1, line=1.5, text=x.label, cex=0.9)
  axis(1, at=seq(0, max.ax, 1)+0.5, tick=F, line=-0.8, labels=seq(0, max.ax, 1))
  axis(1, at=max.ax+1.5, tick=F, line=-0.8, 
       labels=paste(">", max.ax, sep=""))
  mtext(2, line=1.5, text=y.label, cex=0.9)
  axis(2, at=seq(0, max.ax, 1)+0.5, tick=F, line=-0.8, labels=seq(0, max.ax, 1), las=2)
  axis(2, at=max.ax+1.5, tick=F, line=-0.8, las=2, 
       labels=paste(">", max.ax, sep=""))
  #Add heatmap
  col.pal <- colorRampPalette(c("#440154", "#365C8C", "#25A584", "#FDE725"))(length(seq(0, ceiling(1000*max.col), 1)))
  sapply(0:(max.ax+1), function(x){
    sapply(0:(max.ax+1), function(y){
      rect(xleft=x, xright=x+1, 
           ybottom=y, ytop=y+1, 
           border=NA, col=col.pal[plot.dat[y+1, x+1]+1])
    })
  })
}
#Key for per gene heatmaps
plot.perGeneHeatmapsKey <- function(max.col=0.2){
  col.pal <- colorRampPalette(c("#440154", "#365C8C", "#25A584", "#FDE725"))(length(seq(0, ceiling(1000*max.col), 1)))
  par(mar=rep(0.1, 4))
  plot(x=c(0, 1), y=c(0, length(col.pal)), type="n", 
       xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="")
  rect(xleft=0, xright=1, ybottom=(1:length(col.pal))-1, ytop=1:length(col.pal), 
       border=col.pal, col=col.pal)
  box()
}
#Scaled dotplots of counts per gene for two functional classes vs permuted expectation
plot.perGeneDotplots <- function(gene.data, x.category, y.category, max.ax=3, max.size=5000, 
                                 x.label=NULL, y.label=NULL, perm.times=1000){
  #Get empirical data
  x.idx <- which(colnames(gene.data)==x.category)
  y.idx <- which(colnames(gene.data)==y.category)
  plot.dat <- sapply(0:(max.ax+1), function(x.ct){
    sapply(0:(max.ax+1), function(y.ct){
      if(x.ct>max.ax & y.ct>max.ax){
        length(which(gene.data[, x.idx]>max.ax &
                       gene.data[, y.idx]>max.ax))
      }else if(x.ct>max.ax & y.ct<=max.ax){
        length(which(gene.data[, x.idx]>max.ax &
                       gene.data[, y.idx]==y.ct))
      }else if(x.ct<=max.ax & y.ct>max.ax){
        length(which(gene.data[, x.idx]==x.ct &
                       gene.data[, y.idx]>max.ax))
      }else{
        length(which(gene.data[, x.idx]==x.ct &
                       gene.data[, y.idx]==y.ct))
      }
    })
  })
  plot.dat[which(plot.dat>max.size)] <- max.size
  plot.dat <- sqrt(plot.dat)/(2*sqrt(max.size))
  #Get permuted expectations
  plot.exp <- permute.perGeneCounts(gene.data, x.category, y.category, times=perm.times, max.ax)
  plot.exp[which(plot.exp>max.size)] <- max.size
  plot.exp <- sqrt(plot.exp)/(2*sqrt(max.size))
  #Prep plot area
  par(mar=c(3, 3, 0.5, 0.5), bty="n")
  plot(x=c(0, max.ax+2), y=c(0, max.ax+2), type="n", 
       xlab="", ylab="", xaxt="n", yaxt="n")
  abline(h=par("usr")[1], v=par("usr")[3])
  mtext(1, line=1.5, text=x.label, cex=0.9)
  axis(1, at=seq(0, max.ax, 1)+0.5, tick=F, line=-0.8, labels=seq(0, max.ax, 1))
  axis(1, at=max.ax+1.5, tick=F, line=-0.8, 
       labels=paste(">", max.ax, sep=""))
  mtext(2, line=1.5, text=y.label, cex=0.9)
  axis(2, at=seq(0, max.ax, 1)+0.5, tick=F, line=-0.8, labels=seq(0, max.ax, 1), las=2)
  axis(2, at=max.ax+1.5, tick=F, line=-0.8, las=2, 
       labels=paste(">", max.ax, sep=""))
  #Add circles
  sapply(0:(max.ax+1), function(x){
    sapply(0:(max.ax+1), function(y){
      #Grey shaded circle for expectation
      draw.circle(x=x+0.5, y=y+0.5, radius=plot.exp[x+1, y+1], 
                  border=NA, col="gray80")
      #Black open circle for empirical observation
      draw.circle(x=x+0.5, y=y+0.5, radius=plot.dat[y+1, x+1], 
                  col=NA, lwd=2)
    })
  })
}
#Scaled dotplot key
plot.perGeneDotplotKey <- function(max.size=5000, marks=c(10, 100, 1000, 5000)){
  marks <- sqrt(marks)/(2*sqrt(max.size))
  par(mar=rep(0.1, 4), bty="n")
  plot(x=c(0, 1), y=c(0, length(marks)), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  sapply(1:length(marks), function(i){
    draw.circle(x=0.5, y=(i/1.5)-0.5, 
                radius=marks[i], 
                lwd=2, col="gray80")
  })
}


##################################
###CODE FOR CONSTRAINT COMPARISONS
##################################
#Prep covariates matrics for regression analysis
prepCovariates <- function(gene.data, y){
  cov <- data.frame("gene"=gene.data$gene, 
                    "y.SV_count"=gene.data[which(colnames(gene.data)==y)], 
                    "gene_length"=scale(log10(as.numeric(gene.data$gene_length)), center=T, scale=T), 
                    "exon_count"=scale(log10(as.numeric(gene.data$exon_count)), center=T, scale=T), 
                    # "exon_min"=scale(log10(as.numeric(gene.data$exon_min)), center=T, scale=T), 
                    # "exon_max"=scale(log10(as.numeric(gene.data$exon_min)), center=T, scale=T), 
                    # "exon_mean"=scale(log10(as.numeric(gene.data$exon_mean)), center=T, scale=T), 
                    "exon_median"=scale(log10(as.numeric(gene.data$exon_median)), center=T, scale=T), 
                    # "exon_harm.mean"=scale(log10(as.numeric(gene.data$exon_harm.mean)), center=T, scale=T), 
                    "exon_sum"=scale(log10(as.numeric(gene.data$exon_mean)), center=T, scale=T), 
                    "intron_count"=scale(log10(as.numeric(gene.data$intron_count)), center=T, scale=T), 
                    # "intron_min"=scale(log10(as.numeric(gene.data$intron_min)), center=T, scale=T), 
                    # "intron_max"=scale(log10(as.numeric(gene.data$intron_min)), center=T, scale=T), 
                    # "intron_mean"=scale(log10(as.numeric(gene.data$intron_mean)), center=T, scale=T), 
                    "intron_median"=scale(log10(as.numeric(gene.data$intron_median)), center=T, scale=T), 
                    # "intron_harm.mean"=scale(log10(as.numeric(gene.data$intron_harm.mean)), center=T, scale=T), 
                    "intron_sum"=scale(log10(as.numeric(gene.data$intron_mean)), center=T, scale=T), 
                    "segdup"=round(gene.data$segdup, 0))
  colnames(cov)[2] <- "y.SV_count"
  chrom.dummy.mat <- as.data.frame(sapply(unique(gene.data$chrom), function(chr){
    v <- rep(0, times=nrow(gene.data))
    v[which(gene.data$chrom==chr)] <- 1
    return(v)
  }))
  cov <- cbind(cov, chrom.dummy.mat)
  # rownames(cov) <- gene.data$gene
  return(cov)
}

#Fit model for predicting # of rare functional SV
fitConstraintModel <- function(cov, train.genes=which(gene.data$ptv_oe_dec>=5 & gene.data$ptv_oe_dec<=9)){
  #Fit Poisson model
  # outlier.gene.cutoff <- quantile(cov$y.SV_count, 0.99)
  # outlier.gene.cutoff <- 3
  # train.genes <- intersect(train.genes, which(cov$y.SV_count<=outlier.gene.cutoff))
  # glm.fit <- glm(y.SV_count ~ ., data=cov[train.genes, -1], family="poisson")
  glm.fit <- glm.nb(y.SV_count ~ ., data=cov[train.genes, -1])
  #Evaluate variance explained
  # af <- anova(glm.fit)
  # af$PctExp <- 100*af$`Sum Sq`/sum(af$`Sum Sq`)
  # pct.expl <- sum(af$PctExp[-length(af$PctExp)])
  #Apply model to all genes
  glm.fit.vals <- predict.glm(glm.fit, newdata=cov[, -(1:2)], type="response")
  #Return output data frame
  res <- data.frame("gene"=cov$gene, 
                    "SV_count_raw"=cov$y.SV_count, 
                    "SV_count_glm_exp"=glm.fit.vals)
  return(res)
}

#Make summed dotplots
plotSummedDots <- function(deciles, SV_count_obs, SV_count_exp, 
                           color, title=NULL, ymax=NULL, 
                           xlabel="SNV pLoF Constraint Percentile", 
                           ax.labels=TRUE, cex.labels=0.8, 
                           tck=NULL, yline=-0.4, conf.int=F, parmar=c(2, 3.25, 1.75, 2)){
  require(zoo, quietly=T)
  d <- sort(unique(as.numeric(deciles)))
  #Compute summed obs/exp
  means <- sapply(d, function(i){
    exp.sum <- sum(SV_count_exp[which(deciles==i)])
    obs.sum <- sum(SV_count_obs[which(deciles==i)])
    return(obs.sum/exp.sum)
  })
  means[which(means<0)] <- NA
  #Compute CI
  if(conf.int==T){
    cis <- t(sapply(d, function(i){
      oe.vals <- SV_count_obs[which(deciles==i)]/SV_count_exp[which(deciles==i)]
      get.mean <- function(vals, indices){mean(vals[indices])}
      set.seed(i)
      boot.obj <- boot(data=oe.vals, statistic=get.mean, R=1000)
      ci <- boot.ci(boot.obj, conf=0.9, type="basic")$basic[4:5]
      return(ci)
    }))
  }
  
  #Prep plot area
  if(is.null(ymax)){
    ymax <- max(means, na.rm=T)
  }
  par(mar=parmar)
  plot(x=c(0, length(d)), y=c(0, ymax), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  # axis(1, at=seq(0, 100, 10), labels=NA)
  sapply(seq(0, 100, 20), function(p){
    axis(1, at=p, tick=F, cex.axis=cex.labels, line=-0.9, 
         labels=bquote(.(p)^'th'))
  })
  # text(x=(1:10)-0.75, y=par('usr')[3]-(0.115*(par("usr")[4]-par("usr")[3])), 
  #      labels=paste(seq(0, 90, 10), "-", seq(10, 100, 10), "%", sep=""), 
  #      cex=0.7, srt=45)
  axis(2, at=axTicks(2), labels=NA, tck=tck)
  if(ax.labels==T){
    mtext(1, line=0.9, text=xlabel, cex=cex.labels)
    mtext(2, line=2.15, text="Rare SV Obs/Exp", cex=cex.labels)
  }
  axis(2, at=axTicks(2), line=yline, cex.axis=cex.labels, 
       labels=paste(100*round(axTicks(2), 2), "%", sep=""), tick=F, las=2)
  mtext(3, text=title, line=0.2, cex=0.9)
  abline(h=1, lwd=2, col="gray80")
  
  #Add points & rolling mean/ci
  # points(x=d-0.5, y=means, col=color, type="l")
  if(conf.int==T){
    polygon(x=c(d-0.5, rev(d-0.5)), 
            y=c(rollapply(cis[, 1], 21, mean, na.rm=T, partial=T), 
                rev(rollapply(cis[, 2], 21, mean, na.rm=T, partial=T))), 
            border=NA, bty="n", col=adjustcolor(color, alpha=0.1))
  }
  points(x=d-0.5, y=means, pch=21, col=color, cex=0.4)
  points(x=d-0.5, y=rollapply(means, 21, mean, na.rm=T, partial=T), 
         lwd=2, type="l", col=color)
  
  #Add correlation coefficient
  cor.res <- cor.test(x=d-0.5, y=means, method="spearman")
  text(x=par("usr")[1]-(0.025*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[3]+(0.885*(par("usr")[4]-par("usr")[3])), 
       pos=4, cex=cex.labels, 
       labels=bquote(rho == .(format(round(cor.res$estimate, 2), nsmall=2))))
  # cor.res <- cor.test(x=d-0.5, y=means, method="pearson")
  cor.p <- cor.res$p.value
  if(cor.p>10^-100){
    cor.p <- format(cor.p, scientific=T)
    cor.p.base <- format(as.numeric(strsplit(cor.p, split="e")[[1]][1]), nsmall=2)
    cor.p.exp <- format(as.numeric(strsplit(cor.p, split="e")[[1]][2]), nsmall=0)
    text(x=par("usr")[2]+(0.025*(par("usr")[2]-par("usr")[1])), 
         y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])), 
         pos=2, cex=cex.labels, 
         labels=bquote(italic(P) == .(format(round(as.numeric(cor.p.base), 2), nsmall=2))*"x"*10^.(cor.p.exp)))
  }else{
    text(x=par("usr")[2]+(0.025*(par("usr")[2]-par("usr")[1])), 
         y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])), 
         pos=2, cex=cex.labels, 
         labels=bquote(italic(P) < 10^-100))
  }
  box()
}

#Wrapper for all constraint analyses for a single class of SV
constraintModelWrapper <- function(gene.data, y, snv.class="ptv", color, title, 
                                   ymax=NULL, return=FALSE, ax.labels=TRUE, 
                                   conf.int=F, cex.labels=0.8, 
                                   tck=NULL, yline=-0.4, parmar=c(2, 3.25, 1.75, 2)){
  if(snv.class=="ptv"){
    # train.genes <- which(gene.data$pli<=0.1)
    train.genes <- which(gene.data$ptv_oe_dec>=5 & gene.data$ptv_oe_dec<=10)
    deciles <- gene.data$ptv_oe_cent
    xlabel <- "SNV pLoF Constraint Percentile"
  }else{
    train.genes <- which(gene.data$mis_oe_dec>=5 & gene.data$mis_oe_dec<=10)
    deciles <- gene.data$mis_oe_cent
    xlabel <- "Missense Constraint Percentile"
  }
  cov <- prepCovariates(gene.data=gene.data, y=y)
  res <- fitConstraintModel(cov=cov, train.genes=train.genes)
  print(paste("VARIANCE EXPLAINED: ", 100*(cor(res[, 2], res[, 3], use="complete.obs")^2), "%", sep=""))
  plotSummedDots(deciles=deciles, 
                 SV_count_obs=res$SV_count_raw, 
                 SV_count_exp=res$SV_count_glm_exp, 
                 color=color, title=title, ymax, 
                 xlabel=xlabel, 
                 ax.labels=ax.labels, 
                 cex.labels=cex.labels, 
                 conf.int=conf.int, 
                 parmar=parmar, 
                 tck=tck, yline=yline)
  if(return==T){
    return(res)
  }
}

# #Plot per-gene SV-only constraint data (boxplot)
# 
# #Wrapper for SV-only constraint analysis (per-gene significance)
# perGeneSVConstraintWrapper <- function(gene.data, y, color, return=F){
#   #Gather data
#   train.genes <- which(gene.data$ptv_oe_dec>=5 & gene.data$ptv_oe<10)
#   cov <- prepCovariates(gene.data=gene.data, y=y)
#   res <- fitConstraintModel(cov=cov, train.genes=train.genes)
#   res$oe <- res[, 2]/res[, 3]
#   res$p <- apply(res[, 2:3], 1, function(vals){
#     obs <- as.numeric(vals[1])
#     exp <- as.numeric(vals[2])
#     ppois(q=obs, lambda=exp)
#   })
#   res$q <- p.adjust(res$p, method="fdr")
#   res$bonf <- p.adjust(res$p, method="bonf")
#   #Plot data
#   
# }

#Mini histogram of # of SV per gene by constraint category
minihist.byConstraint <- function(gene.data, min.pli=0, max.pli=1, 
                                  title=NULL, x.max=10){
  #Get hist data
  counts <- sapply(0:x.max, function(k){
    length(which(gene.data[, which(colnames(gene.data) == "rare.lof.any")]==k
                 & gene.data$pli>=min.pli & gene.data$pli<=max.pli))
  })
  counts <- c(counts, length(which(gene.data[, which(colnames(gene.data) == "rare.lof.any")]>x.max
                                   & gene.data$pli>=min.pli & gene.data$pli<=max.pli)))
  counts <- counts/1000
  #Plot hist
  base.col <- "gray70"
  par(mar=c(2, 3, 2.5, 0.5), bty="n")
  plot(x=c(0, x.max+2), y=c(0, max(counts)), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=(0:(x.max+1))+0.1, xright=c(0:(x.max+1))+0.9, 
       ybottom=0, ytop=counts, col=base.col)
  sapply((0:x.max), function(k){
    axis(1, at=k+0.5, labels=k, cex.axis=0.8, tick=F, line=-1.3)
  })
  axis(1, at=x.max+1.5, labels=paste("> ", x.max, sep=""), cex.axis=0.8, tick=F, line=-1.3)
  mtext(1, text="Rare pLoF SV per Gene", line=0.75, cex=0.75)
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, labels=paste(axTicks(2), "k", sep=""), 
       line=-0.4, las=2, cex.axis=0.8)
  mtext(2, text="Genes", line=1.75, cex=0.75)
  mtext(3, text=title, cex=0.8)
  text(x=mean(c(par("usr")[1], par("usr")[2])), 
       y=par("usr")[4], pos=1, cex=0.7, 
       labels=paste("n = ", prettyNum(length(which(gene.data$pli>=min.pli & gene.data$pli<=max.pli)), 
                                      big.mark=","), 
                    " genes", sep=""))
  text(x=mean(c(par("usr")[1], par("usr")[2])), 
       y=0.9*par("usr")[4], pos=1, cex=0.7, 
       labels=paste("~", format(round(mean(gene.data$rare.lof.any[which(gene.data$pli>=min.pli & gene.data$pli<=max.pli)]), 2), nsmall=2), 
                    " pLoF SV / gene", sep=""))
}


##########################
###HUMAN KNOCKOUT ANALYSIS
##########################
#Gather all data for KO analysis
getKOdata <- function(dat, gene.data){
  #Subset to biallelic, autosomal, resolved pLoF SV
  subdat <- dat[which(!(dat$chrom %in% c("X", "Y")) & !is.na(dat$PROTEIN_CODING__LOF)), ]
  rows.to.exclude <- grep("MULTIALLELIC", subdat$FILTER, fixed=T)
  if(length(rows.to.exclude)>0){
    subdat <- subdat[-rows.to.exclude, ]
  }
  
  #Helper function to return total number of nonredundant genes and total number of variants
  countKO <- function(subdat, gene.data){
    n.variants <- nrow(subdat)
    f.variants <- n.variants/nrow(dat[which(!(dat$chrom %in% c("X", "Y"))), ])
    genes <- unique(sort(unlist(strsplit(subdat$PROTEIN_CODING__LOF, split=","))))
    genes <- genes[which(genes %in% gene.data$gene)]
    n.genes <- length(genes)
    f.genes <- n.genes/nrow(gene.data[which(!(gene.data$chrom %in% c("chrX", "chrY"))), ])
    return(c(n.variants, f.variants, n.genes, f.genes))
  }
  
  #Gather counts
  a <- countKO(subdat, gene.data)
  subdat <- subdat[which(subdat$N_HOMALT>0), ]
  b <- countKO(subdat, gene.data)
  subdat <- subdat[which(subdat$POPMAX_AF<0.01), ]
  c <- countKO(subdat, gene.data)
  subdat <- subdat[which(subdat$N_HOMALT>1), ]
  d <- countKO(subdat, gene.data)
  
  #Prep & return output data frame
  out.res <- t(data.frame(a, b, c, d))
  colnames(out.res) <- c("n.variants", "f.variants", "n.genes", "f.genes")
  rownames(out.res) <- c("All pLoF SV", "pLoF SV with > 0 Homozygotes", 
                         "Rare (AF < 1%) pLoF SV with > 0 Homozygotes", 
                         "Rare (AF < 1%) pLoF SV with > 1 Homozygote")
  return(out.res)
}
#Human KO plot
plotKO <- function(KO.data){
  #Helper function for generic barplot
  genericBarplot <- function(vals, label.format="count"){
    par(mar=c(0.5, 0.5, 2.5, 1.75))
    vals <- rev(vals)
    plot(x=c(0, 1.1*max(vals)), y=c(0, length(vals)), type="n", 
         xlab="", xaxt="n", ylab="", yaxt="n")
    rect(xleft=0, xright=vals, ybottom=(1:length(vals))-0.8, ytop=(1:length(vals))-0.2, 
         col=rev(c(adjustcolor(svtypes$color[which(svtypes$svtype=="DEL")], alpha=0.4), 
                   adjustcolor(svtypes$color[which(svtypes$svtype=="DEL")], alpha=0.6), 
                   adjustcolor(svtypes$color[which(svtypes$svtype=="DEL")], alpha=0.8), 
                   svtypes$color[which(svtypes$svtype=="DEL")])), 
         lwd=0.75)
    abline(v=0, lwd=0.8)
    if(length(axTicks(3))>5){
      x.at <- axTicks(3)[seq(1, length(axTicks(3)), 2)]
      x.at <- c(x.at, x.at[length(x.at)]+(x.at[2]-x.at[1]))
    }else{
      x.at <- axTicks(3)
    }
    if(label.format=="count"){
      text(x=vals-(0.05*(par("usr")[2]-par("usr")[1])), y=(1:length(vals))-0.55, pos=4, cex=0.7, 
           labels=prettyNum(vals, big.mark=","), xpd=T)
      axis(3, at=x.at, labels=NA, tck=-0.05, lwd=0.75)
      axis(3, at=x.at, tick=F, cex.axis=0.7, line=-0.6, 
           labels=paste(format(round(x.at/1000, 1), nsmall=1), "k", sep=""))
    }else{
      text(x=vals-(0.05*(par("usr")[2]-par("usr")[1])), y=(1:length(vals))-0.55, pos=4, cex=0.7, 
           labels=paste(format(round(100*vals, 1), nsmall=1), "%", sep=""), xpd=T)
      axis(3, at=x.at, labels=NA, tck=-0.05, lwd=0.75)
      axis(3, at=x.at, tick=F, cex.axis=0.7, line=-0.6, 
           labels=paste(round(100*x.at, 1), "%", sep=""))
    }
  }
  
  #Prep plot layout
  layout(matrix(1:5, nrow=1, byrow=T), 
         widths=c(2, 1, 1, 1, 1))
  #Panel 1: labels
  par(mar=c(0.5, 0.5, 2.5, 1.75), bty="n")
  plot(x=c(0, 1), y=c(0, nrow(KO.data)), type="n", 
       xlab="", xaxt="n", ylab="", yaxt="n")
  text(x=1.2*par("usr")[2], y=(1:nrow(KO.data))-0.55, pos=2, cex=0.75, 
       labels=rev(rownames(KO.data)), xpd=T)
  #Panel 2: count of variants
  genericBarplot(KO.data[, 1], label.format="count")
  mtext(3, line=1.25, text="SV (Count)", cex=0.7)
  #Panel 3: pct of variants
  genericBarplot(KO.data[, 2], label.format="pct")
  mtext(3, line=1.25, text="SV (Pct.)", cex=0.7)
  #Panel 4: count of genes
  genericBarplot(KO.data[, 3], label.format="count")
  mtext(3, line=1.25, text="Genes (Count)", cex=0.7)
  #Panel 5: pct of genes
  genericBarplot(KO.data[, 4], label.format="pct")
  mtext(3, line=1.25, text="Genes (Pct.)", cex=0.7)
}
#Gather genes from 3-way comparison of existing KO studies
intersect.KOstudies <- function(SV.kos, narasimhan.kos, saleheen.kos){
  #Get intersections
  a.set <- SV.kos
  b.set <- narasimhan.kos
  c.set <- saleheen.kos
  a <- a.set[which(!(a.set %in% b.set) & !(a.set %in% c.set))]
  b <- b.set[which(!(b.set %in% a.set) & !(b.set %in% c.set))]
  c <- c.set[which(!(c.set %in% a.set) & !(c.set %in% b.set))]
  ab <- a.set[which((a.set %in% b.set) & !(a.set %in% c.set))]  
  ac <- a.set[which(!(a.set %in% b.set) & (a.set %in% c.set))]
  bc <- b.set[which(!(b.set %in% a.set) & (b.set %in% c.set))]
  abc <- a.set[which((a.set %in% b.set) & (a.set %in% c.set))]
  #Return list
  return(list("sv"=a, "nara"=b, "sale"=c, 
              "sv_nara"=ab, "sv_sale"=ac, "nara_sale"=bc, 
              "sv_nara_sale"=abc))
}
#Upset plot of three-way comparison vs existing KO studies
plot.KOstudyComparison <- function(ko.intersections){
  #Prep layout
  layout(matrix(1:2, nrow=2), heights=c(3, 1))
  #Top panel: upset bars
  h <- rev(sort(unlist(lapply(ko.intersections, length))))
  par(mar=c(0.1, 12, 0.5, 0.5), bty="n")
  plot(x=c(0, length(ko.intersections)), y=c(0, 1.05*max(h)), type="n", 
       xlab="", xaxt="n", ylab="", yaxt="n")
  col.bars <- rep("gray70", times=length(h))
  col.bars[grep("sv", names(h))] <- svtypes$color[which(svtypes$svtype=="DEL")]
  rect(xleft=(1:length(h))-0.8, xright=(1:length(h))-0.2, 
       ybottom=0, ytop=h, col=col.bars)
  text(x=(1:length(h))-0.5, y=h-(0.025*(par("usr")[4]-par("usr")[3])), pos=3, cex=0.7, 
       labels=prettyNum(h, big.mark=","))
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, line=-0.4, cex.axis=0.7, las=2, 
       labels=prettyNum(axTicks(2), big.mark=","))
  mtext(2, line=2, text="Genes with\nHomozygous pLoF")
  
  #Bottom panel: points & study labels
  par(mar=c(0.1, 12, 0.1, 0.5), bty="n")
  plot(x=c(0, length(h)), y=c(0, 3), type="n", 
       xlab="", xaxt="n", ylab="", yaxt="n")
  point.cex <- 1.25
  sapply(1:3, function(y){
    points(x=(1:length(h))-0.5, y=rep(y-0.5, length(h)), 
           pch=19, col="gray95", cex=point.cex)
  })
  studies <- rev(c("sv", "nara", "sale"))
  seg.lwd <- 2
  segments(x0=((1:length(h))-0.5)[grep("sv_nara", names(h))], 
           x1=((1:length(h))-0.5)[grep("sv_nara", names(h))], 
           y0=rep(2-0.5, length(h))[grep("sv_nara", names(h))], 
           y1=rep(3-0.5, length(h))[grep("sv_nara", names(h))], 
           lwd=seg.lwd)
  segments(x0=((1:length(h))-0.5)[grep("nara_sale", names(h))], 
           x1=((1:length(h))-0.5)[grep("nara_sale", names(h))], 
           y0=rep(1-0.5, length(h))[grep("nara_sale", names(h))], 
           y1=rep(2-0.5, length(h))[grep("nara_sale", names(h))], 
           lwd=seg.lwd)
  segments(x0=((1:length(h))-0.5)[intersect(grep("sv", names(h)), grep("sale", names(h)))], 
           x1=((1:length(h))-0.5)[intersect(grep("sv", names(h)), grep("sale", names(h)))], 
           y0=rep(1-0.5, length(h))[intersect(grep("sv", names(h)), grep("sale", names(h)))], 
           y1=rep(3-0.5, length(h))[intersect(grep("sv", names(h)), grep("sale", names(h)))], 
           lwd=seg.lwd)
  sapply(1:3, function(i){
    if(studies[i]=="sv"){
      pt.col <- svtypes$color[which(svtypes$svtype=="DEL")]
    }else{
      pt.col <- "black"
    }
    points(x=((1:length(h))-0.5)[grep(studies[i], names(h))], 
           y=rep(i-0.5, length(h))[grep(studies[i], names(h))], 
           pch=21, cex=point.cex, bg=pt.col)
  })
  axis(2, at=(1:3)-0.5, tick=F, line=-0.9, las=2, cex.axis=0.75, 
       labels=c(expression("Saleheen"~italic("et al., Nature")~"(2017)"), 
                expression("Narasimhan"~italic("et al., Science")~"(2016)"), 
                "gnomAD-SV (this study)"))
}
#Plot boxplots of PTV or missense obs/exp across a list of gene sets
boxplots.oe <- function(genesets.list, gene.data, snv.category="ptv", base.colors, ymax=2){
  #Gather oe data
  oe.idx <- which(colnames(gene.data)==paste(snv.category, "_oe", sep=""))
  oe.vals <- lapply(genesets.list, function(genes){
    oe.vals <- gene.data[which(gene.data$gene %in% genes), oe.idx]
    oe.vals <- oe.vals[which(!is.na(oe.vals))]
    return(oe.vals)
  })
  
  #Prep plot area
  if(snv.category=="ptv"){
    label="pLoF SNV\nObs/Exp"
  }else{
    label="Missense\nObs/Exp"
  }
  par(mar=c(0.5, 4.5, 0.5, 0.5), bty="n")
  plot(x=c(0.75, length(oe.vals)+0.25), y=c(0, ymax), type="n", 
       xlab="", xaxt="n", ylab="", yaxt="n")
  abline(h=median(oe.vals[[1]]), lwd=2, col=base.colors[1])
  boxplot(oe.vals, add=T, outline=F, staplewex=0, lty=1, col=base.colors, 
          xaxt="n", yaxt="n")
  axis(2, at=seq(0, 2, 0.5), labels=NA)
  axis(2, at=seq(0, 2, 0.5), tick=F, las=2, cex.axis=0.7, line=-0.4, 
       labels=paste(seq(0, 200, 50), "%", sep=""))
  mtext(2, line=2.25, text=label, cex=1.2)
}
boxplots.size <- function(genesets.list, gene.data, base.colors){
  #Gather oe data
  size.idx <- which(colnames(gene.data)=="exon_sum")
  size.vals <- lapply(genesets.list, function(genes){
    size.vals <- gene.data[which(gene.data$gene %in% genes), size.idx]
    size.vals <- size.vals[which(!is.na(size.vals))]
    return(size.vals)
  })
  
  #Prep plot area
  label <- "Gene CDS\nLength (kb)"
  ymax <- 1.02*max(boxplot(size.vals, plot=F)$stats, na.rm=T)
  par(mar=c(0.5, 4.5, 0.5, 0.5), bty="n")
  plot(x=c(0.75, length(size.vals)+0.25), y=c(0, ymax), type="n", 
       xlab="", xaxt="n", ylab="", yaxt="n")
  abline(h=median(size.vals[[1]]), lwd=2, col=base.colors[1])
  boxplot(size.vals, add=T, outline=F, staplewex=0, lty=1, col=base.colors, 
          xaxt="n", yaxt="n")
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.7, line=-0.4, 
       labels=round(axTicks(2)/1000, 0))
  mtext(2, line=2.25, text=label, cex=1.2)
}
getGeneListFraction <- function(comparison.set, genesets.list){
  #Helper: compute fraction of genes in list X that are present in list Y
  XinY <- function(X, Y){
    length(which(X %in% Y))/length(X)
  }
  vals <- unlist(lapply(genesets.list, XinY, Y=comparison.set))
  #Normalize to all genes (assumes all genes is first list
  vals <- vals/vals[1]
  return(vals)
}
getGeneListFractionMultiple <- function(comparison.sets.list, genesets.list, 
                                        comparison.set.names, geneset.names){
  res <- t(do.call("rbind", lapply(comparison.sets.list, getGeneListFraction, genesets.list=genesets.list)))
  colnames(res) <- comparison.set.names
  rownames(res) <- geneset.names
  return(res)
}
#Generate dotplots of KO gene set enrichments
dotplots.KOenrichments <- function(KO.enrichment.data, colors){
  #Round down enrichments to max of 2
  KO.enrichment.data[which(KO.enrichment.data>2)] <- 2
  #Prep plot area
  par(mar=c(0.5, 7, 2.5, 1), bty="n", xpd=F)
  plot(x=c(0, 2), y=c(0, -ncol(KO.enrichment.data)), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  #Add bars
  buffer <- 0.2
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=-(1:ncol(KO.enrichment.data))+buffer, 
       ytop=-(1:ncol(KO.enrichment.data))+1-buffer, 
       border="gray90", col="gray95")
  abline(v=1)
  sapply(1:ncol(KO.enrichment.data), function(i){
    top <- -i+1-buffer; bottom <- -i+buffer
    breaks <- rev(seq(bottom, top, abs(top-bottom)/nrow(KO.enrichment.data)))
    # rect(xleft=0, xright=KO.enrichment.data[, i], 
    #      ybottom=breaks[-length(breaks)], ytop=breaks[-1], 
    #      col=colors)
    points(x=KO.enrichment.data[, i], 
           y=(breaks[-length(breaks)]+breaks[-1])/2, 
           bg=colors, pch=21)
    axis(2, at=-i+0.5, tick=F, line=-0.8, las=2, cex.axis=0.8, 
         labels=colnames(KO.enrichment.data)[i])
    axis(4, at=c(top, bottom), labels=NA, tck=0.02)
    axis(2, at=c(top, bottom), labels=NA, tck=0.02)
  })
  axis(3, at=seq(0, 2, 0.5), labels=NA)
  axis(3, at=seq(0, 2, 0.5), tick=F, labels=paste(seq(0, 200, 50), "%", sep=""), 
       cex.axis=0.75, line=-0.5)
  mtext(3, line=1.4, text="Gene Set Enrichment", cex=1.2)
}


############################
###RARE LOF CARRIER ANALYSIS
############################
#Gather fraction of individuals from a given population with a rare pLoF in a given gene list
countRareLoFCarrierRate <- function(pop="ALL", genelist, dat, maxAF=0.001, mode="ALL"){
  #Get index of all variants with pLoF of at least one gene in list
  lof.in.genelist <- which(unlist(lapply(strsplit(dat$PROTEIN_CODING__LOF, split=","), function(genes){any(genes %in% genelist)})))
  subdat <- dat[lof.in.genelist, ]
  #Set filtering indexes
  if(pop=="ALL"){
    prefix <- ""
  }else{
    prefix <- paste(pop, "_", sep="")
  }
  AF.idx <- which(colnames(subdat)==paste(prefix, "AF", sep=""))
  genos.idx <- which(colnames(subdat)==paste(prefix, "N_BI_GENOS", sep=""))
  het.idx <- which(colnames(subdat)==paste(prefix, "N_HET", sep=""))
  homalt.idx <- which(colnames(subdat)==paste(prefix, "N_HOMALT", sep=""))
  #Count non-ref individuals by SVTYPE
  svtypes.forCarrierAnalysis <- c("DEL", "DUP", "INS", "INV", "CPX")
  fracs <- sapply(svtypes.forCarrierAnalysis, function(svtype){
    hits <- subdat[which(subdat[, AF.idx]<maxAF & subdat$SVTYPE==svtype), ]
    if(mode=="ALL"){
      sum(hits[, het.idx]+hits[, homalt.idx])/max(hits[, genos.idx], na.rm=T)
    }else if(mode=="HET"){
      sum(hits[, het.idx])/max(hits[, genos.idx], na.rm=T)
    }else if(mode=="HOM"){
      sum(hits[, homalt.idx])/max(hits[, genos.idx], na.rm=T)
    }else{
      stop("mode must be one of HOM, HET, ALL")
    }
  })
  return(fracs)
}
#Plot rare pLoF carrier rates for a list of gene lists
plot.rareLoFCarrierByList <- function(pop="ALL", genelists.list, dat, genelist.labels=NULL, barlabels=T, small=F, xmax=NULL, xlabel=NULL){
  res <- do.call("rbind", lapply(genelists.list, countRareLoFCarrierRate, pop=pop, dat=dat))
  #Prep plot area
  plot.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                   svtypes$color[which(svtypes$svtype=="DUP")], 
                   svtypes$color[which(svtypes$svtype=="INS")], 
                   svtypes$color[which(svtypes$svtype=="INV")], 
                   svtypes$color[which(svtypes$svtype=="CPX")])
  par(bty="n")
  if(is.null(xmax)){
    xmax <- max(apply(res, 1, sum))
  }
  plot(x=c(0, 1.1*xmax), y=c(0, -nrow(res)), type="n", 
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  #Plot bars
  sapply(1:nrow(res), function(i){
    rect(xleft=c(0, cumsum(res[i, ])[-ncol(res)]), 
         xright=cumsum(res[i, ]), 
         ybottom=-i+0.2, ytop=-i+0.8, 
         bty="n", border=NA, col=plot.colors)
    rect(xleft=0, xright=sum(res[i, ]), 
         ybottom=-i+0.2, ytop=-i+0.8, 
         col=NA)
    if(barlabels==T){
      if(small==T){
        text(x=sum(res[i, ]), y=-i+0.5, cex=0.6, 
             pos=4, labels=paste(round(100*sum(res[i, ]), 1), "%", sep=""))
      }else{
        text(x=sum(res[i, ]), y=-i+0.5, cex=0.8, 
             pos=4, labels=paste(round(100*sum(res[i, ]), 1), "%", sep=""))
      }
    }
    if(!is.null(genelist.labels)){
      if(small==T){
        axis(2, at=-i+0.5, tick=F, line=-0.8, 
             labels=genelist.labels[i], cex.axis=0.7, las=2)
      }else{
        axis(2, at=-i+0.5, tick=F, line=-0.8, 
             labels=genelist.labels[i], cex.axis=0.85, las=2)
      }
    }
  })
  #Clean up
  if(small!=T){
    axis(1, at=axTicks(1), labels=NA)
    axis(1, at=axTicks(1), tick=F, line=-0.5, cex.axis=0.7, 
         labels=paste(round(100*axTicks(1), 1), "%", sep=""))
    if(is.null(xlabel)){
      xlabel <- "Very Rare (AF<0.1%) SV pLoF Carrier Rate"
    }
    mtext(1, line=1.5, text=xlabel)
  }
}

#Plot rare pLoF carrier rates for a list of gene lists, split by population
plot.rareLoFCarrierByListAndPop <- function(pops, dat, genelists.list, barlabels=T, xmax=NULL, xlabel=NULL){
  res <- lapply(genelists.list, function(g){
    sapply(pops, function(pop){
      countRareLoFCarrierRate(pop=pop, genelist=g, dat=dat)
    })
  })
  #Prep plot area
  plot.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                   svtypes$color[which(svtypes$svtype=="DUP")], 
                   svtypes$color[which(svtypes$svtype=="INS")], 
                   svtypes$color[which(svtypes$svtype=="INV")], 
                   svtypes$color[which(svtypes$svtype=="CPX")])
  if(is.null(xmax)){
    xmax <- max(unlist(lapply(res, function(df){apply(df, 2, sum)})))
  }
  par(bty="n", mfrow=c(length(res), 1), mar=c(1.5, 2, 1.5, 1))
  sapply(1:length(res), function(i){
    plot(x=c(0, 1.125*xmax), y=c(0, -ncol(res[[i]])), type="n", 
         xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
    #Plot bars
    plot.df <- res[[i]]
    sapply(1:ncol(plot.df), function(k){
      rect(xleft=c(0, cumsum(plot.df[, k])[-nrow(plot.df)]), 
           xright=cumsum(plot.df[, k]), 
           ybottom=-k+0.2, ytop=-k+0.8, 
           bty="n", border=NA, col=plot.colors)
      rect(xleft=0, xright=sum(plot.df[, k]), 
           ybottom=-k+0.2, ytop=-k+0.8, 
           col=NA)
      if(barlabels==T){
        text(x=sum(plot.df[, k]), y=-k+0.5, cex=0.9, 
             pos=4, labels=paste(round(100*sum(plot.df[, k]), 1), "%", sep=""))
      }
      axis(2, at=-k+0.5, tick=F, labels=pops[k], cex.axis=0.9, las=2, line=-0.8)
    })
    axis(1, at=axTicks(1), labels=NA, tck=-0.04)
    axis(1, at=axTicks(1), tick=F, line=-0.75, cex.axis=0.8, 
         labels=paste(round(100*axTicks(1), 1), "%", sep=""))
  })
}

#Plot vertical barplot for rare pLoF carrier rates for a list of gene lists
plot.rareLoFCarrierByList.vertical <- function(pop="ALL", genelists.list, dat, modes, genelist.labels=NULL, 
                                               ymax=NULL, y.break=NULL, ylabel=NULL,
                                               parmar=c(5.75, 3.75, 1, 0.5)){
  res <- do.call("rbind", lapply(1:length(genelists.list), function(i){
    countRareLoFCarrierRate(genelist=genelists.list[[i]], pop=pop, dat=dat, mode=modes[i])
  }))
  #Prep plot area
  plot.colors <- c(svtypes$color[which(svtypes$svtype=="DEL")], 
                   svtypes$color[which(svtypes$svtype=="DUP")], 
                   svtypes$color[which(svtypes$svtype=="INS")], 
                   svtypes$color[which(svtypes$svtype=="INV")], 
                   svtypes$color[which(svtypes$svtype=="CPX")])
  par(bty="n")
  if(is.null(ymax)){
    ymax <- max(apply(res, 1, sum))
  }
  if(!is.null(y.break)){
    ymax <- ymax-(y.break[2]-y.break[1])
  }
  par(mar=parmar, bty="n", xpd=F)
  plot(x=c(0, nrow(res))+1, y=c(0, ymax), type="n", 
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  par(xpd=T)
  rect(xleft=1, xright=par("usr")[2], 
       ybottom=par("usr")[3]-(0.06*(par("usr")[4]-par("usr")[3])), 
       ytop=par("usr")[3]-(0.01*(par("usr")[4]-par("usr")[3])), 
       bty="n", border=NA, col="gray90")
  par(xpd=F)
  abline(h=0)
  #Plot bars
  sapply(1:nrow(res), function(i){
    rect(ybottom=c(0, cumsum(res[i, ])[-ncol(res)]), 
         ytop=cumsum(res[i, ]), 
         xleft=i+0.2, xright=i+0.8, 
         bty="n", border=NA, col=plot.colors)
    rect(ybottom=0, ytop=sum(res[i, ]), 
         xleft=i+0.2, xright=i+0.8, 
         col=NA)
    text(y=sum(res[i, ])-(0.02*(par("usr")[4]-par("usr")[3])), 
         x=i+0.6, cex=0.85, xpd=T, 
         pos=3, labels=paste(round(100*sum(res[i, ]), 2), "%", sep=""))
    if(!is.null(genelist.labels)){
      par(xpd=T)
      text(x=i+0.5, y=par("usr")[3]-(0.035*(par("usr")[4]-par("usr")[3])), 
           labels=prettyNum(length(genelists.list[[i]]), big.mark=","), 
           cex=0.7)
      text(x=i+0.85, y=par("usr")[3]-(0.1*(par("usr")[4]-par("usr")[3])), 
           labels=genelist.labels[i], srt=40, pos=2, cex=0.9)
      par(xpd=F)
    }
  })
  #Clean up bottom panel
  axis(2, at=c(par("usr")[3], ymax), labels=NA, tck=0)
  if(!is.null(y.break)){
    axis(2, at=axTicks(2)[which(axTicks(2)<y.break[1])], labels=NA, tck=-0.025)
    axis(2, at=axTicks(2)[which(axTicks(2)<y.break[1])], tick=F, line=-0.6, cex.axis=0.9, las=2, 
         labels=paste(round(100*axTicks(2)[which(axTicks(2)<y.break[1])], 1), "%", sep=""))
  }
  mtext(2, line=2.3, text=ylabel, cex=1.4)
  #Add top panel, if necessary
  if(!is.null(y.break)){
    #Clear all data above future axis break
    rect(xleft=par("usr")[1]+0.2, xright=par("usr")[2]-0.2, 
         ybottom=y.break[1], 
         ytop=par("usr")[4], 
         col="white", border="white", lwd=4)
    #Add stacked rectangles as required
    sapply(which(apply(res, 1, sum)>y.break[1]), function(i){
      revised.vals <- cumsum(res[i, ])-(y.break[2]-y.break[1])
      rect(xleft=i+0.2, xright=i+0.8, 
           ybottom=c(y.break[1], revised.vals[-length(revised.vals)]), 
           ytop=revised.vals, 
           bty="n", border=NA, col=plot.colors)
      rect(xleft=i+0.2, xright=i+0.8, 
           ybottom=0, ytop=max(revised.vals), 
           col=NA)
      #Add bar labels
      text(y=max(revised.vals)-(0.02*(par("usr")[4]-par("usr")[3])), 
           x=i+0.6, cex=0.85, xpd=T, 
           pos=3, labels=paste(round(100*sum(res[i, ]), 2), "%", sep=""))
    })
    #Add top axis
    axis(2, at=axTicks(2)[which(axTicks(2)>y.break[1])], labels=NA, tck=-0.025)
    axis(2, at=axTicks(2)[which(axTicks(2)>y.break[1])], tick=F, line=-0.6, cex.axis=0.9, las=2, 
         labels=paste(round(100*(axTicks(2)[which(axTicks(2)>y.break[1])]+(y.break[2]-y.break[1])), 1), "%", sep=""))
    #Add axis break
    rect(xleft=par("usr")[1], xright=par("usr")[2], 
         ybottom=y.break[1]-(0.01*(par("usr")[4]-par("usr")[3])), 
         ytop=y.break[1]+(0.01*(par("usr")[4]-par("usr")[3])), 
         col="white", bty="n", border=NA)
    axis.break(2, breakpos=y.break[1], brw=0.04)
    # abline(h=c(y.break[1]-(0.01*(par("usr")[4]-par("usr")[3])), 
    #            y.break[1]+(0.01*(par("usr")[4]-par("usr")[3]))), 
    #        lty=2, lwd=0.75)
  }
}


#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot chromosome-wide SV density for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read density pileups
loadPileup <- function(dat.in,centromeres.in,window=11,smooth=T){
  #Read data
  dat <- read.table(dat.in,header=T,comment.char="")
  colnames(dat)[1] <- "chr"
  cents <- read.table(centromeres.in,header=T,comment.char="")
  colnames(cents)[1] <- "chr"
  #Exclude bins overlapping centromeres
  drop.rows <- unlist(sapply(unique(dat$chr),function(k){
    cent.pos <- as.numeric(cents[which(cents$chr==k),2:3])
    which(dat$chr==k & dat$start<cent.pos[2] & dat$end>cent.pos[1])
  }))
  if(length(drop.rows)>0){
    dat <- dat[-drop.rows,]
  }
  #Sum counts
  dat$ALL <- apply(dat[,-c(1:3)],1,sum,na.rm=T)
  #Smooth counts per chromosome
  require(zoo,quietly=T)
  if(smooth==T){
    for(k in unique(dat$chr)){
      dat[which(dat$chr==k),-c(1:3)] <- apply(dat[which(dat$chr==k),-c(1:3)],2,function(vals){
        rollapply(as.numeric(vals),width=window,FUN=mean,partial=T)
      })
    }
  }
  #Assign normalized distance from centromere
  #Range: -1 = at pter, 0 = in centromere, +1 = at qter
  dat$cdist.norm <- as.numeric(unlist(sapply(unique(dat$chr),function(k){
    midpoints <- (dat$end[which(dat$chr==k)]+dat$end[which(dat$chr==k)])/2
    c.mid <- mean(as.numeric(cents[which(cents$chr==k),2:3]))
    c.dists <- midpoints-c.mid
    plen <- c.mid
    qlen <- max(midpoints)-c.mid
    c.dists.n <- c.dists
    c.dists.n[which(c.dists<0)] <- c.dists.n[which(c.dists<0)]/plen
    c.dists.n[which(c.dists>0)] <- c.dists.n[which(c.dists>0)]/qlen
    c.dists.n[which(c.dists.n< -1)] <- -1
    c.dists.n[which(c.dists.n>1)] <- 1
    return(c.dists.n)
  })))
  return(dat)
}
#Gather meta-chromosome average for a single SVTYPE
metaAverage <- function(dat,SVTYPE="ALL",n.bins=250){
  p.bins <- seq(-1,0,by=1/n.bins)
  q.bins <- seq(0,1,by=1/n.bins)
  col.idx <- which(colnames(dat)==SVTYPE)
  p.means <- sapply(1:(length(p.bins)-1),function(i){
    mean(dat[which(dat$cdist.norm>=p.bins[i] & dat$cdist.norm<p.bins[i+1]),col.idx],na.rm=T)
  })
  q.means <- sapply(1:(length(q.bins)-1),function(i){
    mean(dat[which(dat$cdist.norm>q.bins[i] & dat$cdist.norm<=q.bins[i+1]),col.idx],na.rm=T)
  })
  means <- c(p.means,q.means)
  means.norm <- means/mean(means,na.rm=T)
  out.df <- data.frame("norm.pos"=c(p.bins[-length(p.bins)],q.bins[-1]),
                       "mean"=means,"mean.norm"=means.norm)
  return(out.df)
}
#Split densities into terminal, interstitial, and pericentromeric bins
calc.meanByContext <- function(dat,meta.svtypes,ter.buf=0.05,cen.buf=0.05){
  #Helper function to calculate mean, 95% CI, and p-value that the true mean isn't 1
  get.ci <- function(vals){
    vals <- vals[which(!is.na(vals))]
    k <- 1.96*(sd(vals,na.rm=T)/sqrt(length(vals)))
    p.less <- t.test(vals,mu=1,alternative="less")$p.value
    p.greater <- t.test(vals,mu=1,alternative="greater")$p.value
    return(c(log2(c(mean(vals)-k,mean(vals),mean(vals)+k)),p.less,p.greater))
  }
  res <- lapply(meta.svtypes,function(svtype){
    #Mean-normalize all values
    vals <- as.numeric(dat[,which(colnames(dat)==svtype)])
    vals <- vals/mean(vals,na.rm=T)
    #Calculate stats
    ter.idx <- which(dat$cdist.norm<=-1+ter.buf | dat$cdist.norm>=1-ter.buf)
    cen.idx <- which(dat$cdist.norm>=-cen.buf & dat$cdist.norm<=ter.buf)
    int.idx <- which(!(1:nrow(dat) %in% c(ter.idx,cen.idx)))
    ter.stats <- get.ci(vals[ter.idx])
    int.stats <- get.ci(vals[int.idx])
    cen.stats <- get.ci(vals[cen.idx])
    return(data.frame("ter"=ter.stats,
                      "int"=int.stats,
                      "cen"=cen.stats))
  })
  names(res) <- meta.svtypes
  return(res)
}


#####################
###PLOTTING FUNCTIONS
#####################
#Master function to plot raw coverage from all four samples per contig
generateCovPlotsPerChrom <- function(mat,
                                     colors,
                                     contigs.top,
                                     contigs.middle,
                                     contigs.bottom,
                                     labels.on.top=F,
                                     fill=T,norm=F){
  #Normalize data, if optioned
  if(norm==T){
    mat <- lapply(mat,function(chr.mat){
      chr.mat[,-c(1:3)] <- apply(as.data.frame(chr.mat[,-c(1:3)]),2,function(vals){
        vals <- vals/mean(vals,na.rm=T)
        return(vals)
      })
      return(chr.mat)
    })
  }
  
  #Set spacer
  spacer <- 30000000
  
  #Determine total length to be plotted
  chr.lengths <- as.numeric(unlist(lapply(mat,function(chr.mat){
    return(max(chr.mat[,3]))
  })))
  #Add spacer distance between contigs
  chr.lengths[1:22] <- chr.lengths[1:22]+spacer
  
  #Mini helper function to plot a row of summary coverage values
  plotMiniCovSummary <- function(contigs,ymax=NULL){
    #Prep plot area
    if(labels.on.top==T){
      par(mar=c(0.2,2,1.5,1),bty="n",bg="white")
      lpos=3
    }else{
      par(mar=c(1.5,2,0.2,1),bty="n",bg="white")
      lpos=1
    }
    if(is.null(ymax)){
      ymax <- as.integer(ceiling(quantile(unlist(lapply(mat,function(chr.mat){chr.mat[,4]})),probs=0.995)))
    }
    plot(x=c(-0.01*sum(chr.lengths[contigs]),sum(chr.lengths[contigs])-spacer),
         y=c(-0.1*ymax,1.1*ymax),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i")
    #Add background shading rectangle
    rect(xleft=par("usr")[1],xright=par("usr")[2],
         ybottom=par("usr")[3],ytop=par("usr")[4],
         border="gray99",col="gray99")
    box(which="plot",col="white",lwd=3)
    #Add contig positions & labels
    sapply(1:length(contigs),function(i){
      axis(lpos,at=c(c(0,cumsum(chr.lengths[contigs]))[i],
                     c(0,cumsum(chr.lengths[contigs])-spacer)[i+1]),
           labels=NA,col="gray30",tck=0,line=0.1)
      axis(lpos,at=mean(c(c(0,cumsum(chr.lengths[contigs]))[i],
                          c(0,cumsum(chr.lengths[contigs])-spacer)[i+1])),
           tick=F,labels=paste("chr",contigs[i],sep=""),line=-0.5)
    })
    #Add y axis & gridlines
    y.at <- axTicks(2)[seq(1,length(axTicks(2)),by=2)]
    # y.at <- c(y.at,max(y.at)+(y.at[2]-y.at[1]))
    axis(2,at=axTicks(2),labels=NA,tck=-0.05,col="gray40")
    axis(2,at=y.at,labels=NA,tck=-0.1)
    axis(2,at=y.at,tick=F,las=2,cex.axis=1,line=-0.4)
    # mtext(2,line=2.25,text="Copy Number")
    abline(h=axTicks(2),col="gray80",lwd=0.7,lty=3)
    # abline(h=y.at,col="gray80")
    abline(h=0)
    #Add coverage values
    sapply(1:length(contigs),function(i){
      sapply(1:(ncol(mat[[1]])-3),function(s){
        # cov.vals <- rollmean(cov[[s]][[contigs[i]]][,4],k=11,na.pad=T)
        vals <- mat[[contigs[i]]][,s+3]
        plot.vals <- smooth.spline(x=mat[[contigs[i]]][,2]+c(0,cumsum(chr.lengths[contigs]))[i],
                                   y=vals,spar=0.32*mean(c(1,rep(chr.lengths[contigs[1]]/chr.lengths[contigs[i]],4))))
        if(fill==T){
          polygon(x=c(plot.vals$x,rev(plot.vals$x)),
                  y=c(plot.vals$y,rep(0,times=length(plot.vals$y))),
                  border=NA,col=adjustcolor(colors[s],alpha=0.7))
        }
        points(plot.vals$x,plot.vals$y,type="l",lwd=1.25,col=colors[s])
      })
    })
    #Add cleanup rectangles
    rect(xleft=par("usr")[1],xright=par("usr")[2],
         ybottom=c(par("usr")[3],ymax),ytop=c(0,par("usr")[4]),
         border="white",lty=1,lwd=1,col="white")
    abline(h=par("usr")[3:4],lwd=1,col="white")
    # abline(h=c(0,ymax),col="gray80")
    box(lwd=2,col="white")
    rect(xleft=(cumsum(chr.lengths[contigs])-spacer),
         xright=cumsum(chr.lengths[contigs]),
         ybottom=par("usr")[3],ytop=par("usr")[4],
         border="white",lty=0,lwd=0,col="white")
    rect(xleft=c(par("usr")[1],tail((cumsum(chr.lengths[contigs])-spacer),1)),
         xright=c(0,par("usr")[2]),
         ybottom=par("usr")[3],ytop=par("usr")[4],
         border="white",lty=0,lwd=0,col="white")
    #Add N-mask rectangles
    sapply(1:length(contigs),function(i){
      rect(xleft=Nmasks[which(Nmasks[,1]==contigs[i]),2]+c(0,cumsum(chr.lengths[contigs]))[i],
           xright=Nmasks[which(Nmasks[,1]==contigs[i]),3]+c(0,cumsum(chr.lengths[contigs]))[i],
           ybottom=0,ytop=ymax,
           border=NA,col="gray80")
      sapply(1:nrow(Nmasks[which(Nmasks[,1]==contigs[i]),]),function(k){
        if((Nmasks[which(Nmasks[,1]==contigs[i]),3]-Nmasks[which(Nmasks[,1]==contigs[i]),2])[k]>10000000){
          text(x=mean(c(Nmasks[which(Nmasks[,1]==contigs[i]),2][k]+c(0,cumsum(chr.lengths[contigs]))[i],
                        Nmasks[which(Nmasks[,1]==contigs[i]),3][k]+c(0,cumsum(chr.lengths[contigs]))[i])),
               y=ymax/2,labels="N",col="gray72",font=3,cex=1.2)
        }
      })
    })
  }
  
  #####Plot stacked panels
  par(mfrow=c(3,1))
  plotMiniCovSummary(contigs.top)
  plotMiniCovSummary(contigs.middle)
  plotMiniCovSummary(contigs.bottom)
}
###Plot meta-chromosome density for a single SVTYPE
plot.metaDist <- function(meta.dat, color, fill=T, norm=F, xlabel=NULL){
  #Clean meta dat
  meta.dat <- meta.dat[which(!is.na(meta.dat$mean)), ]
  if(norm==T){
    meta.dat$mean <- meta.dat$mean.norm
  }
  plot.vals <- rollapply(data=meta.dat$mean.norm, width=6, mean, partial=T)
  #Set parameters
  ymax <- 1.075*max(max(plot.vals, na.rm=T), 
                    2*mean(plot.vals, na.rm=T))
  lpos <- 1
  #Prep plot area
  par(bty="n", bg="white")
  plot(x=1.015*range(meta.dat$norm.pos), y=c(-0.05*ymax, ymax), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i")
  #Add background shading rectangle
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4], 
       border="gray99", col="gray99")
  box(which="plot", col="white", lwd=3)
  #Add contig positions & labels
  axis(lpos, at=par("usr")[1:2], 
       labels=NA, col="gray30", tck=0, line=0.1)
  if(is.null(xlabel)){
    axis(lpos, at=mean(par("usr")[1:2]), 
         tick=F, labels="Meta-chromosome", line=-0.6)
  }else{
    axis(lpos, at=mean(par("usr")[1:2]), 
         tick=F, labels=xlabel, line=-0.6, cex.axis=0.8)
  }
  #Add y axis & gridlines
  y.at <- axTicks(2)
  y.at <- c(y.at, max(y.at)+(y.at[2]-y.at[1]))
  axis(2, at=y.at, labels=NA, tck=-0.05, col="gray40")
  # axis(2, at=y.at, labels=NA, tck=-0.1)
  axis(2, at=y.at, tick=F, las=2, cex.axis=1, line=-0.4, labels=round(y.at, 2))
  # mtext(2, line=2.25, text="Copy Number")
  abline(h=y.at, col="gray80", lwd=0.7, lty=3)
  # abline(h=y.at, col="gray80")
  # abline(v=0, col="#963231")
  if(norm==T){
    abline(h=1, col="gray50")
  }
  #Add coverage values
  # plot.vals <- smooth.spline(x=meta.dat$norm.pos, y=meta.dat$mean, spar=0.2)
  # plot.vals <- data.frame("x"=meta.dat$norm.pos, 
  # "y"=meta.dat$mean)
  if(fill==T){
    polygon(x=c(meta.dat$norm.pos, rev(meta.dat$norm.pos)), 
            y=c(plot.vals, rep(0, times=length(plot.vals))), 
            border=NA, col=adjustcolor(color, alpha=0.7))
  }
  points(meta.dat$norm.pos, plot.vals, type="l", lwd=1.25, col=color)
  #Cleanup
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=0, 
       bty="n", border=NA, col="white")
  abline(h=0)
}
###Plot point estimates and CIs for ter/int/cen averages by class
plot.metaByContext <- function(dat,meta.svtypes,colors){
  #Get point estimates and CIs
  plot.vals <- calc.meanByContext(dat,meta.svtypes)
  #Prep plot area
  ylims <- c(1.1*min(as.numeric(unlist(lapply(plot.vals,range,na.rm=T)))),
             1.1*max(as.numeric(unlist(lapply(plot.vals,range,na.rm=T)))))
  par(bty="n",mar=c(0.25,2.8,1.25,0.25))
  plot(x=c(0,length(plot.vals)-0.4),y=ylims,type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")
  abline(h=0,col="gray50")
  #Add points per svtype
  sapply(1:length(plot.vals),function(i){
    segments(x0=i-c(1,0.8,0.6),
             x1=i-c(1,0.8,0.6),
             y0=as.numeric(plot.vals[[i]][1,]),
             y1=as.numeric(plot.vals[[i]][3,]),
             col=colors[i],lwd=2)
    points(x=i-c(1,0.8,0.6),
           y=as.numeric(plot.vals[[i]][2,]),
           pch=19,col=colors[i],cex=1.25)
    text(x=i-c(1,0.8,0.6),
         y=as.numeric(plot.vals[[i]][2,]),
         labels=c("T","I","C"),cex=0.6,font=2,col="white")
    #Add category label
    axis(3,at=i-c(1,0.6),tck=0,labels=NA,line=0.1)
    axis(3,at=i-0.8,tick=F,line=-0.9,cex.axis=0.75,labels=meta.svtypes[i],col.axis=colors[i])
    #Add p-values
    par(xpd=T)
    sapply(1:3,function(k){
      if(plot.vals[[i]][4,k]<0.05/(3*length(plot.vals))){
        text(x=(i-c(1,0.8,0.6))[k],
             y=plot.vals[[i]][2,k],
             pos=1,labels="*")
      }
      if(plot.vals[[i]][5,k]<0.05/(3*length(plot.vals))){
        text(x=(i-c(1,0.8,0.6))[k],
             y=plot.vals[[i]][2,k]-(0.04*(par("usr")[4]-par("usr")[3])),
             pos=3,labels="*")
      }
    })
    par(xpd=F)
  })
  
  #Clean up
  axis(2,at=axTicks(2),labels=NA,tck=-0.03)
  sapply(axTicks(2),function(y){
    axis(2,at=y,tick=F,las=2,cex.axis=0.9,line=-0.5,
         labels=bquote("2"^.(y)))
  })
  mtext(2,line=1.75,text="SV Fold-Enrichment")
}

#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Plot per-individual summary data (both sites and functional effects) for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Gather waterfall data from an input data frame
getWaterfallData <- function(dat){
  #Get ancestry data
  dat$pop[which(dat$pop=="SAS")] <- "OTH"
  ancestries <- sort(unique(as.character(dat$pop)))
  ancestries <- c(ancestries[which(ancestries!="OTH")],
                  ancestries[which(ancestries=="OTH")])
  #Get global & ancestry medians and means
  global.medians <- apply(dat[,-c(1:2)],2,median,na.rm=T)
  ancestry.medians <- lapply(ancestries,function(anc){
    apply(dat[which(dat$pop==anc),-c(1:2)],2,median,na.rm=T)
  })
  names(ancestry.medians) <- ancestries
  global.means <- apply(dat[,-c(1:2)],2,mean,na.rm=T)
  ancestry.means <- lapply(ancestries,function(anc){
    apply(dat[which(dat$pop==anc),-c(1:2)],2,mean,na.rm=T)
  })
  names(ancestry.means) <- ancestries
  #Get list of data frames, one df per ancestry, sorted by total N per sample
  ancestry.perSamp <- lapply(ancestries,function(anc){
    subdat <- dat[which(dat$pop==anc),-c(1:2)]
    sums <- apply(subdat,1,sum,na.rm=T)
    subdat[order(-sums),]
  })
  names(ancestry.perSamp) <- ancestries
  #Get number of samples per ancestry
  nsamp.perAncestry <- unlist(lapply(ancestry.perSamp,nrow))
  #Prep output list of values
  out.list <- list("global.medians"=global.medians,
                   "global.means"=global.means,
                   "ancestries"=ancestries,
                   "nsamp.perAncestry"=nsamp.perAncestry,
                   "ancestry.medians"=ancestry.medians,
                   "ancestry.means"=ancestry.means,
                   "ancestry.perSamp"=ancestry.perSamp)
  return(out.list)
}
#Plot waterfall data
plotWaterfall <- function(wdat,colors,ylabel=NULL,y.at.nmax=5,titles=F,MCNV=T,cat.labels=NULL,cex.ylabs=0.8){
  #Set graphical parameters
  MCNV.dens <- 60
  
  #Drop extreme outliers for visualization purposes
  wdat$ancestry.perSamp <- lapply(wdat$ancestry.perSamp,function(subdat){
    totals <- apply(subdat,1,sum,na.rm=T)
    q3 <- quantile(totals,0.75)
    iqr <- IQR(totals)
    outliers <- which(totals>q3+(3*iqr))
    print(length(outliers)/nrow(subdat))
    if(length(outliers)>0){
      subdat <- subdat[-outliers,]
    }
    return(subdat)
  })
  
  #Prep plot area & scaling
  n.samps <- sum(wdat$nsamp.perAncestry)
  x.buffer.small <- 0.02*n.samps
  x.buffer.large <- 0.05*n.samps
  y.max <- 1.025*max(unlist(lapply(wdat$ancestry.perSamp,function(d){apply(d,1,sum,na.rm=T)})))
  n.ancestries <- length(wdat$ancestries)
  xlims <- c(0,n.samps+(n.ancestries*(x.buffer.large+(1.5*x.buffer.small)))+((2*x.buffer.large)+x.buffer.small))
  plot(x=xlims,y=c(0,y.max),type="n",
       xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="",yaxs="i")
  
  #Plot master median of all samples
  rect(xleft=xlims[2]-(x.buffer.large+x.buffer.small),
       xright=xlims[2]-x.buffer.small,
       ybottom=-par("usr")[4],ytop=sum(wdat$global.medians),
       col=NA,lwd=2)
  # segments(x0=xlims[2]-x.buffer.small,x1=xlims[1]+x.buffer.large,
  #          y0=sum(wdat$global.medians),
  #          y1=sum(wdat$global.medians),
  #          col="gray50")
  rect(xleft=xlims[2]-(x.buffer.large+x.buffer.small),
       xright=xlims[2]-x.buffer.small,
       ybottom=c(0,cumsum(wdat$global.medians)[-length(wdat$global.medians)]),
       ytop=cumsum(wdat$global.medians),
       bty="n",border=NA,col=colors)
  if(MCNV==T){
    rect(xleft=xlims[2]-(x.buffer.large+x.buffer.small),
         xright=xlims[2]-x.buffer.small,
         ybottom=c(0,cumsum(wdat$global.medians)[-length(wdat$global.medians)])[grep("MCNV",names(wdat$global.medians),fixed=T)],
         ytop=cumsum(wdat$global.medians)[grep("MCNV",names(wdat$global.medians),fixed=T)],
         bty="n",border=NA,lend="butt",lwd=0.5,col="white",density=MCNV.dens)
  }
  text(x=xlims[2]-((0.5*x.buffer.large)+x.buffer.small),
       y=sum(wdat$global.medians)-(0.025*(par("usr")[4]-par("usr")[3])),
       font=2,pos=3,cex=0.7,xpd=T,
       labels=prettyNum(round(sum(wdat$global.medians),0),big.mark=","))
  #Add title, if optioned
  if(titles==T){
    axis(3,at=c(xlims[2]-(x.buffer.large+x.buffer.small),xlims[2]-x.buffer.small),
         tck=0,labels=NA,line=0.25)
    axis(3,at=mean(c(xlims[2]-(x.buffer.large+x.buffer.small),xlims[2]-x.buffer.small)),
         tick=F,labels="Median",cex.axis=0.7,xpd=T,padj=0,font=2,line=-0.7)
  }
  #Add connector to right axis
  right.axis.at <- seq(par("usr")[3],par("usr")[4],by=(par("usr")[4]-par("usr")[3])/length(wdat$global.medians[which(wdat$global.medians>0)]))
  right.axis.at <- right.axis.at-(0.5*(right.axis.at[2]-right.axis.at[1]))
  right.axis.at <- right.axis.at[-1]
  sapply(1:length(wdat$global.medians[which(wdat$global.medians>0)]),function(i){
    axis(4,at=right.axis.at[i],labels=NA,lwd=0.75,lend="butt")
    axis(4,at=right.axis.at[i],tick=F,cex.axis=0.7,las=2,line=-0.4,
         labels=paste(prettyNum(round(wdat$global.medians[which(wdat$global.medians>0)][i],0),big.mark=","),
                      cat.labels[which(wdat$global.medians>0)][i],sep=" "),
         col.axis=colors[which(wdat$global.medians>0)][i])
    axis(4,at=right.axis.at[i],tick=F,cex.axis=0.7,las=2,line=-0.4,
         labels=prettyNum(round(wdat$global.medians[which(wdat$global.medians>0)][i],0),big.mark=","),
         col.axis="black")
  })
  segments(x0=rep(xlims[2]-x.buffer.small,length(wdat$global.medians[which(wdat$global.medians>0)])),
           x1=rep(par("usr")[2],length(wdat$global.medians[which(wdat$global.medians>0)])),
           y0=as.numeric((c(0,cumsum(wdat$global.medians)[-length(wdat$global.medians)])+cumsum(wdat$global.medians))/2)[which(wdat$global.medians>0)],
           y1=right.axis.at,
           lwd=0.75,lend="butt")
  
  #Iterate over ancestries & plot
  sapply(1:length(wdat$ancestries),function(i){
    #Determine starting coordinate
    if(i>1){
      left.pos <- sum(wdat$nsamp.perAncestry[1:(i-1)])
    }else{
      left.pos <- 0
    }
    left.pos <- left.pos+(x.buffer.large*i)+(1.5*x.buffer.small*i)
    #Plot individual sample bars
    sapply(1:nrow(wdat$ancestry.perSamp[[i]]),function(k){
      rect(xleft=left.pos+k-1,xright=left.pos+k,
           ybottom=c(0,cumsum(as.numeric(wdat$ancestry.perSamp[[i]][k,]))[-ncol(wdat$ancestry.perSamp[[i]])]),
           ytop=cumsum(as.numeric(wdat$ancestry.perSamp[[i]][k,])),
           bty="n",border=NA,col=adjustcolor(colors,alpha=0.5))
      if(MCNV==T){
        rect(xleft=left.pos+k-1,xright=left.pos+k,
             ybottom=c(0,cumsum(as.numeric(wdat$ancestry.perSamp[[i]][k,]))[-ncol(wdat$ancestry.perSamp[[i]])])[grep("MCNV",names(wdat$ancestry.medians[[i]]),fixed=T)],
             ytop=cumsum(as.numeric(wdat$ancestry.perSamp[[i]][k,]))[grep("MCNV",names(wdat$ancestry.medians[[i]]),fixed=T)],
             bty="n",border=NA,lend="butt",lwd=0.5,col="white",density=MCNV.dens)
      }
    })
    #Plot ancestry median
    segments(x0=left.pos-(1.5*x.buffer.small),x1=left.pos+wdat$nsamp.perAncestry[i],
             y0=sum(wdat$ancestry.medians[[i]]),y1=sum(wdat$ancestry.medians[[i]]),
             col="gray40",lwd=0.75)
    rect(xleft=left.pos-(1.5*x.buffer.small),xright=left.pos-(0.5*x.buffer.small),
         ybottom=-par("usr")[4],ytop=sum(wdat$ancestry.medians[[i]]),
         col=NA,lwd=2)
    rect(xleft=left.pos-(1.5*x.buffer.small),xright=left.pos-(0.5*x.buffer.small),
         ybottom=c(0,cumsum(wdat$ancestry.medians[[i]])[-length(wdat$ancestry.medians[[i]])]),
         ytop=cumsum(wdat$ancestry.medians[[i]]),
         bty="n",border=NA,col=colors)
    if(MCNV==T){
      rect(xleft=left.pos-(1.5*x.buffer.small),xright=left.pos-(0.5*x.buffer.small),
           ybottom=c(0,cumsum(wdat$ancestry.medians[[i]])[-length(wdat$ancestry.medians[[i]])])[grep("MCNV",names(wdat$ancestry.medians[[i]]),fixed=T)],
           ytop=cumsum(wdat$ancestry.medians[[i]])[grep("MCNV",names(wdat$ancestry.medians[[i]]),fixed=T)],
           bty="n",border=NA,lend="butt",lwd=0.5,col="white",density=MCNV.dens)
    }
    text(x=left.pos-(1.5*x.buffer.small),
         y=sum(wdat$ancestry.medians[[i]])-(0.025*(par("usr")[4]-par("usr")[3])),
         pos=3,cex=0.7,xpd=T,
         labels=prettyNum(round(sum(wdat$ancestry.medians[[i]]),0),big.mark=","))
    #Add title, if optioned
    if(titles==T){
      axis(3,at=c(left.pos-(1.5*x.buffer.small),left.pos+wdat$nsamp.perAncestry[i]),tck=0,labels=NA,line=0.25)
      # axis(3,at=mean(c(left.pos-(1.5*x.buffer.small),left.pos+wdat$nsamp.perAncestry[i])),
      #      tick=F,labels=paste(pops$name[which(pops$pop==wdat$ancestries[i])],
      #                          "\n(n=",prettyNum(wdat$nsamp.perAncestry[i],big.mark=","),")",sep=""),
      #      cex.axis=0.7,xpd=T,line=-0.8,padj=0)
      axis(3,at=mean(c(left.pos-(1.5*x.buffer.small),left.pos+wdat$nsamp.perAncestry[i])),
           tick=F,labels=pops$pop[which(pops$pop==wdat$ancestries[i])],
           cex.axis=0.7,xpd=T,line=-0.7,padj=0)
    }
  })
  #Clean up
  axis(1,at=c(par("usr")[1],par("usr")[2]),labels=NA,tck=0)
  if(length(axTicks(2))>y.at.nmax){
    y.at <- axTicks(2)[seq(1,length(axTicks(2)),2)]
    y.at <- c(y.at,y.at[length(y.at)]+(y.at[2]-y.at[1]))
  }else{
    y.at <- axTicks(2)
  }
  axis(2,at=y.at,labels=NA)
  axis(2,at=y.at,tick=F,labels=prettyNum(y.at,big.mark=","),
       line=-0.4,cex.axis=cex.ylabs,las=2)
  mtext(2,line=2.5,text=ylabel,cex=0.6)
}

#Generate jitter residuals for a vector of values based on their density
sina.jitter <- function(vals){
  d <- density(vals)
  dv <- approx(d$x,d$y,xout=vals)
  dv <- dv$y/max(dv$y)
  dv.j <- sapply(1:length(vals),function(i){
    jitter(0,amount=dv[i])
  })
  return(dv.j)
}

#Generate sina points to add to existing plot
sina.plot <- function(vals,y.at,color,width=0.1,horiz=T,cex=0.25){
  j <- (width*sina.jitter(vals))+y.at
  if(horiz==T){
    points(x=vals,y=j,pch=21,cex=cex,col=color,bg=adjustcolor(color,alpha=0.3))
    segments(x0=median(vals),x1=median(vals),
             y0=y.at-width,y1=y.at+width,lwd=3,lend="butt")
  }else{
    points(x=j,y=vals,pch=21,cex=cex,col=color,bg=adjustcolor(color,alpha=0.3))
    segments(x0=y.at-width,x1=y.at+width,y0=median(vals),y1=median(vals),lwd=3,lend="butt")
  }
}

#Plot single panel vertical sina plot for a given svtype, by ancestry
plot.sina.singleSvtypeByPop <- function(wdat,svtype,ordered.pops,title=NULL){
  #Get plot values
  plot.vals <- lapply(ordered.pops,function(pop){
    subdf <- wdat$ancestry.perSamp[[which(names(wdat$ancestry.perSamp)==pop)]]
    if(svtype!="ALL"){
      subdf[,grep(paste(".",svtype,sep=""),colnames(subdf),fixed=T)]
    }else{
      apply(subdf,1,sum,na.rm=T)
    }
  })
  if(is.null(title)){
    title <- svtype
  }
  #Prep plot area
  par(mar=c(0.5,2,1.5,0.5),bty="n")
  plot(x=c(0,length(plot.vals)),y=range(unlist(plot.vals),na.rm=T),
       type="n",xlab="",xaxt="n",ylab="",yaxt="n")  
  #Add jitter & color by pop
  sapply(1:length(ordered.pops),function(i){
    pop <- ordered.pops[i]
    pop.col <- pops$color[which(pops$pop==pop)]
    sina.plot(vals=plot.vals[[i]],y.at=i-0.5,color=pop.col,
              width=0.4,horiz=F,cex=0.1)
  })
  #Add axes
  axis(3,at=c(0,length(plot.vals)),tck=0,labels=NA)
  mtext(3,line=0.2,text=title,cex=0.7)
  y.ax.at <- axTicks(2)
  if(length(y.ax.at)>5){
    y.ax.at <- y.ax.at[seq(1,length(y.ax.at),2)]
  }
  axis(2,at=y.ax.at,labels=NA)
  if(max(y.ax.at)>1000){
    axis(2,at=y.ax.at,tick=F,
         labels=paste(round(y.ax.at/1000,2),"k",sep=""),
         las=2,cex.axis=0.8,line=-0.4)
  }else{
    axis(2,at=y.ax.at,tick=F,labels=prettyNum(y.ax.at,big.mark=","),
         las=2,cex.axis=0.8,line=-0.4)
  }
}

#Wrapper to plot horizontal grid of counts of SV per population by SVTYPE
wrapper.gridSinaSvtypeByPop <- function(wdat){
  par(mfrow=c(1,8))
  sapply(c("ALL","DEL","DUP","INV"),function(svtype){
    if(svtype=="ALL"){
      title <- "All SVs"
    }else if(svtype=="MCNV_DEL"){
      title <- "MCNV (Loss)"
    }else if(svtype=="MCNV_DUP"){
      title <- "MCNV (Gain)"
    }else{
      title <- svtype
    }
    plot.sina.singleSvtypeByPop(wdat=wdat,svtype=svtype,
                                ordered.pops=pops$pop,
                                title=title)
  })
}


#Plot single panel vertical sina plot for a set of functional category, restricted by population
plot.sina.MultiFunc <- function(wdat,funcs,func.labels,func.cols,
                                pop="ALL",ymax=NULL,boxes=F,vline=NULL,
                                pt.labels=F){
  # #DEV:
  # wdat <- all.func.wdat.genes
  # funcs <- c("PROTEIN_CODING__LOF.genes",
  #            "PROTEIN_CODING__DUP_LOF.genes",
  #            "PROTEIN_CODING__COPY_GAIN.genes",
  #            "PROTEIN_CODING__MCNV_LOSS.genes",
  #            "PROTEIN_CODING__MCNV_GAIN.genes")
  # func.labels <- c("pLoF","IED","CG","MCNV\n(pLoF)","MCNV\n(CG)")
  # func.cols <- c(svtypes$color[which(svtypes$svtype=="DEL")],
  #                svtypes$color[which(svtypes$svtype=="MCNV")],
  #                svtypes$color[which(svtypes$svtype=="DUP")],
  #                svtypes$color[which(svtypes$svtype=="DEL")],
  #                svtypes$color[which(svtypes$svtype=="DUP")])
  # pop <- "EUR"
  #Get plot values
  if(pop=="ALL"){
    d.tmp <- do.call("rbind",wdat$ancestry.perSamp)
    plot.dat <- as.data.frame(do.call("cbind", lapply(funcs,function(func){
      d.tmp[,grep(func,colnames(d.tmp),fixed=T)]
    })))
  }else{
    d.tmp <- wdat$ancestry.perSamp[[which(names(wdat$ancestry.perSamp)==pop)]]
    plot.dat <- as.data.frame(do.call("cbind", lapply(funcs,function(func){
      d.tmp[,grep(func,colnames(d.tmp),fixed=T)]
    })))
  }
  colnames(plot.dat) <- funcs
  if(is.null(ymax)){
    ymax <- max(apply(plot.dat,2,quantile,probs=0.999))
  }
  #Prep plot area
  par(mar=c(2.5,2.75,0.5,0.5),bty="n")
  plot(x=c(0.1,ncol(plot.dat)+0.25),y=c(0,ymax),
       type="n",xlab="",xaxt="n",ylab="",yaxt="n",yaxs="i")
  if(!is.null(vline)){
    abline(v=vline,lty=2,col="gray80")
  }
  # axis(1,at=c(-100,100))
  #Add jitter & color by class
  sapply(1:ncol(plot.dat),function(i){
    if(boxes==T){
      boxplot(plot.dat[,i],at=i-0.5,col=func.cols[i],
              outline=F,add=T,lty=1,staplewex=0,
              xaxt="n",yaxt="n",xlab="",ylab="")
      if(pt.labels==T){
        text(x=i-0.625,y=median(plot.dat[,i]),pos=4,cex=0.65,xpd=T,
             labels=prettyNum(median(plot.dat[,i]),big.mark=","))
      }
    }else{
      sina.plot(vals=plot.dat[,i],y.at=i-0.5,color=func.cols[i],
                width=0.3,horiz=F,cex=0.025)
      if(pt.labels==T){
        text(x=i-0.525,y=median(plot.dat[,i]),pos=4,cex=0.65,xpd=T,
             labels=prettyNum(median(plot.dat[,i]),big.mark=","))
      }
    }
    axis(1,at=i-0.5,las=2,labels=func.labels[i],srt=2,tick=F,cex.axis=0.85,line=-0.8)
  })
  #Add axes
  y.ax.at <- axTicks(2)
  if(length(y.ax.at)>5){
    y.ax.at <- y.ax.at[seq(1,length(y.ax.at),2)]
  }
  axis(2,at=c(-100,10000))
  axis(2,at=y.ax.at,labels=NA)
  if(max(y.ax.at)>1000){
    axis(2,at=y.ax.at,tick=F,
         labels=paste(round(y.ax.at/1000,2),"k",sep=""),
         las=2,cex.axis=0.8,line=-0.4)
  }else{
    axis(2,at=y.ax.at,tick=F,labels=prettyNum(y.ax.at,big.mark=","),
         las=2,cex.axis=0.8,line=-0.4)
  }
  if(pop!="ALL"){
    mtext(2,line=1.9,text=paste("Genes per",pop,"Genome",sep=" "))
  }else{
    mtext(2,line=1.9,text="Genes per Genome")
  }
}


#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# gnomAD v2 SV analysis script

# Generate VCF-wide distribution plots for formal gnomAD analysis


###Set master parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & clean vcf2bed input
read.vcf2bed <- function(vcf2bed.in, aps.in=NULL){
  #Read data
  dat <- read.table(vcf2bed.in, comment.char="", header=T, sep="\t")
  colnames(dat)[1] <- "chrom"
  #Restrict to sites with at least one observed alternative allele
  dat <- dat[union(which(dat$AC>0), grep("MULTIALLELIC", dat$FILTER, fixed=T)), ]
  #Drop columns not being used (to save memory)
  cols.to.drop <- c("start", "end", "CHR2", "CPX_INTERVALS", 
                    "END", "SOURCE", "STRANDS", "UNRESOLVED_TYPE", 
                    "LINCRNA__LOF", "LINCRNA__DUP_LOF", "LINCRNA__COPY_GAIN", 
                    "LINCRNA__DUP_PARTIAL", "LINCRNA__MSV_EXON_OVR", 
                    "LINCRNA__INTRONIC", "LINCRNA__INV_SPAN", "LINCRNA__UTR")
  dat <- dat[, -which(colnames(dat) %in% cols.to.drop)]
  #Convert numeric columns
  numeric.columns <- sort(unique(c(grep("FREQ", colnames(dat), fixed=T), 
                                   grep("AN", colnames(dat), fixed=T), 
                                   grep("AC", colnames(dat), fixed=T), 
                                   grep("AF", colnames(dat), fixed=T))))
  numeric.columns <- setdiff(numeric.columns, grep("SPAN", colnames(dat), fixed=T))
  dat[, numeric.columns] <- apply(dat[, numeric.columns], 2, as.numeric)
  #Read & add APS, if optioned
  if(!is.null(aps.in)){
    aps <- read.table(aps.in, header=T, sep="\t", comment.char="")
    colnames(aps)[1] <- "VID"
    dat <- merge(dat, aps, by.x="name", by.y="VID", all.x=T, all.y=F, sort=F)
    dat$APS[which(dat$chrom %in% c("X", "Y"))] <- NA
  }
  return(dat)
}

#Plot samples in study vs prior SV studies
plot.sampleSize <- function(priors, pops, pop.assignments){
  ordered.pops <- c("OTH", "SAS", "EUR", "EAS", "AMR", "AFR")
  #Format input data
  plot.dat <- as.data.frame(t(sapply(priors$Study, function(s){
    s.dat <- sapply(ordered.pops, function(p){
      priors[which(priors$Study==s), which(colnames(priors)==paste(p, ".pop", sep=""))]
    })
  })))
  g.dat <- sapply(ordered.pops, function(p){
    length(which(pop.assignments$pop==p))
  })
  # g.dat[which(names(g.dat)=="OTH")] <- g.dat[which(names(g.dat)=="OTH")]+(14891-sum(g.dat))
  plot.dat <- rbind(plot.dat, "gnomAD-SV"=g.dat)
  plot.dat <- apply(plot.dat, 2, as.numeric)
  rownames(plot.dat) <- c(priors$Study, "gnomAD-SV")
  plot.dat <- plot.dat[order(-apply(plot.dat, 1, sum)), ]
  #Prep plot values
  plot.colors <-sapply(ordered.pops, function(p){
    if(p=="SAS"){
      pops$color[which(pops$pop=="OTH")]
    }else{
      pops$color[which(pops$pop==p)]
    }
  })
  #Plot
  par(bty="n")
  plot(x=c(0, 1.15*15000), y=c(0, -nrow(plot.dat)), type="n", 
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="")
  x.at <- axTicks(1)
  if(length(x.at)>5){
    x.at <- x.at[seq(1, length(x.at), 2)]
    x.at <- c(x.at, x.at[length(x.at)]+(x.at[2]-x.at[1]))
  }
  axis(1, at=x.at, labels=NA, tck=-0.03)
  sapply(x.at, function(k){
    axis(1, at=k, tick=F, line=-0.9, cex.axis=0.75, 
         labels=paste(round(k/1000, 0), "k", sep=""))
  })
  mtext(1, text="Sample Size", line=1, cex=0.85)
  sapply(1:nrow(plot.dat), function(i){
    rect(ybottom=-i+0.2, ytop=-i+0.8, 
         xleft=c(0, cumsum(as.numeric(plot.dat[i, -ncol(plot.dat)]))), 
         xright=cumsum(as.numeric(plot.dat[i, ])), 
         border=NA, bty="n", col=plot.colors)
    rect(ybottom=-i+0.2, ytop=-i+0.8, 
         xleft=0, xright=sum(as.numeric(plot.dat[i, ])), 
         col=NA)
    par(xpd=T)
    text(x=sum(as.numeric(plot.dat[i, ]))-(0.05*(par("usr")[2]-par("usr")[1])), 
         y=-i+0.5, pos=4, cex=0.6, 
         labels=prettyNum(sum(as.numeric(plot.dat[i, ])), big.mark=","))
    par(xpd=F)
    if(rownames(plot.dat)[i]=="gnomAD-SV"){
      axis(2, at=-i+(2/3), tick=F, line=-0.8, cex.axis=0.75, las=2, 
           labels="gnomAD-SV")
      axis(2, at=-i+(1/3), tick=F, line=-0.8, cex.axis=0.5, las=2, 
           labels="This study")
    }else{
      year <- paste("(", priors$Year[which(priors$Study==rownames(plot.dat)[i])], ")", sep="")
      axis(2, at=-i+0.5, tick=F, line=-0.8, cex.axis=0.75, las=2, 
           labels=priors$Cohort[which(priors$Study==rownames(plot.dat)[i])])
      # if(rownames(plot.dat)[i]=="HehirKwa"){
      #   axis(2, at=-i+(1/3), tick=F, line=-0.8, cex.axis=0.6, las=2, 
      #        labels=bquote("Hehir-Kwa" ~ italic("et al.") ~  .(year)))
      # }else{
      #   axis(2, at=-i+(1/3), tick=F, line=-0.8, cex.axis=0.6, las=2, 
      #        labels=bquote(.(priors$Author[which(priors$Study==rownames(plot.dat)[i])]) ~ italic("et al.") ~ .(year)))
      # }
    }
  })
}

#Plot total SV sites discovered per study vs prior SV studies
plot.SVcountsByStudy <- function(priors, dat, svtypes){
  #Format input data
  plot.dat <- as.data.frame(t(sapply(priors$Study, function(s){
    s.dat <- sapply(c("DEL", "DUP", "MCNV", "INS", "INV", "CPX", "OTH"), function(p){
      if(p!="CPX"){
        priors[which(priors$Study==s), which(colnames(priors)==p)]
      }else{
        return(0)
      }
    })
  })))
  g.dat <- sapply(c("DEL", "DUP", "MCNV", "INS", "INV", "CPX"), function(p){
    length(which(dat$SVTYPE==p))
  })
  g.dat <- c(g.dat, "OTH"=length(which(!(dat$SVTYPE %in% c("DEL", "DUP", "MCNV", "INS", "INV", "CPX")))))
  plot.dat <- rbind(plot.dat, "gnomAD-SV"=g.dat)
  plot.dat <- apply(plot.dat, 2, as.numeric)
  rownames(plot.dat) <- c(priors$Study, "gnomAD-SV")
  # plot.dat <- plot.dat[order(-apply(plot.dat, 1, sum)), ]
  plot.dat <- plot.dat[match(c("gnomAD-SV", "Sudmant", "HehirKwa", "Chiang"), rownames(plot.dat)), ]
  #Prep plot values
  plot.colors <-sapply(c("DEL", "DUP", "MCNV", "INS", "INV", "CPX", "OTH"), function(p){
    svtypes$color[which(svtypes$svtype==p)]
  })
  #Plot
  par(bty="n")
  plot(x=c(0, 1.15*max(apply(plot.dat, 1, sum))), y=c(0, -nrow(plot.dat)), type="n", 
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="")
  # x.at <- axTicks(1)
  x.at <- seq(0, 1000000, 200000)
  if(length(x.at)>6){
    x.at <- x.at[seq(1, length(x.at), 2)]
    x.at <- c(x.at, x.at[length(x.at)]+(x.at[2]-x.at[1]))
  }
  axis(1, at=x.at, labels=NA, tck=-0.03)
  sapply(x.at, function(k){
    axis(1, at=k, tick=F, line=-0.9, cex.axis=0.75, 
         labels=paste(round(k/1000, 0), "k", sep=""))
  })
  mtext(1, text="SVs Discovered", line=1, cex=0.85)
  sapply(1:nrow(plot.dat), function(i){
    rect(ybottom=-i+0.2, ytop=-i+0.8, 
         xleft=c(0, cumsum(as.numeric(plot.dat[i, -ncol(plot.dat)]))), 
         xright=cumsum(as.numeric(plot.dat[i, ])), 
         border=NA, bty="n", col=plot.colors)
    rect(ybottom=-i+0.2, ytop=-i+0.8, 
         xleft=0, xright=sum(as.numeric(plot.dat[i, ])), 
         col=NA)
    par(xpd=T)
    text(x=sum(as.numeric(plot.dat[i, ]))-(0.05*(par("usr")[2]-par("usr")[1])), 
         y=-i+0.5, pos=4, cex=0.6, 
         labels=prettyNum(sum(as.numeric(plot.dat[i, ])), big.mark=","))
    par(xpd=F)
    if(rownames(plot.dat)[i]=="gnomAD-SV"){
      axis(2, at=-i+(2/3), tick=F, line=-0.8, cex.axis=0.75, las=2, 
           labels="gnomAD-SV")
      axis(2, at=-i+(1/3), tick=F, line=-0.8, cex.axis=0.5, las=2, 
           labels="This study")
    }else{
      year <- paste("(", priors$Year[which(priors$Study==rownames(plot.dat)[i])], ")", sep="")
      axis(2, at=-i+0.5, tick=F, line=-0.8, cex.axis=0.75, las=2, 
           labels=priors$Cohort[which(priors$Study==rownames(plot.dat)[i])])
      # if(rownames(plot.dat)[i]=="HehirKwa"){
      #   axis(2, at=-i+(1/3), tick=F, line=-0.8, cex.axis=0.6, las=2, 
      #        labels=bquote("Hehir-Kwa" ~ italic("et al.") ~  .(year)))
      # }else{
      #   axis(2, at=-i+(1/3), tick=F, line=-0.8, cex.axis=0.6, las=2, 
      #        labels=bquote(.(priors$Author[which(priors$Study==rownames(plot.dat)[i])]) ~ italic("et al.") ~ .(year)))
      # }
    }
  })
}

#Plot simple log-scaled bars of SV by count
plot.totalCounts <- function(dat, svtypes, thousandG=F){
  #Gather data
  counts <- lapply(svtypes$svtype[which(svtypes$svtype!="OTH")], function(svtype){
    return(log10(length(which(dat$SVTYPE==svtype))))
  })
  names(counts) <- svtypes$svtype[which(svtypes$svtype!="OTH")]
  #counts$OTH <- log10(length(which(!(dat$SVTYPE %in% svtypes$svtype[which(svtypes$svtype!="OTH")]))))
  counts <- lapply(counts, function(val){if(is.infinite(val)){val <- 0}; return(val)})
  counts <- unlist(counts)
  sudmant.counts <- log10(c(42279, 6025, 2929, 16631+168, 786, NA, NA))
  names(sudmant.counts) <- names(counts)
  #Plot
  log.minor <- log10(as.numeric(sapply(0:10, function(i){(1:8)*(10^i)})))
  log.major <- 1:9
  par(mar=c(3, 3.5, 3, 0.5), bty="n")
  plot(x=c(0, 2*nrow(svtypes)), y=c(0, 1.02*max(counts)), type="n", 
       xaxt="n", xlab="", yaxt="n", ylab="", yaxs="i")
  if(thousandG==T){
    rect(xleft=seq(1, 2*length(counts), 2)-0.8, xright=seq(1, 2*length(counts), 2), 
         ybottom=0, ytop=sudmant.counts, 
         lwd=0.5, col="gray70")
    rect(xleft=seq(2, 2*length(counts), 2)-1, xright=seq(2, 2*length(counts), 2)-0.2, 
         ybottom=0, ytop=counts, 
         lwd=0.5, col=svtypes$color)
  }else{
    rect(xleft=seq(1, 2*length(counts), 2)-0.8, xright=seq(2, 2*length(counts), 2)-0.2, 
         ybottom=0, ytop=counts, 
         lwd=0.5, col=svtypes$color)
  }
  segments(x0=seq(1, 2*length(counts), 2)-0.8, 
           x1=seq(2, 2*length(counts), 2)-0.2, 
           y0=par("usr")[3], y1=par("usr")[3], lwd=2)
  axis(1, at=seq(1, 2*length(counts), 2), tick=F, line=-0.9, 
       labels=svtypes$svtype, cex.axis=0.8, las=2)
  axis(2, at=log.minor, labels=NA, tck=-0.02, lwd=0.9)
  axis(2, at=log.major, labels=NA, tck=-0.04, lwd=1.1)
  axis(2, at=1:6, tick=F, line=-0.5, cex.axis=0.8, las=2, 
       labels=c("10", "100", "1k", "10k", "100k", "1M"))
  mtext(2, text="SVs Discovered", line=2)
  sapply(1:length(counts), function(i){
    if(thousandG==T){
      if(is.na(sudmant.counts[i])){
        sudmant.y <- 0
        sudmant.lab <- 0
      }else{
        sudmant.y <- sudmant.counts[i]
        sudmant.lab <- 10^sudmant.counts[i]
      }
      text(x=2*i-1.4, y=sudmant.y+(0.025*(par("usr")[4]-par("usr")[3])), 
           labels=prettyNum(sudmant.lab, big.mark=","), 
           srt=90, adj=0, xpd=T, cex=0.7, col="gray60")
      text(x=2*i-0.55, y=counts[i]+(0.025*(par("usr")[4]-par("usr")[3])), 
           labels=prettyNum(10^counts[i], big.mark=","), 
           srt=90, adj=0, xpd=T, cex=0.7)
    }else{
      text(x=2*i-1, y=counts[i]+(0.025*(par("usr")[4]-par("usr")[3])), 
           labels=prettyNum(10^counts[i], big.mark=","), 
           srt=90, adj=0, xpd=T, cex=0.7)
    }
  })
}

#Build summary table of counts,  sizes,  frequencies
build.table <- function(dat.wrelateds, dat.all, dat){
  svtypes.fortable <- c("DEL", "DUP", "MCNV", "INS", "INV", "CPX", "CTX", "BND")
  #Total count of sites
  sites.count <- sapply(svtypes.fortable, function(svtype){
    length(which(dat.wrelateds$SVTYPE==svtype))
  })
  sites.count <- c(sum(sites.count), sites.count)
  #Pct of PASS sites
  pass.count <- sapply(svtypes.fortable, function(svtype){
    length(which(dat.wrelateds$SVTYPE==svtype & (dat.wrelateds$FILTER=="PASS" | dat.wrelateds$FILTER=="MULTIALLELIC")))
  })
  pass.count <- c(sum(pass.count), pass.count)
  pct.pass <- pass.count/sites.count
  #Count of PASS sites in unrelated genomes
  pass.count.unrelated <- sapply(svtypes.fortable, function(svtype){
    length(which(dat.all$SVTYPE==svtype & (dat.all$FILTER=="PASS" | dat.all$FILTER=="MULTIALLELIC")))
  })
  pass.count.unrelated <- c(sum(pass.count.unrelated), pass.count.unrelated)
  #Count of variants by freq bin
  count.byfreq <- t(sapply(svtypes.fortable, function(svtype){
    singletons <- length(which(dat$SVTYPE==svtype & dat$AC==1))
    rare <- length(which(dat$SVTYPE==svtype & dat$AC>1 & dat$AF<0.01))
    common <- length(which(dat$SVTYPE==svtype & (dat$AF>=0.01 | dat$FILTER=="MULTIALLELIC")))
    c(singletons, rare, common)
  }))
  colnames(count.byfreq) <- c("Singleton", "Rare", "Common")
  count.byfreq <- rbind(c(length(which(dat$AC==1)), 
                          length(which(dat$AC>1 & dat$AF<0.01)), 
                          length(which(dat$AF>=0.01))), 
                        count.byfreq)
  #Median size
  med.size <- sapply(svtypes.fortable, function(svtype){
    median(dat$SVLEN[which(dat$SVTYPE==svtype)])
  })
  med.size <- c(median(dat$SVLEN, na.rm=T), med.size)
  med.size[which(med.size<1)] <- NA
  #Format output table
  table.out <- cbind(data.frame("Abbrev."=c("ALL", svtypes.fortable), 
                                "Total Sites Discovered"=sites.count, 
                                "Pct.Pass"=pct.pass, 
                                "Passing Sites"=pass.count.unrelated), 
                     count.byfreq, 
                     data.frame("Med.Size"=med.size))
  return(table.out)
}

# Plot rolling mean size distribution lines
plot.sizes <- function(dat, svtypes, step=0.02, xlims=c(50, 10000000)){
  #Iterate over SVTYPES and get log-scaled count of SV by class per bin
  require(zoo, quietly=T)
  size.bins <- seq(log10(50), log10(50000000), step)
  plot.dat <- as.data.frame(sapply(svtypes$svtype[which(svtypes$svtype!="OTH")], function(svtype){
    res <- sapply(1:length(size.bins), function(i){
      nrow(dat[which(dat$SVLEN>=10^size.bins[i] & dat$SVLEN<10^(size.bins[i]+step) & dat$SVTYPE==svtype), ])
    })
    return(res)
  }))
  plot.dat <- apply(plot.dat, 2, function(vals){rollapply(vals, sum, width=5, partial=T)})
  plot.dat <- apply(plot.dat, 2, log10)
  plot.dat[which(is.infinite(plot.dat))] <- 0
  plot.dat <- as.data.frame(plot.dat)
  
  #Prep plot area
  par(bty="n", mar=c(2, 3.5, 0.5, 1))
  L1.peak <- log10(6000)
  SVA.peak <- log10(1200)
  ALU.peak <- log10(280)
  log.minor <- log10(as.numeric(sapply(0:10, function(i){(1:8)*(10^i)})))
  log.mid <- log10(as.numeric(sapply(0:10, function(i){c(1, 5)*(10^i)})))
  log.major <- 1:9
  plot(x=log10(xlims), y=c(0, 1.05*max(plot.dat)), type="n", 
       xaxt="n", xaxs="i", xlab="", yaxt="n", ylab="", yaxs="i")
  # abline(v=c(L1.peak, SVA.peak, ALU.peak), lty=2, col="gray80")
  # sapply(c(L1.peak, SVA.peak, ALU.peak), function(pos){
  #   axis(3, at=pos, labels=NA)
  # })
  # abline(h=log.minor, v=log.minor, col="gray97", lwd=0.5)
  # abline(h=log.mid, v=log.mid, col="gray94", lwd=0.75)
  # abline(h=log.mid, v=log.major, col="gray91")
  # rect(xleft=c(L1.peak, SVA.peak, ALU.peak)-(0.01*(par("usr")[2]-par("usr")[1])), 
  #      xright=c(L1.peak, SVA.peak, ALU.peak)+(0.01*(par("usr")[2]-par("usr")[1])), 
  #      ybottom=par("usr")[3], ytop=par("usr")[4], 
  #      border=NA, bty="n", col=adjustcolor("gray85", alpha=0.5))
  
  #Plot lines per SV class
  sapply(1:ncol(plot.dat), function(i){
    points(x=size.bins, y=plot.dat[, i], lwd=2.5, col=svtypes$color[i], type="l")
  })
  
  #Add axes
  axis(1, at=log.minor, labels=NA, tck=-0.0175, lwd=0.9)
  # axis(1, at=log.mid, labels=NA, tck=-0.02)
  axis(1, at=log.major, labels=NA, tck=-0.035, lwd=1.1)
  sapply(1:6, function(i){
    axis(1, at=i+1, tick=F, line=-0.9, cex.axis=0.65, 
         labels=c("100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb")[i])
  })
  
  mtext(1, text="SV Size", line=1)
  axis(2, at=log.minor, labels=NA, tck=-0.0175, lwd=0.9)
  # axis(2, at=log.mid, labels=NA, tck=-0.02)
  axis(2, at=log.major, labels=NA, tck=-0.035, lwd=1.1)
  axis(2, at=1:6, tick=F, line=-0.6, cex.axis=0.8, las=2, 
       labels=c("10", "100", "1k", "10k", "100k", "1M"))
  mtext(2, text="SVs Discovered", line=2)
  
  #Legend
  # legend("topright", pch=19, col=svtypes$color, 
  #        legend=svtypes$svtype[which(svtypes$svtype!="OTH")], 
  #        border=NA, bty="n", cex=0.8)
}

#Plot distribution of allele counts by SVTYPE
plot.freqByType <- function(dat, ymin=NULL, axlabel.cex=1){
  #Gather data
  exclude <- grep("MULTIALLELIC", dat$FILTER, fixed=T)
  if(length(exclude)>0){
    dat <- dat[-exclude, ]
  }
  dat <- dat[which(dat$SVTYPE!="MCNV" & !(dat$SV.chrom %in% c("X", "Y"))), ]
  AC.cutoffs <- c(1:9, seq(10, 90, 10), seq(100, 900, 100), seq(1000, 25000, 1000))
  AF.dat <- as.data.frame(sapply(c("DUP", "DEL", "INV"), function(svtype){
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
  sapply(c("DUP", "DEL", "INV"), function(svtype){
    points(x=log10(AC.cutoffs), y=AF.dat[, which(colnames(AF.dat)==svtype)], 
           col=svtypes$color[which(svtypes$svtype==svtype)], type="l", lwd=3)
  })
  sapply(c("DUP", "DEL", "INV"), function(svtype){
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

#Plot distribution of allele counts by coding context (for a single svtype)
plot.freqByContext <- function(dat, svtype, ymin=NULL, axlabel.cex=1,
                               parmar=c(2, 3.25, 0.5, 0.5)){
  #Gather data
  exclude <- grep("MULTIALLELIC", dat$FILTER, fixed=T)
  if(length(exclude)>0){
    dat <- dat[-exclude, ]
  }
  dat <- dat[which(dat$SVTYPE==svtype & !(dat$`SV.chrom` %in% c("X", "Y"))), ]
  AC.cutoffs <- c(1:9, seq(10, 90, 10), seq(100, 900, 100), seq(1000, 25000, 1000))
  #Helper function to compute cdf
  get.af.cdf <- function(dat, AC.cutoffs){
    sapply(AC.cutoffs, function(AC){
      length(which(dat$AC[which(dat$SVTYPE==svtype)]<=AC))/length(which(dat$SVTYPE==svtype))
    })
  }
  coding.AFs <- get.af.cdf(dat[which(!is.na(dat$LOF)
                                     | !is.na(dat$COPY_GAIN)
                                     | !is.na(dat$DUP_LOF)), ],
                           AC.cutoffs)
  coding.nonfunc.AFs <- get.af.cdf(dat[which(is.na(dat$LOF)
                                             & is.na(dat$COPY_GAIN)
                                             & is.na(dat$DUP_LOF)
                                             & (!is.na(dat$INTRONIC))
                                             | !is.na(dat$DUP_PARTIAL)), ],
                                   AC.cutoffs)
  noncoding.AFs <- get.af.cdf(dat[which(dat$INTERGENIC=="True"), ], 
                              AC.cutoffs)
  AF.dat <- data.frame("coding"=coding.AFs,
                       "coding.nonfunc"=coding.nonfunc.AFs,
                       "noncoding"=noncoding.AFs)
  #Plot
  if(is.null(ymin)){
    ymin <- log10(floor(100*min(AF.dat, na.rm=T))/100)
  }else{
    ymin <- log10(ymin)
  }
  AF.dat <- as.data.frame(apply(AF.dat, 2, log10))
  xrange <- c(-0.25, max(log10(AC.cutoffs)))
  common.threshold <- min(as.numeric(dat$AC[which(as.numeric(dat$AF)>=0.01 & as.numeric(dat$AN)==max(dat$AN, na.rm=T))]), na.rm=T)
  par(mar=parmar, bty="n")
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
  base.color <- svtypes$color[which(svtypes$svtype==svtype)]
  colors <- sapply(c(1, 2/3, 1/3), function(a){adjustcolor(base.color, alpha=a)})
  sapply(1:3, function(i){
    points(x=log10(AC.cutoffs), y=AF.dat[, i], 
           col=colors[i], type="l", lwd=3)
  })
  sapply(1:3, function(i){
    # points(x=log10(AC.cutoffs)[1], y=AF.dat[, which(colnames(AF.dat)==svtype)][1], 
    #        col="black", pch="-", cex=1.4)
    # points(x=log10(AC.cutoffs)[1], y=AF.dat[, which(colnames(AF.dat)==svtype)][1], 
    #        col=svtypes$color[which(svtypes$svtype==svtype)], pch="-", cex=1.2)
    rect(xleft=log10(AC.cutoffs)[1]-(0.03*(par("usr")[2]-par("usr")[1])), 
         xright=log10(AC.cutoffs)[1]+(0.03*(par("usr")[2]-par("usr")[1])), 
         ybottom=AF.dat[, i][1]-(0.01*(par("usr")[4]-par("usr")[3])), 
         ytop=AF.dat[, i][1]+(0.01*(par("usr")[4]-par("usr")[3])), 
         col="white")
    rect(xleft=log10(AC.cutoffs)[1]-(0.03*(par("usr")[2]-par("usr")[1])), 
         xright=log10(AC.cutoffs)[1]+(0.03*(par("usr")[2]-par("usr")[1])), 
         ybottom=AF.dat[, i][1]-(0.01*(par("usr")[4]-par("usr")[3])), 
         ytop=AF.dat[, i][1]+(0.01*(par("usr")[4]-par("usr")[3])), 
         col=colors[i])
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
  legend("right", legend=c("Genic (gene-altering)",
                           "Genic (no alteration)",
                           "Strictly intergenic"),
         cex=0.8*axlabel.cex, fill=colors, bg="white")
}

#Get proportion of singletons and 95% CI for a set of variants
calc.fracSingletons.singleClass <- function(dat, boot.n=100, conf=0.95, aps=F, min.sv=10){
  dat <- dat[which(!(dat$chrom %in% c("X", "Y"))), ]
  ACs <- dat$AC
  v.aps <- dat$APS
  ACs <- ACs[which(!is.na(ACs) & ACs>0)]
  v.aps <- v.aps[which(!is.na(v.aps))]
  if(length(ACs) >= min.sv & length(v.aps) >= min.sv){
    if(aps==T){
      helper.getFracSingle <- function(v.aps, indices){sum(v.aps[indices])/length(v.aps[indices])}
      point.est <- helper.getFracSingle(v.aps, indices=1:length(v.aps))
      calc.ci <- function(v.aps, n, conf){
        set.seed(0)
        boot.obj <- boot(data=v.aps, statistic=helper.getFracSingle, R=n)
        ci <- boot.ci(boot.obj, conf=conf, type="basic")$basic[4:5]
        return(ci)
      }
      ci <- calc.ci(v.aps, n=boot.n, conf=conf)
    }else{
      helper.getFracSingle <- function(ACs, indices){length(which(ACs[indices]==1))/length(ACs[indices])}
      point.est <- helper.getFracSingle(ACs, indices=1:length(ACs))
      calc.ci <- function(ACs, n, conf){
        set.seed(0)
        boot.obj <- boot(data=ACs, statistic=helper.getFracSingle, R=n)
        ci <- boot.ci(boot.obj, conf=conf, type="basic")$basic[4:5]
        return(ci)
      }
      ci <- calc.ci(ACs, n=boot.n, conf=conf)
    }
    if(length(ci) != 2){
      ci <- c(NA, NA)
    }
    return(c(point.est, ci))
  }else{
    return(c(NA, NA, NA))
  }
}

#Gather fraction of singletons for a single SV class in bins by size
calc.fracSingletons.singleClass.binned <- function(dat, svtype, 
                                                   include.ins=F, 
                                                   d=c(0, 1000, 10000, 100000, 1000000000),
                                                   aps=F){
  if(svtype!="ALL"){
    tmpdat <- dat[which(dat$SVTYPE==svtype & dat$FILTER=="PASS"), ]
  }else{
    tmpdat <- dat[which(dat$FILTER=="PASS"), ]
  }
  if(include.ins==F){
    tmpdat <- tmpdat[which(tmpdat$SVTYPE!="INS"), ]
  }
  if(is.null(d)){
    d <- quantile(tmpdat$SVLEN, probs=seq(0, 1, 0.1))
  }
  d.singles <- t(sapply(2:length(d), function(i){
    boot.dat <- tmpdat[which(tmpdat$SVLEN>=d[i-1] & tmpdat$SVLEN<d[i]), ]
    if(nrow(boot.dat) > 0){
      calc.fracSingletons.singleClass(dat=boot.dat, aps=aps)
    }else{
      return(c(NA, NA, NA))
    }
  }))
  return(d.singles)
}

#Plot fraction of singletons by size for a single SV class
plot.fracSingletons.bySize <- function(dat, svtype, 
                                       include.ins=T, 
                                       d=c(0, 1000, 10000, 100000, 1000000, 1000000000), 
                                       d.labels=c("<1kb", "1-10kb", "10-100kb", "100kb-1Mb", ">1Mb"), 
                                       color="black", title=NULL, mar=c(3, 2.75, 0.5, 0.5), ylab.cex=1.2, 
                                       ylims=NULL, aps=F, dual=F){
  #Get plot data
  p.dat <- calc.fracSingletons.singleClass.binned(dat, svtype, 
                                                  include.ins, d, 
                                                  aps=aps)
  if(dual==T){
    p.dat.intergenic <- calc.fracSingletons.singleClass.binned(dat[which(dat$PROTEIN_CODING__INTERGENIC=="True" 
                                                                         | (is.na(dat$PROTEIN_CODING__LOF)
                                                                            & is.na(dat$PROTEIN_CODING__DUP_LOF)
                                                                            & is.na(dat$PROTEIN_CODING__COPY_GAIN)
                                                                            & is.na(dat$PROTEIN_CODING__DUP_PARTIAL)
                                                                            & is.na(dat$PROTEIN_CODING__INTRONIC)
                                                                            & is.na(dat$PROTEIN_CODING__PROMOTER)
                                                                            & is.na(dat$PROTEIN_CODING__UTR))), ],
                                                               sd_sr_cov, svtype, 
                                                               max.cov, include.ins, d, 
                                                               aps=aps)
  }
  if(aps==F){
    exclude <- which(apply(p.dat, 1, function(vals){any(vals==0)}))
    if(length(exclude)>0){
      p.dat[exclude, 1:3] <- NA
    }
    p.dat[which(p.dat>1)] <- 1
    p.dat[which(p.dat<0)] <- 0
    if(dual==T){
      exclude.int <- which(apply(p.dat.intergenic, 1, function(vals){any(vals==0)}))
      if(length(exclude.int)>0){
        p.dat.intergenic[exclude.int, 1:3] <- NA
      }
      p.dat.intergenic[which(p.dat.intergenic>1)] <- 1
      p.dat.intergenic[which(p.dat.intergenic<0)] <- 0
      if(is.null(ylims)){
        ylims <- range(rbind(p.dat, p.dat.intergenic), na.rm=T)
      }
    }else{
      if(is.null(ylims)){
        ylims <- range(p.dat, na.rm=T)
      }
    }
  }else{
    if(dual==T){
      if(is.null(ylims)){
        if(max(rbind(p.dat, p.dat.intergenic), na.rm=T) < 0.2){
          ylims <- c(-2, 2)*max(rbind(p.dat, p.dat.intergenic), na.rm=T)
        }else{
          ylims <- c(-1, 1)*max(rbind(p.dat, p.dat.intergenic), na.rm=T)
        }
      }
    }else{
      if(is.null(ylims)){
        if(max(p.dat, na.rm=T) < 0.2){
          ylims <- c(-2, 2)*max(p.dat, na.rm=T)
        }else{
          ylims <- c(-1, 1)*max(p.dat, na.rm=T)
        }
      }
    }
  }
  #Prep plot area
  par(mar=mar, bty="n")
  plot(x=c(0.25, length(d)-1.25), y=ylims, type="n", 
       xlab="", xaxt="n", ylab="", yaxt="n")
  if(aps==T){
    abline(h=0, lty=2, col="gray70")
  }
  #Plot dots and CIs
  if(dual==T){
    segments(x0=(1:(length(d)-1))-0.65, 
             x1=(1:(length(d)-1))-0.65, 
             y0=p.dat.intergenic[, 2], y1=p.dat.intergenic[, 3], 
             lend="round", lwd=1.75, col=adjustcolor(color, alpha=0.5))
    points(x=(2:length(d))-1.65, y=p.dat.intergenic[, 1], pch=19, col="white", lwd=0)
    points(x=(2:length(d))-1.65, y=p.dat.intergenic[, 1], pch=19, col=adjustcolor(color, alpha=0.5), lwd=0)
    segments(x0=(1:(length(d)-1))-0.35, 
             x1=(1:(length(d)-1))-0.35, 
             y0=p.dat[, 2], y1=p.dat[, 3], 
             lend="round", lwd=1.75, col=color)
    points(x=(2:length(d))-1.35, y=p.dat[, 1], pch=19, col=color)
  }else{
    segments(x0=(1:(length(d)-1))-0.5, 
             x1=(1:(length(d)-1))-0.5, 
             y0=p.dat[, 2], y1=p.dat[, 3], 
             lend="round", lwd=2, col=color)
    points(x=(2:length(d))-1.5, y=p.dat[, 1], pch=19, col=color)
  }
  #Add axes
  text(x=(1:(length(d)-1))+0.1, y=par("usr")[3]-(0.025*(par("usr")[4]-par("usr")[3])), 
       xpd=T, srt=45, labels=d.labels, cex=0.85, pos=2)
  if(par("usr")[4]-par("usr")[3] > 1.5){
    y.at <- seq(-1, 1, 0.2)
  }else{
    y.at <- seq(-1, 1, 0.1)
  }
  axis(2, at=y.at, labels=NA, tck=-0.03)
  if(aps==T){
    axis(2, at=y.at, tick=F, las=2, cex.axis=0.7, line=-0.7)
    mtext(2, line=1.7, text="APS", cex=ylab.cex)
  }else{
    axis(2, at=y.at, tick=F, labels=paste(round(100*y.at, 0), "%", sep=""), 
         las=2, cex.axis=0.7, line=-0.7)
    mtext(2, line=1.7, text="Singleton Proportion", cex=ylab.cex)
  }
  mtext(3, line=0, text=paste("     ", title, sep=""), cex=1.2*ylab.cex)
}


# #Plot distribution of allele counts by SVLEN decile
# plot.freqBySize <- function(dat){
#   #Gather data
#   AC.cutoffs <- c(1:9, seq(10, 90, 10), seq(100, 900, 100), seq(1000, 25000, 1000))
#   AF.dat.sub <- dat[which(dat$SVTYPE %in% c("DEL", "DUP", "INV", "CPX")), ]
#   size.cutoffs <- quantile(AF.dat.sub$SVLEN, probs=(0:5)/5)
# 
#   # AF.dat.sub <- dat[which(!(dat$SVTYPE %in% c("BND", "CTX"))), ]
#   # AF.dat.sub <- AF.dat.sub[-grep("INS:ME", AF.dat.sub$svtype, fixed=T), ]
#   AF.dat <- as.data.frame(sapply(1:5, function(i){
#     sapply(AC.cutoffs, function(AC){
#       length(which(AF.dat.sub$AC[which(AF.dat.sub$SVLEN>size.cutoffs[i] & AF.dat.sub$SVLEN<=size.cutoffs[i+1])]<=AC))/length(which(AF.dat.sub$SVLEN>size.cutoffs[i] & AF.dat.sub$SVLEN<=size.cutoffs[i+1]))
#     })
#   }))
#   #Plot
#   ymin <- log10(floor(100*min(AF.dat, na.rm=T))/100)
#   AF.dat <- as.data.frame(apply(AF.dat, 2, log10))
#   xrange <- c(-0.25, max(log10(AC.cutoffs)))
#   common.threshold <- min(dat$AC[which(dat$AF>=0.01 & dat$AN==max(dat$AN))])
#   col.pal <- colorRampPalette(c("#440154", "#365C8C", "#25A584", "#FDE725"))(5)
#   par(mar=c(2.5, 3.25, 0.5, 0.5), bty="n")
#   plot(x=xrange, y=c(ymin, log10(1)), type="n", 
#        xaxt="n", xlab="", yaxt="n", ylab="")
#   # rect(xleft=-0.1, xright=0.1, ybottom=par("usr")[3], ytop=par("usr")[4], 
#   #      border=NA, bty="n", col="gray90")
#   # abline(v=log10(common.threshold), col="gray80")
#   segments(x0=log10(common.threshold), x1=log10(common.threshold), 
#            y0=par("usr")[3], y1=log10(1), col="gray80")
#   text(x=log10(common.threshold), y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])), 
#        labels="Rare\n(AF<1%)", pos=2, cex=0.7)
#   text(x=log10(common.threshold), y=par("usr")[3]+(0.08*(par("usr")[4]-par("usr")[3])), 
#        labels="Common\n(AF>1%)", pos=4, cex=0.7)
#   sapply(1:5, function(i){
#     points(x=log10(AC.cutoffs), y=AF.dat[, i], 
#            col=col.pal[i], type="l", lwd=3)
#   })
#   sapply(c("OTH", "INS", "DUP", "DEL", "INV", "CPX"), function(svtype){
#     points(x=log10(AC.cutoffs)[1], y=AF.dat[, which(colnames(AF.dat)==svtype)][1], 
#            bg=svtypes$color[which(svtypes$svtype==svtype)], pch=21, cex=1.2)
#   })
#   logscale.all <- log10(as.numeric(unlist(sapply(c(0:9), function(i){(1:9)*(10^i)}))))
#   logscale.major <- 0:9
#   axis(1, at=logscale.all, labels=NA, tck=-0.015, lwd=0.7)
#   axis(1, at=logscale.major, labels=NA, tck=-0.03, lwd=1.1)
#   axis(1, at=c(0:4), tick=F, line=-0.7, cex.axis=0.8, 
#        labels=c("1", "10", "100", "1k", "10k"))
#   mtext(1, text="Maximum Allele Count", line=1.3)
#   logscale.pct.all <- log10((1:100)/100)
#   logscale.pct.major <- log10(seq(10, 100, 10)/100)
#   axis(2, at=logscale.pct.all, labels=NA, tck=-0.015, lwd=0.9)
#   axis(2, at=logscale.pct.major, labels=NA, tck=-0.03, lwd=1.1)
#   axis(2, at=logscale.pct.major, tick=F, line=-0.5, cex.axis=0.8, las=2, 
#        labels=paste(seq(10, 100, 10), "%", sep=""))
#   mtext(2, text="Fraction of SVs", line=2.25)
# }

#Gather Hardy-Weinberg data
makeHWEmat <- function(dat, pop=NULL){
  sub.dat <- dat[which(!(dat$`SV.chrom` %in% c("X", "Y"))), ]
  cols.to.exclude <- grep("MULTIALLELIC", sub.dat$FILTER, fixed=T)
  if(length(cols.to.exclude)>0){
    sub.dat <- sub.dat[-cols.to.exclude, ]
  }
  # sub.dat <- sub.dat[-grep("PESR_GT_OVERDISPERSION", sub.dat$FILTER, fixed=T), ]
  # sub.dat <- sub.dat[-grep("UNRESOLVED", sub.dat$FILTER, fixed=T), ]
  if(!is.null(pop)){
    n.genos.idx <- which(colnames(sub.dat)==paste(pop, "_N_BI_GENOS", sep=""))
    n.homref.idx <- which(colnames(sub.dat)==paste(pop, "_N_HOMREF", sep=""))
    n.het.idx <- which(colnames(sub.dat)==paste(pop, "_N_HET", sep=""))
    n.homalt.idx <- which(colnames(sub.dat)==paste(pop, "_N_HOMALT", sep=""))
  }else{
    n.genos.idx <- which(colnames(sub.dat)=="N_BI_GENOS")
    n.homref.idx <- which(colnames(sub.dat)=="N_HOM_REF")
    n.het.idx <- which(colnames(sub.dat)=="N_HET")
    n.homalt.idx <- which(colnames(sub.dat)=="N_HOMALT")
  }
  sub.dat <- sub.dat[which(sub.dat[, n.genos.idx]>0), ]
  sub.dat <- sub.dat[which(apply(sub.dat[, c(n.het.idx, n.homalt.idx)], 1, sum)>0), ]
  HWE.mat <- data.frame("AA"=sub.dat[, n.homref.idx], 
                        "AB"=sub.dat[, n.het.idx], 
                        "BB"=sub.dat[, n.homalt.idx])
  HWE.mat <- HWE.mat[complete.cases(HWE.mat), ]
  return(HWE.mat)
}
#Hardy-Weinberg ternary plot
plot.HWE <- function(dat, pop=NULL, title=NULL, full.legend=F, lab.cex=1){
  require(HardyWeinberg, quietly=T)
  #Gather HW p-values & colors
  HWE.mat <- makeHWEmat(dat=dat, pop=pop)
  HW.p <- HWChisqStats(X=HWE.mat, x.linked=F, pvalues=T)
  HW.cols <- rep("#4DAC26", times=length(HW.p))
  HW.cols[which(HW.p<0.05)] <- "#81F850"
  HW.cols[which(HW.p<0.05/length(HW.p))] <- "#AC26A1"
  
  #Generate HW plot frame
  par(mar=c(1, 1, 1, 1), bty="n")
  plot(x=1.15*c(-1/sqrt(3), 1/sqrt(3)), y=c(-0.15, 1.15), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  segments(x0=c(-1/sqrt(3), 0, 1/sqrt(3)), 
           x1=c(0, 1/sqrt(3), -1/sqrt(3)), 
           y0=c(0, 1, 0), y1=c(1, 0, 0))
  HWTernaryPlot(X=HWE.mat, n=max(HWE.mat, na.rm=T), newframe=F, 
                vbounds=F, mafbounds=F, 
                region=1, vertexlab=NA, 
                alpha=0.05, 
                curvecols=c("#4DAC26", "#81F850", NA, NA), pch=NA)
  
  #Add axes
  text(x=c(-1/sqrt(3), 1/sqrt(3)), y=0, labels=c("0/0", "1/1"), 
       pos=1, cex=0.8, xpd=T, font=2)
  text(x=0, y=1, labels="0/1", pos=3, cex=0.8, xpd=T, font=2)
  
  #Finish HW plot
  HWTernaryPlot(X=HWE.mat, n=max(HWE.mat, na.rm=T), newframe=F, 
                vbounds=F, mafbounds=F, 
                region=1, vertexlab=NA, 
                alpha=0.03/nrow(HWE.mat), 
                curvecols=c("#4DAC26", "#AC26A1", NA, NA), 
                pch=21, cex=0.3, signifcolour=F, markercol=NA, 
                markerbgcol=adjustcolor(HW.cols, alpha=0.25))
  segments(x0=c(-1/sqrt(3), 0, 1/sqrt(3)), 
           x1=c(0, 1/sqrt(3), -1/sqrt(3)), 
           y0=c(0, 1, 0), y1=c(1, 0, 0))
  
  #Add legend
  n.pass <- length(which(HW.p>=0.05))
  print(paste("PASS: ", n.pass/length(HW.p), sep=""))
  n.nom <- length(which(HW.p<0.05 & HW.p>=0.05/nrow(HWE.mat)))
  print(paste("NOMINAL FAILS: ", n.nom/length(HW.p), sep=""))
  n.bonf <- length(which(HW.p<0.05/nrow(HWE.mat)))
  print(paste("BONFERRONI FAILS: ", n.bonf/length(HW.p), sep=""))
  legend("right", pch=19, col=c("#4DAC26", "#81F850", "#AC26A1"), pt.cex=1.3, 
         legend=c(paste(round(100*(n.pass/nrow(HWE.mat)), 0), "%", sep=""), 
                  paste(round(100*(n.nom/nrow(HWE.mat)), 0), "%", sep=""), 
                  paste(round(100*(n.bonf/nrow(HWE.mat)), 0), "%", sep="")), 
         bty="n", bg=NA, cex=0.7)
  text(x=par("usr")[2], y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])), pos=2, cex=0.7, 
       labels=paste(title, "\n \n ", sep=""), font=2)
  text(x=par("usr")[2], y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])), pos=2, cex=0.7, 
       labels=paste(" \n", prettyNum(max(apply(HWE.mat, 1, sum), na.rm=T), big.mark=","), " Samples\n ", sep=""))
  text(x=par("usr")[2], y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])), pos=2, cex=0.7, 
       labels=paste(" \n \n", prettyNum(nrow(HWE.mat), big.mark=","), " SVs", sep=""))
}

#Allele frequency correlation plot between populations
plot.crossPopCorr <- function(dat, pops, popA, popB){
  #Subset to biallelic,  autosomal sites with non-zero AF in both populations
  AN.A.idx <- which(colnames(dat)==paste(popA, "_AN", sep=""))
  AC.A.idx <- which(colnames(dat)==paste(popA, "_AC", sep=""))
  AF.A.idx <- which(colnames(dat)==paste(popA, "_AF", sep=""))
  AN.B.idx <- which(colnames(dat)==paste(popB, "_AN", sep=""))
  AC.B.idx <- which(colnames(dat)==paste(popB, "_AC", sep=""))
  AF.B.idx <- which(colnames(dat)==paste(popB, "_AF", sep=""))
  subdat <- dat[which(dat[, AF.A.idx]>0 & dat[, AF.B.idx]>0 & !(dat$SV.chrom %in% c("X", "Y"))), ]
  rows.to.drop <- sort(unique(c(grep("MULTIALLELIC", subdat$FILTER, fixed=T), 
                                grep("PESR_GT_OVERDISPERSION", subdat$FILTER, fixed=T), 
                                grep("UNRESOLVED", subdat$FILTER, fixed=T))))
  
  if(length(rows.to.drop)>0){
    subdat <- subdat[-rows.to.drop, ]
  }
  
  #Calculate p-values with chi-square test
  pvals <- sapply(1:nrow(subdat), function(i){
    AN.A <- subdat[i, AN.A.idx]
    AC.A <- subdat[i, AC.A.idx]
    AN.B <- subdat[i, AN.B.idx]
    AC.B <- subdat[i, AC.B.idx]
    return(chisq.test(data.frame(c(AN.A-AC.A, AC.A), c(AN.B-AC.B, AC.B)))$p.value)
  })
  pvals.bonf <- p.adjust(pvals, method="bonferroni")
  
  #Prepare plot
  logscale.all <- log10(as.numeric(unlist(sapply(c(0:9), function(i){(1:9)*(10^i)}))))
  logscale.major <- 0:9
  major.labels <- sapply(logscale.major, function(i){expression(paste(i^"th"))})
  par(mar=c(2.6, 2.6, 1.5, 1.5))
  plot(x=log10(c(min(subdat[, AF.A.idx]), 1)), 
       y=log10(c(min(subdat[, AF.B.idx]), 1)), 
       type="n", xaxt="n", yaxt="n", xlab="", ylab="")
  axis(1, at=-logscale.all, labels=NA, tck=-0.015, lwd=0.7)
  axis(1, at=-logscale.major, labels=NA, tck=-0.03, lwd=1.1)
  mtext(1, text=paste(popA, " AF", sep=""), line=1.25)
  axis(2, at=-logscale.all, labels=NA, tck=-0.015, lwd=0.7)
  axis(2, at=-logscale.major, labels=NA, tck=-0.03, lwd=1.1)
  mtext(2, text=paste(popB, " AF", sep=""), line=1.6)
  sapply(-logscale.major, function(i){
    axis(1, at=i, labels=bquote('10'^.(i)), tick=F, line=-0.7, cex.axis=0.8)
    axis(2, at=i, labels=bquote('10'^.(i)), tick=F, line=-0.6, cex.axis=0.8, las=2)
  })
  
  
  
  
  #Add points
  pt.cex <- 0.3
  alpha <- 0.1
  col.A <- pops$color[which(pops$pop==popA)]
  col.B <- pops$color[which(pops$pop==popB)]
  points(x=as.numeric(unlist(log10(subdat[which(pvals.bonf<=0.05 & subdat[, AF.A.idx]>subdat[, AF.B.idx]), AF.A.idx]))), 
         y=as.numeric(unlist(log10(subdat[which(pvals.bonf<=0.05 & subdat[, AF.A.idx]>subdat[, AF.B.idx]), AF.B.idx]))), 
         pch=19, cex=pt.cex, lwd=0, 
         col=adjustcolor(col.A, alpha=alpha))
  points(x=as.numeric(unlist(log10(subdat[which(pvals.bonf<=0.05 & subdat[, AF.B.idx]>subdat[, AF.A.idx]), AF.A.idx]))), 
         y=as.numeric(unlist(log10(subdat[which(pvals.bonf<=0.05 & subdat[, AF.B.idx]>subdat[, AF.A.idx]), AF.B.idx]))), 
         pch=19, cex=pt.cex, lwd=0, 
         col=adjustcolor(col.B, alpha=alpha))
  points(x=as.numeric(unlist(log10(subdat[which(pvals.bonf>0.05), AF.A.idx]))), 
         y=as.numeric(unlist(log10(subdat[which(pvals.bonf>0.05), AF.B.idx]))), 
         pch=19, cex=pt.cex, lwd=0, 
         col=adjustcolor("gray50", alpha=alpha))
  
  #Add stats
  # abline(lm(subdat[, AF.B.idx] ~ subdat[, AF.A.idx]), lty=2, col="gray50")
  frac.A <- length(which(pvals.bonf<=0.05 & subdat[, AF.A.idx]>subdat[, AF.B.idx]))/nrow(subdat)
  text(x=par("usr")[2]+(0.03*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[3]+(0.05*(par("usr")[4]-par("usr")[3])), 
       pos=2, cex=0.8, col=col.A, 
       labels=paste(round(100*frac.A, digits=1), "%", sep=""))
  frac.B <- length(which(pvals.bonf<=0.05 & subdat[, AF.B.idx]>subdat[, AF.A.idx]))/nrow(subdat)
  text(x=par("usr")[1]-(0.03*(par("usr")[2]-par("usr")[1])), 
       y=par("usr")[4]-(0.05*(par("usr")[4]-par("usr")[3])), 
       pos=4, cex=0.8, col=col.B, 
       labels=paste(round(100*frac.B, digits=1), "%", sep=""))
  AB.cor <- format(round(cor(subdat[, AF.A.idx], subdat[, AF.B.idx])^2, 3), nsmall=3)
  mtext(3, line=0.1, text=bquote(italic(R)^2 == .(AB.cor)))
}

###Perform Watterson estimator analysis of mutation rate
#Individual function to get mutation rate estimate for a single pop & svtype
getMu <- function(dat, pop=NULL, svtype=NULL, Ne=10000){
  #Format variables
  if(!is.null(pop)){
    pop <- paste(pop, "_", sep="")
  }
  #Get number of chromosomes assessed
  pop.AN.idx <- which(colnames(dat)==paste(pop, "N_BI_GENOS", sep=""))
  n <- 2*max(dat[, pop.AN.idx], na.rm=T)
  #Get number of autosomal biallelic segregating sites
  pop.AF.idx <- which(colnames(dat)==paste(pop, "AF", sep=""))
  if(!is.null(svtype)){
    K <- length(which(dat[, pop.AF.idx]>0 & dat$SVTYPE==svtype & !(dat$chrom %in% c("X", "Y")) & dat$SVTYPE != "MCNV"))
  }else{
    K <- length(which(dat[, pop.AF.idx]>0 & !(dat$chrom %in% c("X", "Y")) & dat$SVTYPE != "MCNV"))
  }
  #Get n-1th harmonic number
  harmsum <- sum(sapply(1:(n-1), function(k){1/k}))
  #Get Watterson estimator
  theta.hat <- K/harmsum
  #Solve for mutation rate
  mu <- theta.hat/(4*Ne)
  return(mu)
}
#Wrapper function to get mutation rate estimates for all svtypes & populations
getAllMus <- function(dat, Ne){
  mu.svtypes <- c("DEL", "DUP", "INV")
  #Get mutation rates across classes
  mu.AllClasses <- sapply(mu.svtypes, function(svtype){
    getMu(dat, pop=NULL, svtype=svtype, Ne=Ne)
  })
  #Get mutation rate of all SV across populations
  mu.AllPops <- sapply(pops$pop, function(pop){
    getMu(dat, pop=pop, svtype=NULL, Ne=Ne)
  })
  #Get mutation rate by svtype & population
  mu.ClassByPop <- sapply(pops$pop, function(pop){
    sapply(mu.svtypes, function(svtype){
      getMu(dat, pop=pop, svtype=svtype, Ne=Ne)
    })
  })
  #Get mean by class across pops
  mu.PopMeanByClass <- apply(mu.ClassByPop, 1, function(vals){
    lower <- as.numeric(t.test(vals)$conf.int[1])
    mean <- as.numeric(t.test(vals)$estimate)
    upper <- as.numeric(t.test(vals)$conf.int[2])
    return(c(lower, mean, upper))
  })
  mu.PopMeanByClass <- cbind(apply(mu.PopMeanByClass, 1, sum), mu.PopMeanByClass)
  colnames(mu.PopMeanByClass)[1] <- "ALL"
  rownames(mu.PopMeanByClass) <- c("CI.lowerBound", "mean", "CI.upperBound")
  return(list("muByClass"=mu.AllClasses, 
              "muByPop"=mu.AllPops, 
              "muByClassAndPop"=mu.ClassByPop, 
              "mu.Means"=mu.PopMeanByClass))
}
#Plot mutation rates
plotMus <- function(dat, Ne, werling=T){
  mu.dat <- getAllMus(dat, Ne)
  plot.dat <- mu.dat$mu.Means
  mu.svtypes <- c("DEL", "DUP", "INV")
  mu.cols <- c("gray30", sapply(mu.svtypes, function(svtype){
    svtypes$color[which(svtypes$svtype==svtype)]
  }))
  werling.rates <- c(166, 80, 47, 37, 0, 1)/1038
  #Prep plot area
  par(mar=c(1.5, 3.5, 0.5, 0.5), bty="n")
  plot(x=c(-0.025, ncol(plot.dat))+0.2, y=c(0, max(plot.dat)), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  axis(1, at=c(-4, 10), tck=0, labels=NA, lwd=1.5)
  axis(2, at=c(-1, 1), tck=0, labels=NA, lwd=1.5)
  axis(1, at=(1:ncol(plot.dat))-0.5, tick=F, line=-0.8, 
       labels=c("All SV", mu.svtypes))
  axis(2, labels=NA)
  axis(2, tick=F, las=2, cex.axis=0.8, line=-0.4)
  mtext(2, line=2, text="Mutation Rate (SVs per generation)")
  #Add points & CIs
  text.buf <- 0.04*(par("usr")[4]-par("usr")[3])
  sapply(1:ncol(plot.dat), function(i){
    segments(x0=i-0.5, x1=i-0.5, 
             y0=plot.dat[1, i], y1=plot.dat[3, i], 
             lend="round", lwd=2, col=mu.cols[i])
    points(x=i-0.5, y=plot.dat[2, i], pch=19, cex=1.5, col=mu.cols[i])
    par(xpd=T)
    # text(x=i-c(0.6, 0.45, 0.6), y=plot.dat[, i]+c(-text.buf, 0, text.buf), 
    #      cex=c(0.6, 0.7, 0.6), labels=format(round(plot.dat[, i], 4), nsmall=4), pos=4, 
    #      col=mu.cols[i], font=c(3, 2, 3))
    text(x=i-0.45, y=plot.dat[2, i], 
         cex=1, labels=format(round(plot.dat[2, i], 3), nsmall=3), pos=4, 
         col=mu.cols[i], font=2)
    par(xpd=F)
    if(werling==T){
      points(x=i-0.5, y=werling.rates[i], pch=23, lwd=1.5)
    }
  })
  if(werling==T){
    legend("topright", border=NA, bty="n", 
           legend=c(as.expression(bquote(mu ~ "from Watterson" ~ hat(theta[italic(W)]) ~ "in gnomAD,  " ~ N[e] == .(prettyNum(Ne, big.mark=",")))), 
                    expression("Rate of validated" ~ italic("de novo") ~ "SVs from 519 quartets")), 
           pch=c(19, 23), pt.lwd=1.5, lwd=c(2, NA), lty=c(1, NA), cex=0.8, pt.cex=c(1.25, 1))
  }else{
    legend("topright", border=NA, bty="n", 
           legend=as.expression(bquote(mu ~ "from Watterson" ~ hat(theta[italic(W)]) ~ "in gnomAD,  " ~ N[e] == .(prettyNum(Ne, big.mark=",")))), 
           pch=19, pt.lwd=1.5, lwd=2, lty=1, cex=0.8, pt.cex=1.25)
  }
}


###################
###COMPLEX SV TABLE
###################
#Master function to make complex SV table
masterCPXtable <- function(dat, 
                           cpx.types=c("delINV", "INVdel", "delINVdel", 
                                       "dupINV", "INVdup", "dupINVdup", 
                                       "delINVdup", "dupINVdel", 
                                       # "piDUP_FR", "piDUP_RF", 
                                       "dDUP", "dDUP_iDEL", "INS_iDEL")){
  #Set parameters
  filler.color <- "gray50"
  panels <- c("CPX Class", "Involved SV\nSignatures", "Rearrangement Schematic", "Variants", 
              "SV per\nGenome", "Median\nSize", "Size\nDistribution", "Mean AF")
  n.panels <- length(panels)
  panel.widths=c(0.75, 0.55, 2, 0.5, 0.5, 0.5, 1, 0.5)
  panel.heights=c(1, rep(1, times=length(cpx.types)+1), 0.75)
  
  #Empty header text box with single value
  headerTextBox <- function(text, cex=1){
    plot(x=c(0, 1), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    text(x=0.5, y=par("usr")[3], pos=3, labels=text, cex=cex, font=2)
    axis(1, at=c(par("usr")[1], par("usr")[2]), tck=0, labels=NA)
  }
  
  #Empty text box with single value
  textBox <- function(text, cex=1, pos=2){
    plot(x=c(0, 1), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    text(x=par("usr")[2], y=0.5, pos=pos, labels=text, cex=cex)
  }
  
  #Blank box with nothing
  blankBox <- function(){
    plot(x=c(0, 1), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
  }
  
  #Boxes for SV types involved in CPX class
  cpx.signature <- function(cpx.type){
    plot(x=c(0, 4), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    rect(xleft=(0:3)+0.2, xright=(0:3)+0.8, 
         ybottom=0.2, ytop=0.8, col="gray95", border="gray90")
    if(length(grep("DEL", cpx.type, ignore.case=T))>0){
      rect(xleft=0.2, xright=0.8, ybottom=0.2, ytop=0.8, 
           col=svtypes$color[which(svtypes$svtype=="DEL")])
    }
    if(length(grep("DUP", cpx.type, ignore.case=T))>0){
      rect(xleft=1.2, xright=1.8, ybottom=0.2, ytop=0.8, 
           col=svtypes$color[which(svtypes$svtype=="DUP")])
    }
    if(length(c(grep("INS", cpx.type, ignore.case=T), 
                grep("dDUP", cpx.type, ignore.case=T)))>0){
      rect(xleft=2.2, xright=2.8, ybottom=0.2, ytop=0.8, 
           col=svtypes$color[which(svtypes$svtype=="INS")])
    }
    if(length(c(grep("INV", cpx.type, ignore.case=T), 
                grep("piDUP", cpx.type, ignore.case=T)))>0){
      rect(xleft=3.2, xright=3.8, ybottom=0.2, ytop=0.8, 
           col=svtypes$color[which(svtypes$svtype=="INV")])
    }
  }
  
  #Simple barplot
  cpx.simpleBar <- function(value, max.value=1){
    plot(x=c(0, max.value), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    rect(xleft=0, xright=value, 
         ybottom=0.25, ytop=0.75, 
         col=filler.color)
  }
  
  #Simple size density plot
  cpx.SVLEN.dens <- function(dat, cpx.type){
    require(zoo, quietly=T)
    if(cpx.type=="ALL"){
      sizes <- log10(dat$SVLEN[which(dat$SVTYPE=="CPX")])
    }else{
      sizes <- log10(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)])
    }
    h <- hist(sizes, breaks=seq(0, 10, 0.05), plot=F)$counts
    h <- rollapply(h, 11, mean, partial=T)
    h <- h/max(h, na.rm=T)
    out.df <- data.frame("min.size"=10^seq(0, 9.95, 0.05), 
                         "dens"=h)
    plot(x=log10(c(50, 10000000)), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    polygon(x=log10(c(out.df$min.size, rev(out.df$min.size))), 
            y=c(out.df$dens, rep(0, times=nrow(out.df))), 
            col=filler.color)
    axis(1, at=c(par("usr")[1], par("usr")[2]), tck=0, labels=NA)
    abline(v=median(sizes), col=svtypes$color[which(svtypes$svtype=="CPX")], lwd=2)
  }
  
  #Pct axis
  pctAxis <- function(max.value=1){
    plot(x=c(0, max.value), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    axis(3, at=axTicks(3), labels=NA, tck=0.125)
    axis(3, at=axTicks(3), tick=F, line=-2.25, cex.axis=0.7, 
         labels=paste(round(100*axTicks(3), 1), "%", sep=""))
    # mtext(3, text="Pct.", line=-2.5, cex=0.7)
  }
  
  #Size axis
  sizeAxis <- function(){
    plot(x=log10(c(50, 10000000)), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    axis(3, at=as.numeric(unlist(sapply(1:7, function(x){log10((1:9)*(10^x))}))), labels=NA, tck=0.1, lwd=0.7)
    axis(3, at=2:7, labels=NA, tck=0.2)
    labels <- c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb")
    sapply(2:7, function(i){
      axis(3, at=i, tick=F, labels=labels[i], line=-2.25, cex.axis=0.7)
    })
    # mtext(3, text="SV Size", line=-2.5, cex=0.7)
  }
  
  #Prep figure layout
  layout(matrix(1:(n.panels*(length(cpx.types)+3)), byrow=T, nrow=length(cpx.types)+3), 
         widths=panel.widths, heights=panel.heights)
  
  #Headers
  par(bty="n", mar=c(0.1, 1, 0.1, 1), xpd=T)
  sapply(panels, function(title){
    headerTextBox(title)
  })
  
  #Add top row for summary of all CPX
  par(bty="n", mar=c(0.5, 0.5, 0.5, 0.5), xpd=F)
  textBox(text="All CPX SV")
  blankBox()
  blankBox()
  textBox(text=prettyNum(length(which(dat$SVTYPE=="CPX")), big.mark=","))
  textBox(text=prettyNum(round(sum(c(dat$N_HET[which(dat$SVTYPE=="CPX")], 
                                     dat$N_HOMALT[which(dat$SVTYPE=="CPX")]))/max(dat$N_BI_GENOS[which(dat$SVTYPE=="CPX")]), 
                               1), big.mark=","))
  textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE=="CPX")])/1000, 1), nsmall=1), "kb", sep=""))
  cpx.SVLEN.dens(dat, cpx.type="ALL")
  textBox(text=paste(format(round(100*mean(dat$AF[which(dat$SVTYPE=="CPX")], na.rm=T), 2), nsmall=2), "%", sep=""))
  
  #Iterate over cpx types - one row per type
  par(bty="n", mar=c(0.5, 0.5, 0.5, 0.5), xpd=F)
  sapply(1:length(cpx.types), function(i){
    cpx.type <- cpx.types[i]
    #CPX class (text)
    textBox(text=cpx.type)
    #Signatures
    cpx.signature(cpx.type)
    #Placeholder for rearrangement schematic (done in illustrator)
    blankBox()
    #Count of variants
    textBox(text=prettyNum(length(which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)), big.mark=","))
    # #Pct of all CPX
    # cpx.simpleBar(value=length(which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type))/length(which(dat$SVTYPE=="CPX")), 
    #               max.value=max(table(dat$CPX_TYPE[which(dat$SVTYPE=="CPX")])/length(which(dat$SVTYPE=="CPX"))))
    #Mean # of variants per genome
    textBox(text=prettyNum(round(sum(c(dat$N_HET[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)], 
                                       dat$N_HOMALT[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)]))/max(dat$N_BI_GENOS[which(dat$SVTYPE=="CPX")]), 
                                 1), big.mark=","))
    #Median size (text)
    textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)])/1000, 1), nsmall=1), "kb", sep=""))
    #Size distribution
    cpx.SVLEN.dens(dat, cpx.type)
    #Mean AF (text)
    textBox(text=paste(format(round(100*mean(dat$AF[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE==cpx.type)], na.rm=T), 2), nsmall=2), "%", sep=""))
  })
  
  #Axes,  if necessary
  par(bty="n", mar=c(0.1, 0.5, 0.5, 0.1), xpd=T)
  sapply(rep(0, 4), function(i){blankBox()})
  # pctAxis(max.value=max(table(dat$CPX_TYPE[which(dat$SVTYPE=="CPX")])/length(which(dat$SVTYPE=="CPX"))))
  sapply(rep(0, 2), function(i){blankBox()})
  sizeAxis()
}

################################
###MERGED SIMPLE + COMPLEX TABLE
################################
#Master function to make merged canonical + complex SV table
masterSVtableFigure <- function(dat.wrelateds, dat.all, dat, svtypes, med.sitesPerSample){
  #Set parameters
  filler.color <- "gray50"
  simple.rows <- c("All SV", "DEL", "DUP", "MCNV", "INS", "INV", "CTX", "CPX")
  complex.rows <- c("delINV\nINVdel", "delINVdel", 
                    "dupINV\nINVdup", "dupINVdup", 
                    "delINVdup\ndupINVdel", 
                    # "piDUP (FR)\npiDUP (RF)", 
                    "dDUP", "dDUP-iDEL\nINS-iDEL")
  cpx.groups <- list(c("delINV", "INVdel"), "delINVdel", 
                     c("dupINV", "INVdup"), "dupINVdup", 
                     c("delINVdup", "dupINVdel"), 
                     # c("piDUP_FR", "piDUP_RF"), 
                     "dDUP", c("dDUP_iDEL", "INS_iDEL"))
  #Order complex groups by size
  new.cpx.order <- order(-unlist(lapply(cpx.groups, function(classes){
    median(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE %in% classes)])
  })))
  complex.rows <- complex.rows[new.cpx.order]
  cpx.groups <- cpx.groups[new.cpx.order]
  sv.rows <- c(simple.rows, complex.rows, "BND")
  panels <- c("", "SV Class", "Abbrev.", "Mutational\nSignature", "", "Ref. Allele\nStructure", "Alt. Allele\nStructure(s)", 
              "Resolved\nVariants", "SV per\nGenome", "SV Size", "APS", "")
  n.panels <- length(panels)
  panel.widths=c(0.2, 0.7, 0.5, 0.4, 0.2, 0.65, 1, 0.5, 0.5, 0.6, 0.5, 0.2)
  panel.heights=c(0.9, rep(1, times=length(sv.rows)), 0.75)
  
  #Empty header text box with single value
  headerTextBox <- function(text, cex=1){
    plot(x=c(0, 1), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    if(text!=""){
      text(x=0.5, y=par("usr")[3], pos=3, labels=text, cex=cex, font=2)
      axis(1, at=c(par("usr")[1], par("usr")[2]), tck=0, labels=NA)
    }
  }
  
  #Jewel corresponding to SV class
  jewel <- function(svtype, svtypes){
    plot(x=c(0, 1), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    if(svtype=="ALL"){
      color <- "gray30"
    }else if(svtype %in% c("CTX", "BND")){
      color <- svtypes$color[which(svtypes$svtype=="OTH")]
    }else{
      color <- svtypes$color[which(svtypes$svtype==svtype)]
    }
    points(x=0.5, y=0.5, pch=21, bg=color, cex=2.5)
  }
  
  #Empty text box with single value
  textBox <- function(text, cex=1, pos=2, x=1, font=1, color="black"){
    plot(x=c(0, 1), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    text(x=x, y=0.5, pos=pos, labels=text, cex=cex, font=font, col=color, xpd=T)
  }
  
  #Blank box with nothing
  blankBox <- function(){
    plot(x=c(0, 1), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
  }
  
  #Get text corresponding to each sv class
  getClassText <- function(svtype){
    if(svtype=="ALL"){
      text <- "All SV"
    }else if(svtype=="DEL"){
      text <- "Deletion"
    }else if(svtype=="DUP"){
      text <- "Duplication"
    }else if(svtype=="MCNV"){
      text <- "Multiallelic CNV"
    }else if(svtype=="INS"){
      text <- "Insertion"
    }else if(svtype=="INV"){
      text <- "Inversion"
    }else if(svtype=="CTX"){
      text <- "Reciprocal\nTranslocation"
    }else if(svtype=="BND"){
      text <- "Breakend\n(Unresolved)"
    }else if(svtype=="CPX"){
      text <- "All Complex SVs"
    }
    return(text)
  }
  
  #Get text corresponding to each complex group
  getCPXSubclassText <- function(cpx.subclasses){
    if("delINV" %in% cpx.subclasses){
      text <- "Deletion-Flanked\nInversion"
    }else if("delINVdel" %in% cpx.subclasses){
      text <- "Paired-Deletion\nInversion"
    }else if("dupINV" %in% cpx.subclasses){
      text <- "Dup.-Flanked\nInversion"
    }else if("dupINVdup" %in% cpx.subclasses){
      text <- "Paired-Dup.\nInversion"
    }else if("delINVdup" %in% cpx.subclasses){
      text <- "Paired-Del./Dup.\nInversion"
    }else if("piDUP_FR" %in% cpx.subclasses){
      text <- "Palindromic\nInverted Dup."
    }else if("dDUP" %in% cpx.subclasses){
      text <- "Dispersed\nDuplication"
    }else if("dDUP_iDEL" %in% cpx.subclasses){
      text <- "Insertion with\nIns. Site Del."
    }else{
      text <- "N/A"
    }
    return(text)
  }
  
  #Boxes for SV types involved in CPX class
  signatures <- function(svtype="CPX", cpx.subclasses=NULL){
    cpx.type <- paste(cpx.subclasses, collapse="")
    plot(x=c(0, 4), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    rect(xleft=(0:3)+0.2, xright=(0:3)+0.8, 
         ybottom=0.2, ytop=0.8, col="gray95", border="gray90")
    if(length(grep("DEL", svtype, ignore.case=T))>0 
       | svtype=="MCNV"
       | length(grep("DEL", cpx.type, ignore.case=T))>0){
      rect(xleft=0.2, xright=0.8, ybottom=0.2, ytop=0.8, 
           col=svtypes$color[which(svtypes$svtype=="DEL")])
    }
    if(length(grep("DUP", svtype, ignore.case=T))>0 
       | svtype=="MCNV"
       | length(grep("DUP", cpx.type, ignore.case=T))>0){
      rect(xleft=1.2, xright=1.8, ybottom=0.2, ytop=0.8, 
           col=svtypes$color[which(svtypes$svtype=="DUP")])
    }
    if(length(grep("INS", svtype, ignore.case=T))>0 
       | length(c(grep("INS", cpx.type, ignore.case=T), 
                  grep("dDUP", cpx.type, ignore.case=T)))>0){
      rect(xleft=2.2, xright=2.8, ybottom=0.2, ytop=0.8, 
           col=svtypes$color[which(svtypes$svtype=="INS")])
    }
    if(length(grep("INV", svtype, ignore.case=T))>0 
       | length(c(grep("INV", cpx.type, ignore.case=T), 
                  grep("piDUP", cpx.type, ignore.case=T)))>0){
      rect(xleft=3.2, xright=3.8, ybottom=0.2, ytop=0.8, 
           col=svtypes$color[which(svtypes$svtype=="INV")])
    }
  }
  
  #Get count of sites,  pct pass,  and count of final variants
  getCounts <- function(dat.wrelateds, dat, svtype, cpx.subclasses=NULL){
    if(svtype=="ALL"){
      sites <- nrow(dat.wrelateds)
      # sites.pass <- length(which(dat.wrelateds$FILTER %in% c("PASS", "MULTIALLELIC")))
      final.sv <- nrow(dat)
    }else if(svtype=="CPX" & !is.null(cpx.subclasses)){
      sites <- length(which(dat.wrelateds$SVTYPE==svtype & dat.wrelateds$CPX_TYPE %in% cpx.subclasses))
      final.sv <- length(which(dat$SVTYPE==svtype & dat$CPX_TYPE %in% cpx.subclasses))
    }else{
      sites <- length(which(dat.wrelateds$SVTYPE==svtype))
      # sites.pass <- length(which(dat.wrelateds$SVTYPE==svtype & dat.wrelateds$FILTER %in% c("PASS", "MULTIALLELIC")))
      final.sv <- length(which(dat$SVTYPE==svtype))
    }
    # pct.pass <- sites.pass/sites
    pct.pass <- final.sv/sites
    return(c(sites, pct.pass, final.sv))
  }
  
  #Get count of median sites per sample
  getSitesPerSample <- function(med.sitesPerSample, svtype){
    if(svtype %in% c("CTX", "BND")){
      x <- 0
    }else if(svtype=="MCNV"){
      x <- sum(med.sitesPerSample[1, grep("MCNV_", colnames(med.sitesPerSample), fixed=T)])
    }else{
      x <- sum(med.sitesPerSample[1, which(colnames(med.sitesPerSample)==svtype)])
    }
    return(x)
  }
  
  #Simple barplot
  cpx.simpleBar <- function(value, max.value=1){
    plot(x=c(0, max.value), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    rect(xleft=0, xright=value, 
         ybottom=0.25, ytop=0.75, 
         col=filler.color)
  }
  
  #Simple size density plot
  SVLEN.dens <- function(dat, svtype, cpx.subclasses=NULL, font=1){
    require(zoo, quietly=T)
    if(svtype=="ALL"){
      sizes <- log10(dat$SVLEN)
      line.col <- "gray30"
    }else if(svtype=="CPX"){
      line.col <- svtypes$color[which(svtypes$svtype=="CPX")]
      if(is.null(cpx.subclasses)){
        sizes <- log10(dat$SVLEN[which(dat$SVTYPE=="CPX")])
      }else{
        sizes <- log10(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE %in% cpx.subclasses)])
      }
    }else{
      line.col <- svtypes$color[which(svtypes$svtype==svtype)]
      sizes <- log10(dat$SVLEN[which(dat$SVTYPE==svtype)])
    }
    # h <- hist(sizes, breaks=seq(0, 10, 0.05), plot=F)$counts
    # h <- rollapply(h, 11, mean, partial=T)
    # h <- h/max(h, na.rm=T)
    # out.df <- data.frame("min.size"=10^seq(0, 9.95, 0.05), 
    #                      "dens"=h)
    plot(x=log10(c(50, 10000000)), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    rect(xleft=log10(50), xright=log10(10000000), ybottom=0.1, ytop=0.6, col="white", border="gray85")
    segments(x0=0:10, x1=0:10, y0=0.1, y1=0.6, col="gray90")
    # segments(x0=log10(c(50, 10000000)), x1=log10(c(50, 10000000)), y0=0.1, y1=0.6, col="black", lwd=1.5)
    # abline(v=0:10, col="gray90")
    # polygon(x=log10(c(out.df$min.size, rev(out.df$min.size))), 
    #         y=c(out.df$dens, rep(0, times=nrow(out.df))), 
    #         col=filler.color)
    # axis(1, at=c(par("usr")[1], par("usr")[2]), tck=0, labels=NA)
    # abline(v=median(sizes, na.rm=T), col=line.col, lwd=3, lend="round")
    require(vioplot, quietly=T)
    viosizes <- as.numeric(sizes[which(!is.na(sizes) & !is.infinite(sizes) & !is.nan(sizes))])
    if(length(viosizes)>0){
      vioplot(viosizes, horizontal=T, add=T, col=line.col, at=0.35, drawRect=F, h=0.25, wex=0.5)
      segments(x0=median(sizes, na.rm=T), 
               x1=median(sizes, na.rm=T), 
               y0=0.15, y1=0.55, col="black", lwd=2, lend="round")
      text(x=median(sizes, na.rm=T), y=0.5, pos=3, cex=1, font=font, xpd=T, 
           labels=paste(format(round(median((10^sizes)/1000, na.rm=T), 1), n.small=1), "kb", sep=""))
    }
  }
  
  #Calculate proportion of singletons for a single SV class
  calc.fracSingletons.singleClass <- function(dat, svtype, cpx.subclasses=NULL, boot.n=100, conf=0.95, aps=T){
    dat <- dat[which(!(dat$chrom %in% c("X", "Y"))), ]
    if(svtype=="ALL"){
      ACs <- dat$AC
      v.aps <- dat$APS
    }else if(svtype=="CPX" & !is.null(cpx.subclasses)){
      ACs <- dat[which(dat$SVTYPE==svtype & dat$CPX_TYPE %in% cpx.subclasses), ]$AC
      v.aps <- dat[which(dat$SVTYPE==svtype & dat$CPX_TYPE %in% cpx.subclasses), ]$APS
    }else{
      ACs <- dat[which(dat$SVTYPE==svtype), ]$AC
      v.aps <- dat[which(dat$SVTYPE==svtype), ]$APS
    }
    ACs <- ACs[which(!is.na(ACs) & ACs>0)]
    v.aps <- v.aps[which(!is.na(v.aps))]
    if(aps==T){
      if(length(ACs) > 0){
        helper.getFracSingle <- function(v.aps, indices){sum(v.aps[indices])/length(v.aps[indices])}
        point.est <- helper.getFracSingle(v.aps, indices=1:length(v.aps))
        calc.ci <- function(v.aps, n, conf){
          set.seed(0)
          boot.obj <- boot(data=v.aps, statistic=helper.getFracSingle, R=n)
          ci <- boot.ci(boot.obj, conf=conf, type="basic")$basic[4:5]
          return(ci)
        }
        ci <- calc.ci(v.aps, n=boot.n, conf=conf)
        if(length(ci) != 2){
          ci <- c(NA, NA)
        }
      }else{
        point.est <- NA
        ci <- c(NA, NA)
      }
    }else{
      if(length(ACs) > 0){
        helper.getFracSingle <- function(ACs, indices){length(which(ACs[indices]==1))/length(ACs[indices])}
        point.est <- helper.getFracSingle(ACs, indices=1:length(ACs))
        calc.ci <- function(ACs, n, conf){
          set.seed(0)
          boot.obj <- boot(data=ACs, statistic=helper.getFracSingle, R=n)
          ci <- boot.ci(boot.obj, conf=conf, type="basic")$basic[4:5]
          return(ci)
        }
        ci <- calc.ci(ACs, n=boot.n, conf=conf)
      }else{
        point.est <- NA
        ci <- c(NA, NA)
      }
    }
    return(c(point.est, ci))
  }
  
  #Plot point estimate & 95% CI for proportion singletons
  propSingles <- function(dat, svtype, cpx.subclasses=NULL, font=1){
    if(svtype=="ALL"){
      point.col <- "gray30"
    }else if(svtype=="CPX"){
      point.col <- svtypes$color[which(svtypes$svtype=="CPX")]
    }else if(svtype %in% c("CTX", "BND")){
      point.col <- svtypes$color[which(svtypes$svtype=="OTH")]
    }else{
      point.col <- svtypes$color[which(svtypes$svtype==svtype)]
    }
    s.dat <- calc.fracSingletons.singleClass(dat, svtype=svtype, cpx.subclasses=cpx.subclasses, aps=T)
    plot(x=c(-0.25, 0.25), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    # abline(v=seq(0, 1, 0.25), col="gray90")
    # abline(v=c(0, 1))
    rect(xleft=-0.25, xright=0.25, ybottom=0.1, ytop=0.6, col="white", border="gray85")
    segments(x0=seq(-0.25, 0.25, 0.1), x1=seq(-0.25, 0.25, 0.1), 
             y0=0.1, y1=0.6, col="gray90")
    segments(x0=0, x1=0, y0=0.1, y1=0.6)
    # segments(x0=c(-0.25, 0.25), x1=c(-0.25, 0.25), 
    #          y0=0.1, y1=0.6, col="black", lwd=1.5)
    segments(x0=s.dat[2], x1=s.dat[3], y0=-0.255, y1=-0.255, lend="round", lwd=1.5)
    points(x=s.dat[1], y=0.4, pch=21, col="black", bg=point.col, cex=1.5)
    # rect(xleft=s.dat[1]-0.015, xright=s.dat[1]+0.015, 
    #      ybottom=0.15, ytop=0.55, col=point.col, lwd=1)
    # segments(x0=s.dat[2], x1=s.dat[3], y0=-0.255, y1=-0.255, lend="round", lwd=1.5)
    # points(x=s.dat[1], y=-0.255, pch=21, bg=point.col, cex=1.5)
    text(x=s.dat[1], y=0.5, pos=3, cex=1, font=font, xpd=T, 
         labels=formatC(round(s.dat[1], 2), digits=2), n.small=1)
  }
  
  #Pct axis
  pctAxis <- function(range=c(-0.25, 0.25), ticks.at=seq(-0.25, 0.25, 0.1), labels.at=seq(-0.25, 0.25, 0.2)){
    plot(x=range, y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    axis(3, at=ticks.at, labels=NA, tck=0.15)
    sapply(labels.at, function(x){
      axis(3, at=x, tick=F, cex.axis=0.75, las=2, hadj=1, line=-1.6)
    })
    # mtext(3, text="Pct.", line=-2.5, cex=0.7)
  }
  
  #Size axis
  sizeAxis <- function(){
    plot(x=log10(c(50, 10000000)), y=c(0, 1), type="n", 
         xaxt="n", xlab="", xaxs="i", 
         yaxt="n", ylab="", yaxs="i")
    axis(3, at=as.numeric(unlist(sapply(1:7, function(x){log10((1:9)*(10^x))}))), labels=NA, tck=0.075, lwd=0.5)
    axis(3, at=1:8, labels=NA, tck=0.15)
    labels <- c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb")
    sapply(seq(2, 7, by=1), function(i){
      axis(3, at=i, tick=F, labels=labels[i], cex.axis=0.75, las=2, hadj=1, line=-1.6)
    })
    # mtext(3, text="SV Size", line=-2.5, cex=0.7)
  }
  
  #Prep figure layout
  layout(matrix(1:(n.panels*(length(sv.rows)+2)), byrow=T, nrow=length(sv.rows)+2), 
         widths=panel.widths, heights=panel.heights)
  
  #Headers
  par(bty="n", mar=c(0.1, 1, 0.1, 1), xpd=T)
  sapply(panels, function(title){
    headerTextBox(title)
  })
  
  #Add top row for summary of all SV
  par(bty="n", mar=c(0.4, 0.5, 0.4, 0.5), xpd=F)
  jewel("ALL", svtypes)
  textBox(text=getClassText("ALL"), font=2)
  textBox(text="ALL", font=2)
  textBox("Varies", color="gray50", pos=NULL, x=0.5, font=3)
  blankBox()
  textBox("Varies", color="gray50", pos=NULL, x=0.5, font=3)
  textBox("Varies", color="gray50", pos=NULL, x=0.5, font=3)
  # textBox(text=prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype="ALL")[1], big.mark=","), font=2)
  # textBox(paste(format(round(100*getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype="ALL")[2], 1), n.small=1), "%", sep=""), font=2)
  textBox(text=paste(prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype="ALL")[3], big.mark=","), "*", sep=""), font=2)
  # textBox(text=prettyNum(round(sum(c(dat$N_HET, dat$N_HOMALT), na.rm=T)/max(dat$N_BI_GENOS, na.rm=T), 
  #                              1), big.mark=","))
  textBox(prettyNum(getSitesPerSample(med.sitesPerSample, svtype="ALL"), big.mark=","), font=2)
  SVLEN.dens(dat, svtype="ALL", cpx.subclasses=NULL, font=2)
  # textBox(text=paste(format(round(median(dat$SVLEN)/1000, 1), nsmall=1), "kb", sep=""), font=2)
  propSingles(dat, svtype="ALL", font=2)
  blankBox()
  
  #Iterate over canonical classes
  sapply(c("DEL", "DUP", "MCNV", "INS", "INV", "CTX"), function(svtype){
    #Jewel
    jewel(svtype, svtypes)
    #SV class (text)
    textBox(getClassText(svtype))
    #Abbreviation (text)
    textBox(text=svtype)
    #Signatures
    if(svtype %in% c("CTX", "BND")){
      textBox("N/A", color="gray50", pos=NULL, x=0.5)
    }else{
      signatures(svtype)
    }
    blankBox()
    #Placeholder for rearrangement schematic (done in illustrator)
    blankBox()
    blankBox()
    # textBox(text=prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype=svtype)[1], big.mark=","))
    # textBox(paste(format(round(100*getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype=svtype)[2], 1), n.small=1), "%", sep=""))
    textBox(text=prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype=svtype)[3], big.mark=","))
    #Median # of variants per genome
    textBox(prettyNum(getSitesPerSample(med.sitesPerSample, svtype=svtype), big.mark=","))
    #Median size (text) & size distribution
    if(svtype %in% c("CTX", "BND")){
      textBox(text="N/A", color="gray50", pos=NULL, x=0.5)
      # textBox(text="N/A")
    }else{
      SVLEN.dens(dat, svtype=svtype, cpx.subclasses=NULL)
      # textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE==svtype)])/1000, 1), nsmall=1), "kb", sep=""))
    }
    if(svtype %in% c("MCNV", "BND", "CTX")){
      textBox("N/A", color="gray50", pos=NULL, x=0.5)
    }else{
      propSingles(dat, svtype=svtype)
    }
    blankBox()
  })
  
  #One summary row for all CPX SV
  jewel("CPX", svtypes)
  textBox(getClassText("CPX"), font=2)
  textBox(text="CPX", font=2)
  textBox(text="Varies", color="gray50", pos=NULL, x=0.5, font=3)
  blankBox()
  textBox(text="Varies", color="gray50", pos=NULL, x=0.5, font=3)
  textBox(text="Varies", color="gray50", pos=NULL, x=0.5, font=3)
  # textBox(text=prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype="CPX")[1], big.mark=","), font=2)
  # textBox(paste(format(round(100*getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype="CPX")[2], 1), n.small=1), "%", sep=""), font=2)
  textBox(text=prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype="CPX")[3], big.mark=","), font=2)
  textBox(prettyNum(getSitesPerSample(med.sitesPerSample, svtype="CPX"), big.mark=","), font=2)
  SVLEN.dens(dat, svtype="CPX", cpx.subclasses=NULL, font=2)
  # textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE=="CPX")])/1000, 1), nsmall=1), "kb", sep=""))
  propSingles(dat, svtype="CPX", font=2)
  blankBox()
  
  #Iterate over cpx types - one row per group of types
  sapply(1:length(cpx.groups), function(i){
    #Get parameters
    cpx.group <- complex.rows[i]
    cpx.subclasses <- cpx.groups[[i]]
    #Jewel
    jewel(svtype="CPX", svtypes)
    #Category name
    textBox(text=getCPXSubclassText(cpx.subclasses), cex=1)
    #Abbrev.
    if(length(cpx.subclasses)>1){
      textBox(text=cpx.group, cex=1)
    }else{
      textBox(text=cpx.group)
    }
    #Mutation signature
    signatures(svtype="CPX", cpx.subclasses=cpx.subclasses)
    blankBox()
    #Spacer for alt allele schematics
    sapply(rep(0, 2), function(i){blankBox()})
    #Counts of sites
    # textBox(text=prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, 
    #                                  svtype="CPX", cpx.subclasses=cpx.subclasses)[1], big.mark=","))
    # textBox(paste(format(round(100*getCounts(dat=dat, dat.wrelateds=dat.wrelateds, 
    #                                          svtype="CPX", cpx.subclasses=cpx.subclasses)[2], 1), n.small=1), "%", sep=""))
    textBox(text=prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, 
                                     svtype="CPX", cpx.subclasses=cpx.subclasses)[3], big.mark=","))
    #TODO: collect per-sample counts of each complex type
    textBox("TBD")
    #Size distribution
    SVLEN.dens(dat, svtype="CPX", cpx.subclasses=cpx.subclasses)
    # textBox(text=paste(format(round(median(dat$SVLEN[which(dat$SVTYPE=="CPX" & dat$CPX_TYPE %in% cpx.subclasses)])/1000, 1), nsmall=1), "kb", sep=""))
    propSingles(dat, svtype="CPX", cpx.subclasses=cpx.subclasses)
    #Final spacer
    blankBox()
  })
  
  #Add final row for BNDs  
  #Jewel
  jewel(svtype="BND", svtypes)
  #SV class (text)
  textBox(getClassText("BND"))
  #Abbreviation (text)
  textBox(text="BND")
  #Signatures
  textBox("N/A", color="gray50", pos=NULL, x=0.5)
  blankBox()
  #Placeholder for rearrangement schematic (done in illustrator)
  blankBox()
  blankBox()
  # textBox(text=prettyNum(getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype="BND")[1], big.mark=","))
  # textBox(paste(format(round(100*getCounts(dat=dat, dat.wrelateds=dat.wrelateds, svtype="BND")[2], 1), n.small=1), "%", sep=""))
  textBox(text=paste(prettyNum(getCounts(dat=dat.all, dat.wrelateds=dat.wrelateds, svtype="BND")[3], big.mark=","), "*", sep=""))
  #Median # of variants per genome
  # textBox(prettyNum(getSitesPerSample(med.sitesPerSample, svtype="BND"), big.mark=","))
  textBox("N/A", color="gray50", font=3)
  #Median size (text) & size distribution
  textBox(text="N/A", color="gray50", pos=NULL, x=0.5)
  textBox("N/A", color="gray50", pos=NULL, x=0.5)
  blankBox()
  
  #Axes,  if necessary
  par(bty="n", mar=c(0.1, 0.5, 0.1, 0.5), xpd=T)
  sapply(rep(0, 9), function(i){blankBox()})
  # pctAxis(max.value=max(table(dat$CPX_TYPE[which(dat$SVTYPE=="CPX")])/length(which(dat$SVTYPE=="CPX"))))
  # sapply(rep(0, 1), function(i){blankBox()})
  sizeAxis()
  pctAxis()
}

#Doubleton analysis table
doubleton.analysis <- function(dat, sd_sr_cov, max_sd_sr_cov=0.1){
  #Identify doubletons
  doubletons <- dat[which(dat$N_HET==2 & dat$AC==2 & !(dat$chrom %in% c("X", "Y"))), ]
  doubletons <- doubletons[which(doubletons$name %in% sd_sr_cov[which(sd_sr_cov[, 2] <= max_sd_sr_cov), 1]), ]
  nhet.pop.idxs <- which(colnames(doubletons) %in% paste(pops$pop[which(pops$pop!="OTH")], "_N_HET", sep=""))
  doubletons.samepop <- doubletons[which(apply(doubletons[, nhet.pop.idxs], 1, max)==2), ]
  doubletons.crosspop <- doubletons[which(apply(doubletons[, nhet.pop.idxs], 1, max)==1
                                          & apply(doubletons[, nhet.pop.idxs], 1, sum)==2), ]
  doubletons <- doubletons[which(doubletons$name %in% c(doubletons$name,
                                                        doubletons$name)),]
  #Iterate over svtypes and compute fraction of samepop vs crosspop doubletons
  doubleton.table <- sapply(c("DEL", "DUP", "INS", "INV", "CPX"), function(svtype){
    all <- length(which(doubletons$SVTYPE==svtype))
    samepop <- length(which(doubletons.samepop$SVTYPE==svtype))
    crosspop <- length(which(doubletons.crosspop$SVTYPE==svtype))
    frac.cross <- crosspop/all
    return(c(all, samepop, crosspop, frac.cross))
  })
  doubleton.table <- cbind("Metric"=c("all", "samepop", "crosspop", "frac.cross"),
                           doubleton.table)
  return(doubleton.table)
}

#Get table of gross chromosomal abnormalities
gather.bca.table <- function(dat, max.AF=0.01){
  res <- as.data.frame(t(sapply(c("DEL", "DUP", "INV", "CTX", "CPX"), function(svtype){
    idxs <- which(dat$SVTYPE==svtype 
                  & !(dat$chrom %in% c("X", "Y"))
                  & dat$AF<max.AF 
                  & (dat$SVLEN>=1000000 | dat$SVTYPE=="CTX"))
    # if(svtype=="CPX"){
    #   idxs <- unique(c(idxs, 
    #                    which(dat$SVTYPE==svtype 
    #                          & !(dat$chrom %in% c("X", "Y"))
    #                          & dat$AF<max.AF 
    #                          & dat$CPX_TYPE=="CCR")))
    # }
    variants <- length(idxs)
    carriers <- sum(dat$N_HET[idxs])
    #One CCR is actually a complex translocation,  and should be counted
    if(svtype=="CPX"){carriers <- carriers+1; variants <- variants+1}
    # if(svtype=="CTX"){carriers <- 17; variants <- 15}
    denom <- max(dat$N_BI_GENOS[idxs], na.rm=T)
    rate <- carriers/denom
    binom.ci <- binom.test(carriers, denom)$conf.int
    c(variants, carriers, rate, binom.ci)
  })))
  colnames(res) <- c("variants", "carriers", "mean", "lower", "upper")
  res <- rbind(apply(res, 2, sum), res)
  rownames(res)[1] <- "ALL"
  return(res)
}

#Plot basic metadata for gross chromosomal abnormalities
plot.chrom.abnormalities <- function(dat, max.AF=0.01, 
                                     category.labels=c("Deletion", 
                                                       "Duplication", 
                                                       "Inversion", 
                                                       "Reciprocal\nTranslocation", 
                                                       "Complex SV"),
                                     cex.toplabels=0.9){
  #Get plotting data
  plot.dat <- gather.bca.table(dat=dat, max.AF=max.AF)
  plot.dat <- plot.dat[-1, ]
  xmax <- max(plot.dat[, -c(1:2)])
  point.colors <- as.character(sapply(rownames(plot.dat), function(svtype){
    if(svtype=="CTX"){
      svtype <- "OTH"
    }
    svtypes$color[which(svtypes$svtype==svtype)]
  }))
  #Prep plot area
  layout(matrix(c(1:4), nrow=1, byrow=T), 
         widths=c(4, 1.5, 3, 5))
  #Plot class label
  par(mar=c(0.3, 0.3, 3, 0.3), bty="n", xpd=T)
  plot(x=c(0, 1), y=c(0, -nrow(plot.dat)), type="n", 
       xlab="", ylab="", xaxt="n", yaxt="n")
  text(x=1, y=(-1:-nrow(plot.dat))+0.5, pos=2, 
       labels=category.labels)
  mtext(3, line=1.2, text="Rare SVs", cex=cex.toplabels)
  mtext(3, line=0, text=(bquote("" >= "1Mb")), cex=cex.toplabels)
  axis(3, at=c(0, 1), labels=NA, tck=0)
  #Plot number of variants
  par(xpd=T)
  plot(x=c(0, 1), y=c(0, -nrow(plot.dat)), type="n", 
       xlab="", ylab="", xaxt="n", yaxt="n")
  text(x=0.5, y=(-1:-nrow(plot.dat))+0.5, cex=1, 
       labels=plot.dat$variants)
  mtext(3, line=0, text="SVs", cex=cex.toplabels)
  axis(3, at=c(0, 1), labels=NA, tck=0)
  #Plot number of carriers
  plot(x=c(0, 1), y=c(0, -nrow(plot.dat)), type="n", 
       xlab="", ylab="", xaxt="n", yaxt="n")
  text(x=0.5, y=(-1:-nrow(plot.dat))+0.5, cex=0.95, 
       labels=paste(plot.dat$carriers, "/", prettyNum(max(dat$N_BI_GENOS, na.rm=T), big.mark=",")))
  mtext(3, line=0, text="Samples", cex=cex.toplabels)
  axis(3, at=c(0, 1), labels=NA, tck=0)
  #Plot points & CIs
  par(mar=c(0.3, 0.4, 3, 0.3), bty="n", xpd=F)
  plot(x=c(-0.05*xmax, xmax), y=c(0, -nrow(plot.dat)), type="n", 
       xlab="", ylab="", xaxt="n", yaxt="n")
  abline(v=axTicks(3), col="gray90")
  abline(v=0, col="gray50")
  segments(x0=plot.dat$lower, x1=plot.dat$upper, 
           y0=(-1:-nrow(plot.dat))+0.5, 
           y1=(-1:-nrow(plot.dat))+0.5, 
           lend="round", lwd=2)
  points(x=plot.dat$mean, 
         y=(-1:-nrow(plot.dat))+0.5, 
         pch=21, cex=1.5, bg=point.colors)
  par(xpd=T)
  text(x=plot.dat$mean, y=(-1:-nrow(plot.dat))+0.5, 
       pos=3, labels=paste(format(round(100*plot.dat$mean, 2), nsmall=2), "%", sep=""))
  par(xpd=F)
  axis(3, at=axTicks(3), labels=NA, tck=-0.025)
  sapply(axTicks(3)[seq(1, length(axTicks(3)), 2)], function(x){
    axis(3, at=x, tick=F, line=-0.8, cex.axis=0.9, 
         labels=paste(round(100*x, 1), "%", sep=""))
  })
  mtext(3, line=1.1, text="Samples (%)", cex=cex.toplabels)
  # #Plot schematic label
  # plot(x=c(0, 1), y=c(0, -nrow(plot.dat)), type="n", 
  #      xlab="", ylab="", xaxt="n", yaxt="n")
  # mtext(3, line=-0.9, text="Schematic")
}


# Plot QUAL score stacked barplots for a single SV type, stratified by size
plot.qualCDF.bySize.single <- function(qualscores, dat, svtype, 
                                       size.bins=c(-1, 1000, 10000, 500000000),
                                       qual.bins=seq(0, 1000, 50),
                                       legend=F, plot.xlab=T, plot.ylab=T){
  #Helper function to gather count of quals per qual bin
  get.qual.counts <- function(quals){
    sapply(2:length(qual.bins), function(k){
      return(length(which(quals > qual.bins[k-1] & quals <= qual.bins[k])))
    })
  }
  #Get qual counts per size bin
  counts.per.bin <- sapply(2:length(size.bins), function(i){
    quals.sub <- qualscores[which(qualscores$ID %in% dat$name[which(dat$SVTYPE==svtype 
                                                                    & dat$SVLEN >= size.bins[i-1] 
                                                                    & dat$SVLEN < size.bins[i])]), 2]
    get.qual.counts(quals.sub)
  })
  bin.sums <- apply(counts.per.bin, 1, sum)
  #Prep plot area and parameters
  par(bty="n", mar=c(2.8, 4.8, 0.2, 0.2))
  plot(x=range(qual.bins), y=c(0, max(bin.sums)), type="n",
       xaxt="n", xlab="", yaxt="n", ylab="")
  base.color <- svtypes$color[which(svtypes$svtype==svtype)]
  colors <- sapply(c(1/3, 2/3, 1), function(a){adjustcolor(base.color, alpha=a)})
  #Plot bars
  sapply(2:length(qual.bins), function(x){
    rect(xleft=qual.bins[x-1], xright=qual.bins[x],
         ybottom=cumsum(c(0, counts.per.bin[x-1, 1:2])),
         ytop=cumsum(counts.per.bin[x-1, ]),
         col=colors,
         bty="n", border=NA)
  })
  #Add axes
  x.at <- seq(0, 1000, 250)
  axis(1, at=x.at, labels=NA, tck=-0.03)
  axis(1, at=x.at, tick=F, line=-0.6, labels=x.at)
  if(plot.xlab==T){
    mtext(1, text="QUAL Score", line=1.7, cex=0.9)
  }
  if(length(axTicks(2)) > 5){
    y.at <- axTicks(2)[seq(1, length(axTicks(2)), 2)]
  }else{
    y.at <- axTicks(2)
  }
  axis(2, at=y.at, labels=NA, tck=-0.03)
  axis(2, at=y.at, tick=F, las=2, line=-0.6,
       labels=prettyNum(y.at, big.mark=","))
  if(plot.ylab==T){
    mtext(2, text=paste(svtype, "Count"), line=3.4, cex=0.9)
  }
  #Add legend
  if(legend==T){
    legend("top", legend=rev(c("<1kb", "1-10kb", ">10kb")),
           fill=rev(colors), border=NA, bty="n", cex=0.75)
  }
}


#Plot strip of QUAL score distributions, with no prespecified layout
plot.qualDistribs.strip <- function(qualscores, dat, svtype, 
                                    freq.bins=c(0, 0.01, 0.1, 1),
                                    size.bins=c(-1, 1000, 10000, 500000000),
                                    qual.bins=seq(0, 1000, 50),
                                    legend=F, plot.xlab=T){
  #First panel = ALL variants of a given SV type
  plot.qualCDF.bySize.single(qualscores, dat, svtype, 
                             size.bins=size.bins, qual.bins=qual.bins,
                             legend=legend, plot.ylab=T, plot.xlab=plot.xlab)
  #Remaining panels: split by frequency bin
  sapply(2:length(freq.bins), function(i){
    plot.qualCDF.bySize.single(qualscores, dat[which(dat$AF > freq.bins[i-1]
                                                     & dat$AF <= freq.bins[i]), ], 
                               svtype, size.bins=size.bins, qual.bins=qual.bins,
                               legend=legend, plot.ylab=F, plot.xlab=plot.xlab)
  })
}


#Plot QUAL score distributions by SV type, frequency, and size
plot.qualDistribs <- function(qualscores, dat, 
                              v.svtypes=c("DEL", "DUP", "INS", "INV", "CPX"),
                              freq.bins=c(0, 0.01, 0.1, 1),
                              size.bins=c(-1, 1000, 10000, 500000000),
                              qual.bins=seq(0, 1000, 50)){
  #Helper function for svtype-coded legend
  plot.qual.legend <- function(svtype){
    par(bty="n", mar=c(3, rep(0.5, 3)))
    plot(x=c(0, 1), y=c(0, 1), type="n",
         xaxt="n", xlab="", yaxt="n", ylab="")
    base.color <- svtypes$color[which(svtypes$svtype==svtype)]
    colors <- sapply(c(1/3, 2/3, 1), function(a){adjustcolor(base.color, alpha=a)})
    legend("center", legend=rev(c("<1kb", "1-10kb", ">10kb")),
           fill=rev(colors), border=NA, bty="n")
  }
  layout(matrix(seq(1, (length(freq.bins) + 1) * length(v.svtypes)),
                byrow=T, nrow=length(v.svtypes)),
         widths=c(rep(1, length(freq.bins)), 2/3))
  # par(mfrow=c(length(v.svtypes), length(freq.bins) + 1))
  sapply(v.svtypes, function(svtype){
    if(which(v.svtypes==svtype) == length(v.svtypes)){
      plot.xlab <- T
    }else{
      plot.xlab <- F
    }
    plot.qualDistribs.strip(qualscores, dat[which(dat$FILTER=="PASS"), ],
                            svtype=svtype, freq.bins=freq.bins,
                            size.bins=size.bins, qual.bins=qual.bins,
                            plot.xlab=plot.xlab)
    plot.qual.legend(svtype)
  })
}

