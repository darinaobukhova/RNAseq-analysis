#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

process <- args[1]
config_file_run <- args[2] #file with specifics (outputdir, scripts, etc.,)

source(as.character(config_file_run))

#---Loading libraries---#

library(Rsamtools)
library(gplots) 
library(Biostrings)
library(edgeR)

#---Functions---#

countTable_DIR_functions <- paste(script, "functions", sep="/")
cat("Loading sources...")
countTable_functions <- list.files(pattern="[.]R$", path=countTable_DIR_functions);
sapply(paste(countTable_DIR_functions, countTable_functions, sep="/"), FUN=source)

sampleID <- read.csv2(paste(analyses, log, "samplesheet.csv", sep="/"), header = T, stringsAsFactors = F)
  
cat(sampleID[, "Sample_Project"], "\n")

ref_genomes <- unique(sampleID[, "genome"])

for (ref_genome in ref_genomes){
  sampleID2    <- sampleID[sampleID[,"genome"]%in%ref_genome,]
  flsnames    <- sampleID2[,"SampleID"]
  cat(flsnames, "\n")
  colflsnames <- paste(sampleID2[,"Sample_Project"],sampleID2[,"SampleID"],sep=".")
  cat(colflsnames, "\n")
  grps        <- unique(sampleID2[,"Sample_Project"])
  cat(grps, "\n")
  bamfile=NULL
  if ( project=="Overview" | project=="ResearchLetter" | project=="FIB_ORG_comparison"  ) { 
    bamfile <-  paste(analyses,log,"bam",paste(flsnames,".bam",sep=""),sep="/") 
  } else { 
    for (grp in grps){
      bamfile_grp <- as.vector(paste(analyses,grp,"bam",paste(flsnames[grep(grp,flsnames)],".bam",sep=""),sep="/"))
      bamfile     <-c(bamfile,bamfile_grp)
    }
  }

}


groups      <- as.factor(sampleID2[,"Sample_Project"])
groups.lvl  <- levels(groups)
Biorep      <- summary(groups)

ref_gtf <- "RefSeq"
minreads <- 5
CPM <- 0.25
#Minimum samples present (25% of total number of samples)
MIN <- 6

if (project=="Overview" | project="ResearchLetter") {
  setwd(analyses)
  analyses.result <- paste(analyses, log, "results", sep="/")
  dir.create(analyses.result)
} else {
  setwd(analyses)
  analyses.result <- paste(analyses, process, sep="/")
  dir.create(analyses.result)
}

setwd(analyses.result)

#---GTF files---#

gtf <- gtf.selection(ref_genome, ref_gtf)
gtf2 <- gtf.symbol.selection(ref_genome,ref_gtf)
gtf2 <- gtf2[!duplicated(gtf2),]
gtf  <- gtf[order(gtf["chr"],gtf[,"start"],gtf[,"end"],gtf[,"ID"]),] 

#---Analysing regions of interest---#

cat("Analysis started...")

exons.allchr.whole <- gtf[,c("chr","start","end","ID")]

exons.allchr <- exons.allchr.whole
seqln <- getSequenceLengths(bamfile[[1]])

if (length(seqln)<2) {
  exons.allchr <- exons.allchr.whole; exons.allchr$chr <- sapply(exons.allchr[,"chr"],function(w) {unlist(strsplit(w,"chr"))[2]  }); exons.allchr[exons.allchr[,"chr"]=="M","chr"] <- "MT"}

lst.exons.allchr <- lapply(unique(exons.allchr[,"chr"]),function(x,D){ 
  k  <- which(D[,"chr"]==x);
  IRanges(start=D[k,"start"],end=D[k,"end"],names=D[k,"ID"]);}
  ,D=exons.allchr)

names(lst.exons.allchr) <- unique(exons.allchr[,"chr"])

rl.exons.allchr <- IRangesList()
for(x in names(lst.exons.allchr)){rl.exons.allchr[[x]]<-lst.exons.allchr[[x]]}

# genomic ranges
seqln <- getSequenceLengths(bamfile[1])
which2 <- IRangesList()
for( S in 1:length(seqln) ) {
  which2[[ names(seqln)[S] ]] <- IRanges( start=1, end=seqln[S] )
}

# bam parameters
param <- ScanBamParam(what=c('pos'), which=which2 )

# retrieve counts
totalcount <- getCounts(bamfile[[1]], targets=rl.exons.allchr, param=param, which=which2)
cat(bamfile[[1]], "\n")

for(bams in bamfile[-1]){
  counts       <- getCounts(bams, targets=rl.exons.allchr, param=param, which=which2)
  totalcount   <- data.frame(totalcount,counts,stringsAsFactors=FALSE)
  cat(bams, "\n")
}  

for (bam in bamfile){
  print(bam)
}


totalcount               <- data.frame(ID=exons.allchr[,"ID"],totalcount,stringsAsFactors=FALSE)
colnames(totalcount)[-1] <- colflsnames
totalcount               <- totalcount[,c(1,order(colflsnames)+1)]
colflsnames              <- colflsnames[order(colflsnames)]
totalcount2              <- totalcount; 
totalcount2$ID           <- gtf$ID.exon;
write.table(totalcount2,file=paste(analyses.result,"all_exons_with_counts.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE);


totalGenes <- aggregate(totalcount[,colflsnames],list(totalcount[,"ID"]),FUN=sum); 
colnames(totalGenes) <- c("ID",colflsnames); 

save.image(paste(analyses.result,paste(ref_genome,ref_gtf,"loaded_bams.Rdata",sep="_"),sep="/"))

conditions <- c("expr.groups","expr.all","total","high")
list_exp <- list()
list_exp[["total.allRNA"]]       <- totalGenes
list_exp[["total.mRNA"]]         <- list_exp[["total.allRNA"]][grepl("NM",list_exp[["total.allRNA"]][,"ID"]),]
list_exp[["total.ncRNA"]]        <- list_exp[["total.allRNA"]][grepl("NR",list_exp[["total.allRNA"]][,"ID"]),]

list_exp[["high.allRNA"]]        <- totalGenes[rowSums(cpm(totalGenes[,-1])>CPM)==ncol(totalGenes[,-1]),]
list_exp[["high.mRNA"]]          <- list_exp[["high.allRNA"]][grepl("NM",list_exp[["high.allRNA"]][,"ID"]),]
list_exp[["high.ncRNA"]]         <- list_exp[["high.allRNA"]][grepl("NR",list_exp[["high.allRNA"]][,"ID"]),]

list_exp[["expr.all.allRNA"]]    <- totalGenes[rowSums(totalGenes[,-1]>minreads)==ncol(totalGenes[,-1]),]
list_exp[["expr.all.mRNA"]]      <- list_exp[["expr.all.allRNA"]][grepl("NM",list_exp[["expr.all.allRNA"]][,"ID"]),]
list_exp[["expr.all.ncRNA"]]     <- list_exp[["expr.all.allRNA"]][grepl("NR",list_exp[["expr.all.allRNA"]][,"ID"]),]

list_exp[["expr.groups.allRNA"]] <- combine.groups.expression(totalGenes,groups,minreads)
list_exp[["expr.groups.mRNA"]]   <- list_exp[["expr.groups.allRNA"]][grepl("NM",list_exp[["expr.groups.allRNA"]][,"ID"]),]
list_exp[["expr.groups.ncRNA"]]  <- list_exp[["expr.groups.allRNA"]][grepl("NR",list_exp[["expr.groups.allRNA"]][,"ID"]),]

if (length(groups.lvl)>1){
  lvls <- t(combn(groups.lvl,2))
  
  for(condition in conditions){
    condition_names <- names(list_exp)[grepl(condition,names(list_exp))]
    for(nms in condition_names){
      dir.create(paste(analyses.result,nms,sep="/"))
      
      write.table(list_exp[[nms]],file=paste(analyses.result,nms,"all_genes_with_counts.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE)
      write.table(data.frame(ID=list_exp[[nms]],cpm(list_exp[[nms]][,-1], lib.size=NULL, log=TRUE, prior.count=0.25),stringsAsFactors=FALSE),file=paste(analyses.result,nms,"all_genes_with_counts_cpm_log2.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE)
      plots.input(list_exp[[nms]],input=nms,paste(analyses.result,nms,sep="/"),groups)
      
      if (grepl("expr",nms)) {expr.overlap(list_exp[[nms]],input=nms,paste(analyses.result,nms,sep="/"),groups)}
      for(i in 1:(nrow(lvls))){
        dfs   <- list_exp[[nms]][,grepl(lvls[i,1],colnames(list_exp[[nms]]))|grepl(lvls[i,2],colnames(list_exp[[nms]]))]
        dfs   <- data.frame(ID=list_exp[[nms]][,"ID"],dfs,stringsAsFactors=FALSE)
        dfs   <- dfs[!(rowSums(dfs[,-1]==minreads)==ncol(dfs[,-1])),]
        biorep <- c(length(groups[groups%in%lvls[i,1]]),length(groups[groups%in%lvls[i,2]]))
        Groups <- c(groups[groups%in%lvls[i,1]],groups[groups%in%lvls[i,2]])
        input <- paste(nms,lvls[i,1],"vs",lvls[i,2],sep="_")
        #stats.plus.plots.symbol(dfs, input, biorep, Groups,paste(analyses.result,nms,sep="/"),gtf2)
        #stats.plus.plots.batch(dfs, input, biorep, Groups,paste(analyses.result,nms,sep="/"),gtf2)
      } 
    }
  }
}
save.image(paste(analyses.result,paste(ref_genome,ref_gtf,"full_analyses.Rdata",sep="_"),sep="/"))






