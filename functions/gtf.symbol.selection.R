gtf.symbol.selection <- function(ref_genome,ref_gtf){
  
  if (ref_genome=="Hg38"){ 
    if(ref_gtf=="RefSeq"){ 
      gtf2 <- read.table("/ifs/software/inhouse/GSM/mRNA/reference/versions/1.0/Hg38_UCSC_RefSeq.gtf",sep="\t",stringsAsFactors=FALSE);
      gtf2 <- gtf2[,c(1,10)];
      names(gtf2) <- c("ID","symbol")
      } 
    else {gtf2 <- gtf.whole[,c("ID","symbol","RefSeq")]
    }
  } 
  
  if (ref_genome=="Hg19"){ 
    gtf2 <- read.table("/ifs/software/inhouse/GSM/mRNA/reference/versions/1.0/Hg19/EnsemblHg19.RefSeq.refGene.gtf",sep="\t",stringsAsFactors=FALSE);
    names(gtf2) <- c("RefSeq","chr","strand","txStart","txEnd","cdsStart","cdsEnd","exonStart","exonEnd","Symbol")
    gtf2 <- gtf2[!grepl("_",gtf2[,"chr"]),]
    gtf2 <- gtf2[,c(1,10)];
    names(gtf2) <- c("ID","symbol")
    } 
  if (ref_genome=="danRer11"){ 
    gtf2 <- read.table("/ifs/software/inhouse/GSM/mRNA/reference/versions/1.0/danRer11/danRer11_ncbiRefSeq.gtf",sep="\t",stringsAsFactors=FALSE,header=FALSE);
    names(gtf2) <- c("chr","Start","End","exons","exonStart","exonEnd","RefSeq","Symbol")
    gtf2 <- gtf2[!grepl("_",gtf2[,"chr"]),]
    gtf2 <- gtf2[,c("RefSeq","Symbol")]
#    gtf2 <- gtf2[,c(1,10)];
    names(gtf2) <- c("ID","symbol")
    } 
  
  if(ref_genome=="mm10"){
    gtf2 <- read.table("/storage/projects/ngs/analyses/ResearchProjecten/RNAseq/Reference_Sequences/mm10/Mm10_UCSC_RefSeq.gtf",sep="\t",stringsAsFactors=FALSE);
    gtf2 <- gtf2[,c(1,10)];
    names(gtf2) <- c("ID","symbol")
  }

  if(ref_genome=="rn6"){
    gtf2 <- read.table("/storage/projects/ngs/analyses/ResearchProjecten/RNAseq/Reference_Sequences/rn6/UCSC.Rn6.RefSeq.gtf",sep="\t",stringsAsFactors=FALSE);
    gtf2 <- gtf2[,c(1,10)];
    names(gtf2) <- c("ID","symbol")
  }
gtf2
}
