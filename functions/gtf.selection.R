gtf.selection <- function(ref_genome,ref_gtf){
  
if (ref_genome=="Hg38"){ 
  if(ref_gtf=="RefSeq"){ 
    gtf.whole <- read.table("/ifs/software/inhouse/GSM/mRNA/reference/versions/1.0/Hg38_RefSeq_refGene_whole_exonID.gtf",header=TRUE,sep="\t",stringsAsFactors=FALSE)
  } 
  else {gtf.whole <- read.table("/storage/projects/ngs/analyses/ResearchProjecten/RNAseq/Reference_Sequences/Hg38/Hg38_UCSC_Gencode24_exon.gtf",header=TRUE,sep="\t",stringsAsFactors=FALSE)
  }
 } 

if (ref_genome=="Hg19"){ 
  gtf.whole <- read.table("/ifs/software/inhouse/GSM/mRNA/reference/versions/1.0/Hg19/Hg19_RefSeq_refGene_whole_exonID.gtf",header=TRUE,sep="\t",stringsAsFactors=FALSE); 
 } 
  
if(ref_genome=="mm10"){
  gtf.whole <- read.table("/storage/projects/ngs/analyses/ResearchProjecten/RNAseq/Reference_Sequences/mm10/Mm10_RefSeq_refGene_exonID.gtf",sep="\t",header=TRUE,stringsAsFactors=FALSE);
  }
