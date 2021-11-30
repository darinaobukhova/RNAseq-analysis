# Function to divide genome into single base start-positions (single base starting position alignment read prevents reads overlapping several transcripts, several 100kb region)
getSequenceLengths <- function(bmfl) {
  # read the BAM header and obtain the lengths per chromosome name
  txt <- unlist(scanBamHeader(bmfl))
  inm <- grep("SN:",txt )
  
  iln <- inm + 1
  
  nm  <- gsub( "SN:", "", txt[ inm ] )
  #if (length(grep("chr",nm))==0) { nm <- paste("chr",nm,sep="")} 
  ln  <- as.integer( gsub( "LN:","",txt[ iln ] ) )
  
  names(ln) <- nm
  ln <- ln[names(ln)%in%unique(exons.allchr[,"chr"])]
   ln
}

getCounts <- function( bmfl2, targets=NULL, param=NULL, which2=NULL ) {
  if( is.null(targets) ) stop("the targets to count is NULL")
  if( is.null(param) ) stop("param is NULL")
  if( is.null(which2) ) stop("which is NULL")
  # get the start positions of reads on the genome
  bam   <- scanBam( bmfl2, param=param )
  # convert the start positions to single base regions
  ilist <- lapply( bam, function(x) {IRanges(start=x[[1]][!is.na(x[[1]])], width=1) } )

  irng  <- IRangesList()
  for( i in 1:length(ilist) ) { irng[[i]] <- ilist[[i]] }
  names(irng) <- names(which2)
  
  # get the overlaps between the targets and the iranges
  lcnt <- lapply( 1:length(targets), function(i, tgt, rds ){
    chr <- names(tgt)[i]
    countOverlaps( tgt[[chr]], rds[[chr]] )
  }, tgt=targets, rds=irng )
  
  # return the counts vector
  unlist(lcnt)
}
