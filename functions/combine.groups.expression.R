combine.groups.expression <- function(expr_df,groups,minreads){
 groups.total <- data.frame(stringsAsFactors=FALSE)
 groups.lvl   <- levels(groups)
 for(lvl in groups.lvl){
   totalGenes.group    <- expr_df[,c(1,which(grepl(lvl,colnames(expr_df))))]
   if(lvl == groups.lvl[1]) {
     if(ncol(totalGenes.group)==2){groups.total <- totalGenes.group[totalGenes.group[,-1]>=minreads,] } else {
       groups.total <- totalGenes.group[rowSums(totalGenes.group[,-1]>=minreads)==ncol(totalGenes.group[,-1]),]}
     groups.total <- groups.total[order(groups.total[,1],decreasing=TRUE),]} else {
     if(ncol(totalGenes.group)==2){groups.exp <- totalGenes.group[totalGenes.group[,-1]>=minreads,]} 
     else {groups.exp <- totalGenes.group[rowSums(totalGenes.group[,-1]>=minreads)==ncol(totalGenes.group[,-1]),]}
     groups.exp                <- groups.exp[order(groups.exp[,1],decreasing=TRUE),]
     groups.total_unq.all      <- unique(c(groups.total[,"ID"],groups.exp[,"ID"]))
     groups.total.com          <- groups.total_unq.all %in% groups.total[,"ID"]
     groups.exp.com            <- groups.total_unq.all %in% groups.exp[,"ID"]
     groups.total_T            <- data.frame( groups.total.com, groups.exp.com ,stringsAsFactors=FALSE)
     row.names(groups.total_T) <- groups.total_unq.all
     head(groups.total_T)
     groups.total.unq            <- groups.total_unq.all [groups.total_T[,1] & groups.total_T[,2] ]
     groups.total.overlap        <- groups.total[groups.total[,"ID"]%in%groups.total.unq,]
     groups.exp.overlap          <- groups.exp[groups.exp[,"ID"]%in%groups.total.unq,]
     groups.total.Overlap.total  <- data.frame(ID=groups.total.unq,groups.total.overlap[,-1], groups.exp.overlap[,-1],stringsAsFactors=FALSE)
     names(groups.total.Overlap.total) <- c("ID",(names(groups.total.overlap)[-1]),names(groups.exp.overlap)[-1])
    
     groups.total_unq   <- groups.total_unq.all [groups.total_T[,1]==TRUE & groups.total_T[,2]==FALSE ]
     groups.exp_unq     <- groups.total_unq.all [groups.total_T[,1]==FALSE & groups.total_T[,2]==TRUE ]
    
     groups.total.combi <- groups.total[groups.total[,1]%in%groups.total_unq,]
     groups.exp.combi   <- data.frame( minreads , stringsAsFactors=FALSE)
     while(ncol(groups.exp.combi)<(ncol(groups.exp.overlap)-1)) {
       if (ncol(groups.exp.combi)==ncol(groups.exp.overlap)-1) break 
       else groups.exp.combi <- data.frame(groups.exp.combi, minreads , stringsAsFactors=FALSE)}
     groups.total.combi.total <- data.frame(ID=groups.total_unq,groups.total.combi[,-1],groups.exp.combi,stringsAsFactors=FALSE)
     names(groups.total.combi.total) <- c("ID",(names(groups.total.combi)[-1]),names(groups.exp.overlap)[-1])
     
     if(length(groups.exp_unq)==0){  groups.total <- rbind(groups.total.Overlap.total,groups.total.combi.total)
     groups.total <- groups.total[order(groups.total[,1],decreasing=TRUE),]} else {groups.exp.combi1   <- groups.exp[groups.exp[,1]%in%groups.exp_unq,]
           groups.total.combi1 <- data.frame( minreads , stringsAsFactors=FALSE)
      if(ncol(groups.total.overlap)-1 != 1 ) {
	  while(ncol(groups.total.combi1)<(ncol(groups.total.overlap)-1)) { 
        if (ncol(groups.total.combi1)==ncol(groups.total.overlap)-1) break else groups.total.combi1 <- data.frame(groups.total.combi1, minreads , stringsAsFactors=FALSE)}}
      names(groups.total.combi1) <- names(groups.total.overlap)[-1]
      groups.exp.combi.total     <- data.frame(ID=groups.exp_unq,groups.total.combi1,groups.exp.combi1[,-1],stringsAsFactors=FALSE)
      names(groups.exp.combi.total) <- c("ID",names(groups.total.combi1),(names(groups.exp.combi1)[-1]))
    
     groups.total <- rbind(groups.total.Overlap.total,groups.exp.combi.total,groups.total.combi.total)
     groups.total <- groups.total[order(groups.total[,1],decreasing=TRUE),]
     }
   }
 }
groups.total 
}
