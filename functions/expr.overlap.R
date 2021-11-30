expr.overlap <- function(expr_df,input,wkd.result,groups){
  
  expr_list <- list()
  expr.total_unq    <- NULL
  if(length(levels(groups))<=5){
  for(lvl in levels(groups)){
    expr_df_lvl      <- expr_df[,which(grepl(lvl,groups))+1]>5
    expr_df_lvl      <- data.frame(ID=expr_df[,1],expr_df_lvl,stringsAsFactors=FALSE)
    expr_df_lvl_exp  <- expr_df_lvl[rowSums(expr_df_lvl[,-1])==ncol(expr_df_lvl[,-1]),]
    expr_list[[lvl]] <- expr_df_lvl_exp 
    expr.total_unq   <- unique(c(expr.total_unq,expr_list[[lvl]][,"ID"]))
  }
  
  expr_df_expr_overlap <- data.frame(ID=expr.total_unq,stringsAsFactors=FALSE)
  for(lvl in levels(groups)){
    expr_df_overlap      <- expr.total_unq  %in% expr_list[[lvl]][,"ID"]
    colflsnames <- colnames(expr_df_expr_overlap)
    expr_df_expr_overlap <- data.frame( expr_df_expr_overlap, expr_df_overlap,stringsAsFactors=FALSE)
    names(expr_df_expr_overlap) <- c(colflsnames,lvl)
  }
  
  pdf(paste(wkd.result,paste(input,"venn_diagram_overlap.pdf",sep="_"),sep="/"))  
  reads.ven           <- venn(expr_df_expr_overlap[,-1])
  dev.off()
 }
}
