###########################################
## function to plot results from enrichment

plot.enrich <- function(res, s.set=NA, top.x=NA, compress=F){
  
  ## 'res'      -- results from enrichment analysis (including multiple diseases)
  ## 's.set'    -- whether to plot only a subset
  ## 'top.x'    -- whether to only draw the top x components
  ## 'compress' -- whether to compress pathways based on shared genes
  
  ## get all relevant entries (N.B.: this is data.table notation)
  if(!is.na(s.set)){
    tmp <- res[ pheno == s.set]
  }else{
    tmp <- res
  }
  
  ## order by p-value
  tmp <- tmp[ order(p_value)]
  
  ## create new column
  tmp[, log10p := -log10(p_value)]
  
  ## restrict
  if(!is.na(top.x)){
    tmp <- tmp[1:min(c(nrow(tmp), top.x)),]
  }
  
  ## include a compression step
  if(compress == T){
    
    ## get all genes and pathways
    gen     <- unique(unlist(lapply(tmp$intersection, function(c) strsplit(c, ",")[[1]])))
    ## important, path contains pathways ordered by p-value
    path    <- tmp[, term_id]
    ## create matrix
    mat     <- array(data = 0, dim = c(length(gen), length(path)), dimnames = list(gen, path))
    ## fill the matrix
    for(j in 1:ncol(mat)){
      ## get the genes belonging to the relevant pathway
      p.gen         <- strsplit(tmp$intersection[which(tmp$term_id == colnames(mat)[j])], ",")[[1]]
      ## replace ones for matching genes
      mat[p.gen, j] <- 1
    }
    ## store the number of genes belonging to each pathway
    g.path  <- colSums(mat)
    
    ## empty list of pathways to include
    in.path <- c()
    
    ## do as long as any entry remains
    while(nrow(mat) > 0 & ncol(mat) > 0){
      
      ## add the strongest remaining pathway
      in.path <- c(in.path, path[1])
      
      ## drop all genes belonging to this pathway
      mat     <- mat[-which(mat[, path[1]] == 1),, drop=F]
      ## drop empty and depleted pathways
      ii      <- colSums(mat)
      ii      <- ii/g.path
      mat     <- mat[, which(ii > .5), drop=F]
      ## reduce the list of pathways to look at
      path    <- path[which(path %in% colnames(mat))]
      g.path  <- g.path[which(names(g.path) %in% colnames(mat))]
      
    }
    
    ## report only those back
    tmp <- tmp[ term_id %in% in.path ]
    
  }
  
  ## colour gradient for fold enrichment
  col.vec <- colorRampPalette(c("white", "#fd5602"))(20)
  
  ## empty plot
  plot(c(0, max(tmp$log10p)), c(.5, nrow(tmp)+.5), xlab=expression(-log[10]("p-value")),
       type="n", ylab="", xaxt="n", yaxt="n", ylim=rev(c(.5, nrow(tmp)+.5)))
  ## add axis
  axis(1, lwd=.5)
  
  ## add rectangles
  print(tmp$fc)
  rect(0, 1:nrow(tmp)-.4, tmp$log10p, 1:nrow(tmp)+.4, border="grey50", lwd=.2,
       col=ifelse(tmp$fc >= 20, col.vec[20], col.vec[ceiling(tmp$fc)]))
  
  ## add labels
  # text(0, 1:nrow(tmp), pos=4, labels = paste0(stringr::str_to_sentence(tmp$term_name), " (", tmp$term_id, ")"), cex=.4,
  #      font=2, offset=.05)
  ## from here https://stackoverflow.com/questions/2351744/insert-line-breaks-in-long-string-word-wrap
  text(0, 1:nrow(tmp), pos=4, labels = sapply(tmp$term_name, function(x) paste(strwrap(stringr::str_to_sentence(x), 70), collapse="\n")), cex=.4,
       font=2, offset=.05)
  
  ## add final frame
  pm <- par("usr")
  rect(pm[1], pm[3], pm[2], pm[4], lwd=.5)
  
  return(tmp)
}