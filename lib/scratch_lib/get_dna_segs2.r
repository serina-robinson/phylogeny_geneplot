get_dna_segs<-function(my_list){
## Function that reads in a list of DNA sequences
## Returns a DNA seg object
  
  for(i in 1:length(my_list)) {
    tmp <- my_list[[i]]
    name <- unlist(lapply(1:length(tmp), function(x) { strsplit(names(tmp), "\\.1")[[x]][2] }))
    start <- vector(length = length(tmp))
    end <- vector(length = length(tmp))
    for(j in 1:length(tmp)){
      if(j == 1){
        start[j]<-10
        end[j]<-width(tmp)[j]
      }
      else{
        start[j]<-end[j-1] + 10
        end[j]<-(width(tmp)[j] + start[j] + 1)
      }
    }
    strand <- rep(1,times = length(tmp))
    df1 <- data.frame(name = name,
                    start = start,
                    end = end,
                    strand = strand)
    dna_seg1 <- dna_seg(df1)
    num <- length(dna_seg1$col)

    dna_seg1$col <- "gray"
    dna_seg1$fill <-"gray"
    dna_seg1$gene_type <- "headless_arrows"
  }
  return(dna_seg1)
}


