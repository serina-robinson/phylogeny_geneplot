get_dna_segs<-function(nb, pfam_query){
  DL <- list()
  NL <- nb

  for(i in 1:length(nb)) {
    name <- nb[[i]]$query
    start <- vector(length=nrow(nb[[i]]))
    end <- vector(length=nrow(nb[[i]]))
    direct <- nb[[i]]$dir
    strand <- ifelse(direct=="+", 1, -1)
    nb[[i]]<-NL[[i]]
    pfam<-nb[[i]]$pfam_id1
    
    for(j in 1:nrow(nb[[i]])){
     if (nb[[i]]$dir[j] == "+") {
        start[j] <- nb[[i]]$start[j]
        end[j] <- nb[[i]]$end[j]
    }
      else if (nb[[i]]$dir[j] == "-") {
        start[j] <- nb[[i]]$end[j]
        end[j] <- nb[[i]]$start[j]
      }
    }
    
    df1<-data.frame(name=name,
                    start=start,
                    end=end,
                    strand=strand,
                    pfa=pfam)
    dna_seg1<-dna_seg(df1)
    num<-length(dna_seg1$col)
    dna_seg1$col<-"gray"
    dna_seg1$fill<-"gray"
    dna_seg1$gene_type<-"arrows"
    DL[[i]]<-dna_seg1
  }
  return(DL)
}


