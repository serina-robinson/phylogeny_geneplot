get_dna_segs<-function(nb, pfam_query){
  DL <- list()
  NL <- nb

  for(i in 1:length(nb)) {
    name=nb[[i]]$query
    start <- vector(length=nrow(nb[[i]]))
    end <- vector(length=nrow(nb[[i]]))
    direct <- nb[[i]]$dir
    strand <- ifelse(direct=="+", 1, -1)
    # 
    # STRAND <-strand
    # STRAND
    # if(pfam_query %in% nb[[i]]$pfam_id1){
    #   if(na.omit(nb[[i]]$dir[nb[[i]]$pfam_id1 == pfam_query])[1] == "-"){
    #     strand<-strand * (-1)
    #     ord<-sort(1:nrow(nb[[i]]), decreasing=T)
    # 
    #     for(k in 1:length(ord)){
    #       NL[[i]][k,]<-nb[[i]][(ord[k]),]
    #       STRAND[k]<-strand[ord[k]]
    #     }
    #   }
    # }
    
    nb[[i]]<-NL[[i]]
   # strand<-STRAND
    pfam<-nb[[i]]$pfam_id1
    
    for(j in 1:nrow(nb[[i]])){
     # if(j==1){
     if (nb[[i]]$dir[j] == "+") {
        start[j] <- nb[[i]]$start[j]
        end[j] <- nb[[i]]$end[j]
    }
      else if (nb[[i]]$dir[j] == "-") {
        start[j] <- nb[[i]]$end[j]
        end[j] <- nb[[i]]$start[j]
      }
    }
    print(strand)
    
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


