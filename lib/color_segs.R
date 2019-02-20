color_segs<-function(raw_segs, pfams, pal){

  for(i in 1:length(raw_segs)){
    for(j in 1:length(pfams)){
      ind<-grep(pfams[j],x = raw_segs[[i]]$pfa)
      raw_segs[[i]]$col[ind]<-pal[j]
      raw_segs[[i]]$fill[ind]<-pal[j]
    }
  }  
  
  return(raw_segs)
}

