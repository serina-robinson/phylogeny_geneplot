make_nb_barplot<-function(noc){

  # Trim the hmm column
  # noct<-noc[,c(1:2,4:(ncol(noc)-4))]
  noct <- noc
  
  #Group by the OleC neighborhoods
  nb<-list()
  Cacc<-unique(noct$query)
  for(i in 1:length(Cacc)){
    nb[[i]]<-noct[noct$query==Cacc[i],]
  }
  
  #Look at overall PFAM enrichment
  pfam<-noc$pfam_id1
  pfamtab<-data.frame(table(pfam))

  pfamsort<-pfamtab[order(pfamtab$Freq,decreasing=T),]
  pfam10<-pfamsort
  #[1:8,] 

  #Get PFAM ids
  desc1 <- unlist(lapply(1:length(pfam10$pfam), function(x){ noc$description[grep(pfam10$pfam[x],noc$pfam_id1)][1] }))
  desc2 <- unlist(lapply(1:length(pfam10$pfam), function(x){ noc$description[grep(pfam10$pfam[x],noc$pfam_id2)][1] }))
  pfam_desc<-data.frame(as.character(pfam10$pfam), pfam10$Freq, desc1, desc2, stringsAsFactors=F)
  # head(pfam_desc)
  # pfam_desc$desc<-as.character(pfam_desc$desc)
  # pfam_desc$desc[pfam_desc$pfam=="NO PFAM MATCH"]<-"hypothetical protein"

  colnames(pfam_desc)[1:2]<-c("pfam","freq")
  #return(pfam_desc)

  pfam_desc$perc<-(pfam_desc$freq)/length(Cacc)
  # pfam_desc$desc[pfam_desc$pfam=="NO PFAM MATCH"]<-"hypothetical protein"
  # pfam_desc$desc[pfam_desc$desc==""]<-"hypothetical protein"
  
  
  pfam_desc$pfam<-factor(pfam_desc$pfam,levels=pfam_desc$pfam)
  # source("src/theme.publication.r")
  # #pdf("output/cysteine_hydrolase_neighborhoods.pdf",width=25,height=15)
  # pal<-palette(colorRampPalette(colors=brewer.pal(12,"Paired"))(12))
  # gr<-ggplot(data=pfam_desc,aes(y=perc,x=pfam))+
  #   geom_bar(stat="identity",fill=pal) +
  #   theme.publication(pfam) +
  #   scale_x_discrete(name ="PFAM ID") +
  #   scale_y_continuous(name = "Gene frequency") +
  #   geom_text(aes(label=round(perc,2)),size=12,vjust=-0.25)
  # #dev.off()
  # 
  # pdf(paste0("output/",length(Cacc),"_RODEO_barplot.pdf"),width=25,height=15)
  # ggplot(data=pfam_desc,aes(y=perc,x=pfam))+
  #   geom_bar(stat="identity",fill=pal) +
  #   theme.publication(pfam) +
  #   scale_x_discrete(name ="PFAM ID") +
  #   scale_y_continuous(name = "Gene frequency") +
  #   geom_text(aes(label=round(perc,2)),size=12,vjust=-0.25)
  # dev.off()
  # 
  # pdf(paste0("output/",length(Cacc),"_RODEO_legend.pdf"),width=10,height=10)
  # plot.new()
  # legend("bottom",legend=pfam_desc$desc,fill=pal,
  #        border=FALSE, bty="n", title = "Enzyme")
  # dev.off()
  # 
  return(pfam_desc)
}