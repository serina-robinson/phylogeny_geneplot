# Install packages
pacman::p_load("tidyverse", "genoPlotR","readxl", "janitor", "RColorBrewer","rentrez","ggthemes","ggplot2","Biostrings","data.table","ape","phangorn")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/phylogeny_geneplot/")

# Read in the RODEO output (main_co_occur)
co_occur_xl <- read_excel("data/triuret_main_co_occur.xlsx") %>%
  janitor::clean_names()

# Pull organism names
orgs <- co_occur_xl %>%
  group_by(query) %>%
  dplyr::select(genus_species, query) %>%
  mutate(query_fix = gsub("\\.", "_", query)) %>%
  slice(1) 
orgs

# Pull the accession number of the query sequences
accs <- co_occur_xl %>%
  dplyr::select(query) %>%
  distinct() %>%
  pull()

# Pull the sequences
# sqs <- entrez_fetch(id = accs, db = "protein", rettype = "fasta",
#                    api_key="826357f5ff17c7ec62e583909071e94f9d08")

# Write the sequences to file
# write(sqs, "data/124_triuret_hydrolase_seqs.fasta")

# Read in fasta
aa <- readAAStringSet("data/124_triuret_hydrolase_seqs.fasta")
summary(width(aa)) # width distribution looks good

# Align sequences # Only needs to run once
# aa_al <- AlignSeqs(aa)
# aa_adj <- AdjustAlignment(aa_al)
# aa_stag <- StaggerAlignment(aa_adj)
# BrowseSeqs(aa_stag)
# writeXStringSet(aa_stag, "data/124_triuret_hydrolase_seqs_aligned_stag.fasta")

# Phylogeny using FastTree JTT + CAT model
# FastTree < alignment_file > tree_file 
# File is written to 124_triuret_hydrolase_seqs_aligned_stag.fasta

# Read in newick tree
my_tre <- read.table("data/124_triuret_hydrolase_tree.nwk", header = F) # estimated using FastTree 
my_nwk<-newick2phylog(x.tre = my_tre[,1])

# Group by query
nb <- co_occur_xl %>%
  split(.$query)
nb[[1]]

# Get dna seg objects
# Are all colored gray by default
source("lib/get_dna_segs_v3.r")
raw_segs <- get_dna_segs(nb)
raw_segs[[1]]

# Color the dna segs
source("lib/color_segs.r")
colrs <- fread("data/dna_seg_colors.txt", sep = " ", header = F, data.table = F)
head(colrs)
pfams <- colrs$V1
colrs <- colrs$V2
segs <- color_segs(raw_segs, pfams, colrs)


 # Look at the segs
pdf("output/seg_plot.pdf", width = 15, height = 50)
plot_gene_map(segs)
dev.off()


tonam <- unlist(lapply(1:length(segs),function(x){ gsub("\\.","_",segs[[x]]$name[1]) }))
 # 124 names should match the tip labels

# Reorder segment labels so the match the tree leaves
for(i in 1:length(tonam)){
  tmp <- grep(tonam[i], names(my_nwk$leaves))
  names(segs)[i] <- names(my_nwk$leaves)[tmp]
}
names(segs)

source("lib/reorder_list.r")

segs2 <- reorder_list(segs, names(my_nwk$leaves))
numseqs <- length(segs2)
names(segs2)

# Reorder organism names to match leaves
neworgs <-rep(NA, length(segs))
for(i in 1:length(segs)){
  ind <- grep(names(segs)[i], names(my_nwk$leaves))
  print(ind)
  neworgs[ind]<- orgs$genus_species[i]
}
neworgs


annots <- lapply(1:length(my_nwk$leaves), function(x){
  mid <- middle(segs2[[x]])
  annot <- annotation(x1 = mid, text = c(rep("", nrow(segs2[[x]])-1), neworgs[x]))
})


# Plot the tree
pdf("output/triuret_gene_neighborhood_tree_names_fixed.pdf", width = 15, height = 50)
plot_gene_map(segs2, tree = my_nwk, tree_width = 5, 
              tree_branch_labels_cex = 0,
              annotations = annots)
dev.off() 


# Write a legend
pdf(paste0("output/triuret_color_legend.pdf"),width=10,height=10)
plot.new()
legend("bottom", ncol = 1, legend = pfams, fill = colrs,
          border=FALSE, bty="n")
dev.off()
