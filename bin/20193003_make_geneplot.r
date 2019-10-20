# Install packages
pacman::p_load("tidyverse", "genoPlotR", "ggsci", "readxl", "janitor", "RColorBrewer","rentrez","ggthemes","ggplot2","Biostrings","data.table","ape","phangorn")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/phylogeny_geneplot/")

# Read in the Xanthomonas campestris output
# co_occur_xan <- read_csv("../../github/beta-lactone-NP-discovery/r-scripts/ole_genome_blast/data/rodeo_output/20171012_493seqs_RODEO_output/results_folder/493_OleC_neighborhoods.csv") %>%
co_occur_xan <- read_csv("../../github/beta-lactone-NP-discovery/r-scripts/nocardia_1000_BLAST/data/RODEO_output (19)/NltC_output/main_co_occur.csv") %>%
  janitor::clean_names() %>% 
  dplyr::filter(grepl("NZ_LHUK010004", nucleotide_acc))
co_occur_xan

# Read in the RODEO output (main_co_occur)
co_occur_xl <- read_csv("../../github/jgi_oleA/data/43_JGI_Ole_genes.csv") %>%
  janitor::clean_names() %>%
  dplyr::filter(grepl("Mobilicoccus massiliensis|Granulosicoccus antarcticus|Lysobacter tolerans|Actinoplanes atraurantiacus", genus_species)) %>%
  bind_rows(., co_occur_xan) %>%
  dplyr::filter(!query %in% c("WP_097326010", "WP_097326008")) %>%
  group_by(query) %>%
  dplyr::slice(4:13)
co_occur_xl$genus_species[grepl("Lysobacter", co_occur_xl$genus_species)] <- gsub( "Lysobacter", "Luteimonas", co_occur_xl$genus_species[grepl("Lysobacter", co_occur_xl$genus_species)])

# Pull organism names
orgs <- co_occur_xl %>%
  group_by(query) %>%
  dplyr::select(genus_species, query) %>%
  mutate(query_fix = gsub("\\.", "_", query)) %>%
  dplyr::slice(1)
orgs

# Pull the accession number of the query sequences
accs <- co_occur_xl %>%
  dplyr::select(query) %>%
  distinct() %>%
  pull()


accs
orgs

# Group by query
nb <- co_occur_xl %>%
  split(.$query) 

# Get the most frequent PFAM IDs
source("lib/make_nb_barplot.R")
colnames(co_occur_xl)
pfams<-make_nb_barplot(co_occur_xl)
pfams <- rbind(pfams[1:7,], pfams[pfams$pfam == "PF00535",])

# Get dna seg objects
# Are all colored gray by default
source("lib/get_dna_segs_v3.r")
raw_segs <- get_dna_segs(nb)
head(raw_segs[[1]])

# Color the dna segs
pfam_desc$pfam
pal
colrs <- data.frame(cbind(as.character(pfam_desc$pfam), pal), stringsAsFactors = F)
colrs

colnames(colrs) <- c("V1", "V2")
head(colrs)
pfams2 <- colrs$V1
colrs2 <- colrs$V2
source("lib/color_segs.r")
segs <- color_segs(raw_segs, pfams2, colrs2)
names(segs) <- orgs$genus_species

source("lib/reorder_list.r")
segs <- reorder_list(segs, sort(names(segs), decreasing = T))

# Annotations
# annots <- lapply(1:length(segs), function(x){
#   mid <- middle(segs[[x]])
#   #segs[[x]]$name[length(segs[[x]]$name)]<-nms[x] ##Only uncomment if you pull from rentrez
#   annot <- genoPlotR::annotation(x1=mid, text=c(rep("",round(length(segs[[x]]$name)/2)),  
#                                                 # paste0(segs[[x]]$name[length(segs[[x]]$name)], " ", orgs$genus_species[x]),
#                                                 orgs$genus_species[x],
#                                                 rep("", round(length(segs[[x]]$name)/2-1))))
#   
# })
# annots

# Plot the segements
pdf("output/seg_plot.pdf", width = 8, height = 3)
plot_gene_map(segs) #annotations= annots, 
              #annotation_height = 3, annotation_cex = 1.5)
dev.off()

pfams$desc1[is.na(pfams$desc1)] <- pfams$desc2[is.na(pfams$desc1)]
txt <- paste0(pfams$desc1, " (", pfams$pfam, ")")
txt
txt2 <- gsub("AMP-binding enzyme", "OleC", txt)
txt2 <- gsub("OleC (OleBC)", "OleBC fusion (PF00561 + PF00501)", txt2, fixed = T)
txt2 <- gsub("(DUF4156)", "DUF4156", txt2, fixed = T)
txt2 <- gsub("3-beta hydroxysteroid dehydrogenase/isomerase family", "OleD", txt2)
txt2 <- gsub("alpha/beta hydrolase fold", "OleB", txt2)
txt2
# txt2  <- gsub("s group 1", "", txt2)
# txt2 <- gsub(", catalytic domain", "", txt2)
txt2 <- gsub("3-Oxoacyl-[acyl-carrier-protein (ACP)] synthase III C terminal", "OleA", txt2, fixed = T)
segs[[5]]

pdf(paste0("output/Megan_paper_geneplot_legend.pdf"), width = 12, height = 6)
plot.new()
legend("bottom", legend = txt2, fill=pal, ncol = 2)#ncol = 1)
      # border = FALSE, bty = "n", title = "Enzyme")
dev.off()
