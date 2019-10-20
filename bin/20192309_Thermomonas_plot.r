# Actinoplanes 
# Micromonospora peucetia is the same as Granulosicoccus
# Mobilicoccus massiliensis 
# Arthrobacter globiformis
# 

# Install packages
pacman::p_load("tidyverse", "genoPlotR", "ggsci", "readxl", "janitor", "RColorBrewer","rentrez","ggthemes","ggplot2","Biostrings","data.table","ape","phangorn")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/phylogeny_geneplot/")

# Read in the Xanthomonas campestris output
# co_occur_xan <- read_csv("../../github/beta-lactone-NP-discovery/r-scripts/ole_genome_blast/data/rodeo_output/20171012_493seqs_RODEO_output/results_folder/493_OleC_neighborhoods.csv") %>%
# co_occur_xan <- read_csv("../../github/beta-lactone-NP-discovery/r-scripts/nocardia_1000_BLAST/data/RODEO_output (19)/NltC_output/main_co_occur.csv") %>%
co_occur_1 <- read_csv("data/920_OleA/main_co_occur.csv") %>%
  janitor::clean_names() %>% 
  # dplyr::filter(grepl("NZ_LHUK010004", nucleotide_acc)) %>%
  dplyr::filter(grepl(paste0(c("Mobilicoccus massiliensis", "Granulosicoccus antarcticus", "Lysobacter tolerans", "Chromatocurvus", "Arthrobacter globiformis",
                               "Actinoplanes atraurantiacus", #"Silanimonas lenta", "Arenimonas oryziterrae", "Thermomonas haemolytica", "Pseudoxanthomonas",
                               "Micromonospora peucetia", "Kytococcus sedentarious", "Xanthomonas translucens"), collapse = "|"), genus_species))
unique(co_occur_1$genus_species)

co_occur_2 <- read_csv("data/46_JGI_OleA/output/main_co_occur.csv") %>%
  janitor::clean_names() %>% 
  # dplyr::filter(grepl("NZ_LHUK010004", nucleotide_acc)) %>%
  dplyr::filter(grepl(paste0(c("Mobilicoccus massiliensis", "Granulosicoccus antarcticus", "Lysobacter tolerans", "Chromatocurvus", "Arthrobacter globiformis",
                               "Actinoplanes atraurantiacus", "Silanimonas lenta", "Arenimonas oryziterrae", "Thermomonas haemolytica", "Pseudoxanthomonas",
                               "Micromonospora peucetia", "Kytococcus", "Xanthomonas translucens"), collapse = "|"), genus_species))  

unique(co_occur_2$genus_species)

comb <- unique(unique(co_occur_2$genus_species), unique(co_occur_2$genus_species))

# Read in the RODEO output (main_co_occur)
head(co_occur_2)
co_occur_xl <- co_occur_1 %>%
  bind_rows(., co_occur_2) %>%
  group_by(query) 
# dplyr::slice(4:13)
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

# Group by query
nb <- co_occur_xl %>%
  dplyr::mutate(pfam_id1 = case_when(pfam_id1 == "PF08545" ~ "PF08541",
                                     pfam_id1 == "PF01370" ~ "PF01073",
                                     TRUE ~ pfam_id1)) %>%
  split(.$query) 

# Fix the PFAM misannotations 


# Get the most frequent PFAM IDs
source("lib/make_nb_barplot.R")
colnames(co_occur_xl)

pfams <- make_nb_barplot(co_occur_xl)
pfams[1:3,]
pfams$desc1[1:3] <- c("OleC", "OleA",  "OleD")
head(pfams)
# pfams <- rbind(pfams[1:7,], pfams[pfams$pfam == "PF00535",], pfams[3:4,], pfams[pfams$pfam == "PF00067",], pfams[5:7,])
pfams <- rbind(pfams[1:7,]) #pfams[pfams$pfam == "PF00535",], pfams[pfams$pfam == "PF00067",])
pfams
# Get dna seg objects
# Are all colored gray by default
source("lib/get_dna_segs_v3.r")
raw_segs <- get_dna_segs(nb)

# raw_segs[[1]]$pfa <- gsub("PF00501", "OleBC fusion", raw_segs[[1]]$pfa)
# raw_segs[[4]]$pfa <- gsub("PF00501", "OleBC fusion", raw_segs[[1]]$pfa)

# raw_segs[[9]]$pfa <- gsub("PF08545", "PF08541", raw_segs[[9]]$pfa) # check these
# raw_segs[[9]]$pfa <- gsub("PF01370", "PF01073", raw_segs[[9]]$pfa) # check these
# raw_segs[[9]]$pfa <- gsub("PF08545", "PF08541", raw_segs[[9]]$pfa) # check these

# Color the DNA segs
# colrs <- fread("data/dna_seg_colors.txt", sep = " ", header = F, data.table = F)

# pal <- c(colorRampPalette(colors = brewer.pal(12, "Paired"))(12), "black")
# pal <- c(colorRampPalette(colors = brewer.pal(12, "Paired"))(12), "black")
pal <- c(colorRampPalette(colors = brewer.pal(8, "Set1"))(8), "black")
length(pfams)
pfams
# pal <- hue_pal()(nrow(pfams))
#pal[6] <- "cornflowerblue"
# pal[4] <- "orange"
colrs <- data.frame(cbind(as.character(pfams$pfam), pal), stringsAsFactors = F)
colrs
colnames(colrs) <- c("V1", "V2")

head(colrs)
pfams2 <- colrs$V1
colrs2 <- colrs$V2
pfams2
pfams3 <- pfams2[1:7]
colrs3 <- colrs2[1:7]
desc3 <- pfams$desc1[1:7]
desc3

source("lib/color_segs.r")
segs <- color_segs(raw_segs, pfams3, colrs3)
segs
segs[[1]]

nams <- unique(unlist(lapply(1:length(segs), function(x) segs[[x]]$name)))

unique(co_occur_xl$genus_species)

dtf <- data.frame(cbind(co_occur_xl$query, co_occur_xl$genus_species), stringsAsFactors = F)
colnames(dtf) <- c("acc", "taxon")
dtf_ord <- dtf[order(dtf$acc),]
tail(dtf_ord)
nams
genun <- as.character(unique(dtf_ord$taxon))
genun

names(segs) <- genun #unique(co_occur_xl$genus_species)
#specnams <- c("Mobilicoccus massiliensis SIT2", "Xanthomonas campestris pv. campestris str. ATCC 33913",
#              "Luteimonas tolerans UM1", "Granulosicoccus antarcticus IMCC3135", "Actinoplanes atraurantiacus CGMCC 4.6857")
# names(segs) <- specnams

source("lib/reorder_list.r")
length(names(segs))
names(segs)

list_reorder <- names(segs)[c(1, 12, 4, 5, 11, 2, 3, 6, 7, 8, 9, 10)]
list_reorder

segs <- reorder_list(segs, list_reorder)
# segs
# segs[[2]] <- segs[[2]][1:8,]
# segs[[3]] <- segs[[3]][2:7,]
# segs[[4]] <- segs[[4]][5:10,]
# segs[[5]] <- segs[[5]][5:10,]
# segs

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
head(segs)
# Plot the segements
pdf("output/20192309_Thermomonas.pdf", width = 12, height = 5)
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
pal
txt2

pfams3
colrs3

pdf(paste0("output/20192309_Thermomonas_legend.pdf"), width = 20, height = 10)
plot.new()
legend("bottom", legend = paste0(pfams3, " (", desc3, ")"), fill = colrs3, ncol = 2, bty = "n") #ncol = 1)
# border = FALSE, bty = "n", title = "Enzyme")
dev.off()
