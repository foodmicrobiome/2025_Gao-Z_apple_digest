########################
### Apple Microbiome ###
########################

# rblauste@umd.edu #
# last edited: 03/01/2025 #

# load packages
library(SCRuB)
library(ggplot2)
library(vegan)
library(reshape2)
library(gridExtra)
library(stringr)


###################################################################
### Analysis prep: Sequenced read counts to relative abundances ###
###################################################################

###
### load tables & process read counts with SCRuB
###

## metadata for SCRuB to remove negative control counts
meta_scrub = read.csv("sample_data/metadata/metadata_scrub.csv", h=T, row.names=1)
meta_scrub = meta_scrub[order(rownames(meta_scrub)),]
rownames(meta_scrub) = str_replace_all(rownames(meta_scrub), "-", ".")
head(meta_scrub)

## phylum (level 2)
# set table
taxa = read.csv("sample_data/taxa_tables/level-2-counts_phylum.csv", h=T, row.names=1)
taxa = taxa[,order(colnames(taxa))]
colnames(taxa) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(taxa), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_library_prep"))
taxa_phy_s = as.data.frame(t(scr_out$decontaminated_samples))

## class (level 3)
# set table
taxa = read.csv("sample_data/taxa_tables/level-3-counts_class.csv", h=T, row.names=1)
taxa = taxa[,order(colnames(taxa))]
colnames(taxa) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(taxa), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_library_prep"))
taxa_class_s = as.data.frame(t(scr_out$decontaminated_samples))

## order (level 4)
# set table
taxa = read.csv("sample_data/taxa_tables/level-4-counts_order.csv", h=T, row.names=1)
taxa = taxa[,order(colnames(taxa))]
colnames(taxa) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(taxa), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_library_prep"))
taxa_order_s = as.data.frame(t(scr_out$decontaminated_samples))

## family (level 5)
# set table
taxa = read.csv("sample_data/taxa_tables/level-5-counts_family.csv", h=T, row.names=1)
taxa = taxa[,order(colnames(taxa))]
colnames(taxa) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(taxa), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_library_prep"))
taxa_fam_s = as.data.frame(t(scr_out$decontaminated_samples))

## genus (level 6)
# set table
taxa = read.csv("sample_data/taxa_tables/level-6-counts_genus.csv", h=T, row.names=1)
taxa = taxa[,order(colnames(taxa))]
colnames(taxa) == rownames(meta_scrub)
#scrub
scr_out <- SCRuB(t(taxa), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_library_prep"))
taxa_gen_s = as.data.frame(t(scr_out$decontaminated_samples))

## summary stats (genus counts)
counts = apply(taxa_gen_s, 2, sum)
mean(counts)
sd(counts)
sd(counts)/sqrt(length(counts))


###
### convert to relative abundances (first, remove chloroplast/mitochondria and singletons)
###

## call chloroplast (order level) + mitochondria (family level) counts
chloroplast_counts = apply(taxa_order_s[grep("Chloroplast", rownames(taxa_order_s)),], 2, sum)
mitochondria_counts = apply(taxa_fam_s[grep("Mitochondria", rownames(taxa_fam_s)),], 2, sum)

## phylum
taxa_clean = taxa_phy_s
taxa_clean[grep("Cyanob", rownames(taxa_clean)),] = taxa_clean[grep("Cyanob", rownames(taxa_clean)),] - chloroplast_counts
taxa_clean[grep("Proteo", rownames(taxa_clean)),] = taxa_clean[grep("Proteo", rownames(taxa_clean)),] - mitochondria_counts
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]

taxa_clean = as.matrix(taxa_clean)
taxa_clean[which(taxa_clean<0)] = 0
taxa_clean = as.data.frame(taxa_clean)

taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_phy_RA = taxa_clean_RA

## class
taxa_clean = taxa_class_s
taxa_clean[grep("Oxyphotobacteria", rownames(taxa_clean)),] = taxa_clean[grep("Oxyphotobacteria", rownames(taxa_clean)),] - chloroplast_counts
taxa_clean[grep("Alphaproteo", rownames(taxa_clean)),] = taxa_clean[grep("Alphaproteo", rownames(taxa_clean)),] - mitochondria_counts
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]

taxa_clean = as.matrix(taxa_clean)
taxa_clean[which(taxa_clean<0)] = 0
taxa_clean = as.data.frame(taxa_clean)

taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_class_RA = taxa_clean_RA

## order
taxa_clean = taxa_order_s
taxa_clean = taxa_clean[-c(grep("Chloroplast", rownames(taxa_clean))),]
taxa_clean[grep("Rickettsiales", rownames(taxa_clean)),] = taxa_clean[grep("Rickettsiales", rownames(taxa_clean)),] - mitochondria_counts
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]
taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_order_RA = taxa_clean_RA

## family
taxa_clean = taxa_fam_s
taxa_clean = taxa_clean[-c(grep("Chloroplast", rownames(taxa_clean))),] 
taxa_clean = taxa_clean[-c(grep("Mitochondria", rownames(taxa_clean))),] 
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]
taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_fam_RA = taxa_clean_RA

## genus
taxa_clean = taxa_gen_s
taxa_clean = taxa_clean[-c(grep("Chloroplast", rownames(taxa_clean))),] 
taxa_clean = taxa_clean[-c(grep("Mitochondria", rownames(taxa_clean))),] 
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]
taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_gen_RA = taxa_clean_RA

## export relative abundance tables
#write.table(taxa_phy_RA, "sample_data/taxa_tables/rel_abundance-phylum.txt", sep = "\t", quote = FALSE)
#write.table(taxa_class_RA, "sample_data/taxa_tables/rel_abundance-class.txt", sep = "\t", quote = FALSE)
#write.table(taxa_order_RA, "sample_data/taxa_tables/rel_abundance-order.txt", sep = "\t", quote = FALSE)
#write.table(taxa_fam_RA, "sample_data/taxa_tables/rel_abundance-family.txt", sep = "\t", quote = FALSE)
#write.table(taxa_gen_RA, "sample_data/taxa_tables/rel_abundance-genus.txt", sep = "\t", quote = FALSE)

## summary stats
counts = apply(taxa_clean, 2, sum)
mean(counts)
sd(counts)
sd(counts)/sqrt(length(counts))

## check positive control to validate presence of targets
check_CT = taxa_gen_RA[order(taxa_gen_RA$zymo_pos_ct, decreasing = T),]
check_CT$avg = apply(check_CT[,1:50],1,mean)
check_CT[1:10, 51:52]
sum(check_CT[1:8, 51])
# Notes: 
# -positive control 'ZymoBIOMICS Microbial Community DNA Standard'
# -control in line with expected reads/abundances, including Listeria and Salmonella (column 51)
# -pathogens Listeria and Salmonella have 0 reads (0% abundance) in apple microbiome samples (column 52)
# -foodborne pathogens were NOT enriched


###
### set inputs for downstream analysis
###

## taxa table
colnames(taxa_gen_RA)
taxa_table = taxa_gen_RA[,1:50] # remove control
#colnames(taxa_phy_RA)
#taxa_table = taxa_phy_RA[,1:50] # remove control

## metadata table
meta = read.csv("sample_data/metadata/metadata_all.csv", h=T, row.names=1)
meta = meta[order(meta$sample.id),]
meta$sample.id = str_replace_all(meta$sample.id, "-", ".")
meta = meta[grep("ZG", meta$sample.id),]
meta = meta[grep("F", meta$is_control),]
head(meta)

# check that taxa and metadata tables are in the same sample order
meta$sample.id == colnames(taxa_table)


#######################
### ALPHA DIVERSITY ###
#######################

## Shannon and Simpson diversity metrics
diversity_vec = matrix(nrow = dim(taxa_table)[2], ncol = 2)
diversity_vec = as.data.frame(diversity_vec)
for (a in 1:dim(taxa_table)[2]) {
  diversity_vec[a,1] = diversity(taxa_table[,a], index = "shannon")
  diversity_vec[a,2] = diversity(taxa_table[,a], index = "simpson")
}
colnames(diversity_vec) = c("Shannon", "Simpson")
diversity_vec

## add metric to metadata table
meta$Shannon = diversity_vec$Shannon
meta$Simpson = diversity_vec$Simpson
head(meta)
meta$alpha_shape = c(21)
meta$alpha_shape[grep("5", meta$Enrich_pH)] = 24
meta$alpha_shape[grep("7", meta$Enrich_pH)] = 22

## plot alpha for all samples
ggplot(meta, 
       aes(x = factor(Time_Point,
                      levels = c("Native", "Enrich", "Digest")), 
           y = Shannon,
           fill = Time_Point)) +
  geom_boxplot(outlier.size = 0, alpha=0.7, outlier.color = "white") +
  geom_jitter(size = 1, alpha=0.5, pch=meta$alpha_shape, width = 0.15) +
  theme_bw() +
  ylab("Shannon Index") +
  ylim(0, 4.2) +
  scale_fill_manual(values = c("yellow", "red","blue")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14))

## plot alpha, distinguishing apple type and pH
ggplot(meta[grep("Nat|En|Dig", meta$Time_Point),], 
  aes(x = factor(Time_Point,
                 levels = c("Native", "Enrich", "Digest")), 
      fill = factor(Enrich_pH,
                    levels = c("Native", "5", "7")),
      col = factor(Enrich_pH,
                   levels = c("Native", "5", "7")),
      y = Shannon)) +
  geom_boxplot(position = position_dodge(width=0.8), alpha=0.8,
               outlier.size = 0, outlier.color = "white") +
  geom_jitter(position = position_dodge(width=0.8), 
              size = 1.5, alpha=0.5, pch=21, col="black") +
  theme_bw() +
  ylab("Shannon Index") +
  facet_wrap(~Apple_type) +
  theme(legend.text = element_text(size = 11),
        legend.title = element_blank(),
        #legend.position = "none",
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle=45, hjust=1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

## stats

# native vs. enrich vs. digest
summary(aov(Shannon~Time_Point, meta)) # Sample type (native, enrich, digest): SIGNIFICANT
TukeyHSD(aov(Shannon~Time_Point, meta))

# enrich: effect of apple type and pH
summary(aov(Shannon~Enrich_pH*Apple_type, meta[grep("Enrich", meta$Time_Point),])) # Apple type & pH: not sig
TukeyHSD(aov(Shannon~Enrich_pH*Apple_type, meta[grep("Enrich", meta$Time_Point),]))

# digest: effect of apple type and pH
summary(aov(Shannon~Enrich_pH*Apple_type, meta[grep("Digest", meta$Time_Point),])) # Apple type & pH: not sig, though apple type is p=0.095
TukeyHSD(aov(Shannon~Enrich_pH*Apple_type, meta[grep("Digest", meta$Time_Point),]))

# digest: effect of digest matrix
summary(aov(Shannon~Digest_matrix, meta[grep("Digest", meta$Time_Point),])) # Digest matrix: not sig


######################
### BETA DIVERSITY ###
######################

###
### PCoA: All samples
###

# beta-diversity measure
beta <- vegdist(t(taxa_table), 'bray', binary = T)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# add metadata
ord$Apple_type = meta$Apple_type
ord$Sample_type = factor(meta$Time_Point,
                  levels = c("Native", "Enrich", "Digest"))
ord$pH = factor(meta$Enrich_pH,
                levels = c("Native", "5", "7"))
ord$round = meta$Enrich_round
ord$digest = meta$Digest_matrix

# plot
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, 
                       fill = Sample_type,
                       shape = pH)) +
  geom_point(alpha = 0.8, size=4, color = "black") +
  theme_bw() +
  xlab('PCoA1 (15.3%)') +
  ylab('PCoA2 (7.6%)') +
  scale_shape_manual(values = c(21, 24, 22),
                     labels = c("Native", "Enrich-pH-5", "Enrich-pH-7")) +
  scale_fill_manual(values = c("blue", "red", "yellow"),
                    labels = c("Digest_Native", "Digest_Enrich-pH-5", "Digest_Enrich-pH-7")) +
  guides(shape = guide_legend(override.aes = list(fill = c("blue", "red", "red"),
                                                  shape = c(21, 24, 22)))) +
  guides(fill = guide_legend(override.aes = list(fill = c("yellow"),
                                                  shape = c(21, 24, 22)))) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank()) 


# stats
adonis(beta ~ Sample_type, data = ord, permutations = 999)$aov.tab # Sample type (native, enrich, digest): SIGNIFICANT
adonis(beta ~ Apple_type, data = ord, permutations = 999)$aov.tab # Apple type: not sig
adonis(beta ~ pH, data = ord, permutations = 999)$aov.tab # pH: not sig


###
### PCoA: Enrichment conditions (subset enrich samples)
###

# beta-diversity measure
beta <- vegdist(t(taxa_table[,grep("Enrich", meta$Time_Point)]), 'bray', binary = T)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# add metadata
ord$Apple_type = meta$Apple_type[grep("Enrich", meta$Time_Point)]
ord$Sample_type = factor(meta$Time_Point[grep("Enrich", meta$Time_Point)], 
                         levels = c("Native", "Enrich", "Digest"))
ord$pH = factor(meta$Enrich_pH[grep("Enrich", meta$Time_Point)], 
                levels = c("Native", "5", "7"))
ord$round = meta$Enrich_round[grep("Enrich", meta$Time_Point)]

# plot
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, 
                       shape = pH,
                       #fill = Apple, color = Apple
                       fill = round, color = round
                       )) +
  geom_point(alpha = 0.8, size=4, stroke = 1.2) +
  theme_bw() +
  xlab('PCoA1 (27.4%)') +
  ylab('PCoA2 (17.0%)') +
  scale_shape_manual(values = c(21, 22)) +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_color_brewer(palette = "Dark2") +
  theme(axis.text = element_text(size=14, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black")) 

# stats
adonis(beta ~ Apple_type, data = ord, permutations = 999)$aov.tab # Apple type: not sig
adonis(beta ~ pH, data = ord, permutations = 999)$aov.tab # pH: not sig
adonis(beta ~ round, data = ord, permutations = 999)$aov.tab # Round: SIGNIFICANT


## Within each 'round', evaluate effect of pH

# round 1: not sig
beta <- vegdist(t(taxa_table[,intersect(grep("Enrich", meta$Time_Point),
                                        grep("1", meta$Enrich_round))]), 'bray', binary = T)
adonis(beta ~ Enrich_pH,  permutations = 999,
       data = meta[intersect(grep("Enrich", meta$Time_Point),
                             grep("1", meta$Enrich_round)),])$aov.tab

# round 2: not sig
beta <- vegdist(t(taxa_table[,intersect(grep("Enrich", meta$Time_Point),
                                        grep("2", meta$Enrich_round))]), 'bray', binary = T)
adonis(beta ~ Enrich_pH,  permutations = 999,
       data = meta[intersect(grep("Enrich", meta$Time_Point),
                             grep("2", meta$Enrich_round)),])$aov.tab

# round 3: p=0.06
beta <- vegdist(t(taxa_table[,intersect(grep("Enrich", meta$Time_Point),
                                        grep("3", meta$Enrich_round))]), 'bray', binary = T)
adonis(beta ~ Enrich_pH,  permutations = 999,
       data = meta[intersect(grep("Enrich", meta$Time_Point),
                             grep("3", meta$Enrich_round)),])$aov.tab

###
### PCoA: Digest conditions (subset digest samples)
###

# beta-diversity measure
beta <- vegdist(t(taxa_table[,grep("Digest", meta$Time_Point)]), 'bray', binary = T)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# add metadata
ord$Apple_type = meta$Apple_type[grep("Digest", meta$Time_Point)]
ord$pH = factor(meta$Enrich_pH[grep("Digest", meta$Time_Point)], levels = c("Native", "5", "7"))
ord$round = meta$Enrich_round[grep("Digest", meta$Time_Point)]
ord$digest = meta$Digest_matrix[grep("Digest", meta$Time_Point)]

# plot
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, 
                       #shape = Apple,
                       #shape = pH,
                       #fill = pH, color = pH
                       #shape = round,
                       fill = digest, color = digest
                       #fill = Apple, color = Apple
)) +
  geom_hline(yintercept = 0, lty = 2, alpha=0.7) +
  geom_vline(xintercept = 0, lty = 2, alpha=0.7) +
  geom_point(alpha = 0.7,
             size=4, col = "black", pch=21) +
  stat_ellipse(level=0.4, alpha=0.5) +
  theme_bw() +
  xlab('PCoA1 (9.5%)') +
  ylab('PCoA2 (8.3%)') +
  #scale_shape_manual(values = c(21, 22, 1)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text = element_text(size=14, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black")) 

# stats
adonis(beta ~ Apple_type, data = ord, permutations = 999)$aov.tab # Apple type: not sig
adonis(beta ~ pH, data = ord, permutations = 999)$aov.tab # pH: not sig
adonis(beta ~ round, data = ord, permutations = 999)$aov.tab # Round: approach significance (p=0.067)
adonis(beta ~ digest, data = ord, permutations = 999)$aov.tab # Digest matrix: SIGNIFICANT (p=0.002)


#################################
### RELATIVE ABUNDANCE: GENUS ###
#################################

###
### Genus
###

## set table
taxa_plot = taxa_table

## add average and re-order table
taxa_plot$avg = apply(taxa_plot, 1, mean, na.rm=TRUE)
taxa_plot$avg_N = apply(taxa_plot[,grep("N", meta$Time_Point)], 1, mean, na.rm=TRUE)
taxa_plot$avg_E = apply(taxa_plot[,grep("E", meta$Time_Point)], 1, mean, na.rm=TRUE)
taxa_plot$avg_D = apply(taxa_plot[,grep("D", meta$Time_Point)], 1, mean, na.rm=TRUE)
taxa_plot = taxa_plot[order(taxa_plot$avg_E, decreasing = T),]
taxa_plot[1:10, c(51:54)]
sum(taxa_plot$avg_N[1:8])
sum(taxa_plot$avg_E[1:8])
sum(taxa_plot$avg_D[1:8])

# total genera detected
length(which(taxa_plot$avg_N > 0))
length(which(taxa_plot$avg_E > 0))
length(which(taxa_plot$avg_D > 0))

## check for potential pathogen presence (NONE FROM LIST BELOW DETECTED)
taxa_plot[grep("Salmonella|Escherichia|Listeria|Enterococcus|Staphy", rownames(taxa_plot)), 51:54]

## subset key taxa
taxa_plot_select = rbind(apply(taxa_plot, 2, sum),
                         apply(taxa_plot[grep("_Acidobacteria", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Actinobacteria", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Bacteroidetes", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Chloroflexi", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Firmicutes", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Planctomycetes", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Proteobacteria", rownames(taxa_plot)),], 2, sum),
                         taxa_plot[c(1:8),])
rownames(taxa_plot_select)
rownames(taxa_plot_select) = c("D_0__Bacteria", 
                               "D_0__Bacteria;D_1__Acidobacteria",
                               "D_0__Bacteria;D_1__Actinobacteria",
                               "D_0__Bacteria;D_1__Bacteroidetes",
                               "D_0__Bacteria;D_1__Chloroflexi",
                               "D_0__Bacteria;D_1__Firmicutes",
                               "D_0__Bacteria;D_1__Planctomycetes",
                               "D_0__Bacteria;D_1__Proteobacteria",
                               rownames(taxa_plot_select[c(9:16),]))
head(taxa_plot_select)

# clean taxa
rownames(taxa_plot_select)
taxa_plot_select[grep("^D_0__Bacteria;D_1__Actinobacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria;D_1__Actinobacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Microbacteriaceae", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^D_0__Bacteria;D_1__Firmicutes$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria;D_1__Firmicutes$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Bacillus|Lactococcus", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^D_0__Bacteria;D_1__Proteobacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria;D_1__Proteobacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Gammaproteobacteria|Alphaproteobacteria", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^D_0__Bacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep(";", rownames(taxa_plot_select)),], 2, sum)

# check sums
apply(taxa_plot_select,2,sum)
min(taxa_plot_select)

## prep plot
rownames(taxa_plot_select)
rownames(taxa_plot_select) = c("Other", 
                               "p__Acidobacteria", "p__Actinobacteria", "p__Bacteroidetes", 
                               "p__Chloroflexi", "p__Firmicutes", "p__Planctomycetes", "p__Proteobacteria", 
                               "Pantoea", "Gluconobacter", "Bacillus", 
                               "Enterobacteriaceae sp.", "Kosakonia", "Rahnella",
                               "Pseudomonas", "Microbacteriaceae sp.")

# melt data
colnames(taxa_plot_select)[1:50] == meta$sample.id
taxa_plot_select2 = as.data.frame(t(taxa_plot_select[,1:50]))
taxa_plot_select2$Sample_type = meta$Time_Point
taxa_plot_select2$Apple_type = meta$Apple_type
taxa_plot_select2$pH = meta$Enrich_pH
taxa_plot_select2$round = meta$Enrich_round
taxa_plot_select2$label = meta$Composite_label2
rownames(taxa_plot_select2) = substr(rownames(taxa_plot_select2), 
                                     start = 4, stop = nchar(rownames(taxa_plot_select2)))
head(taxa_plot_select2)

taxa_melt = melt(data.frame(taxa_plot_select2,
                            rownames(taxa_plot_select2)))
head(taxa_melt)
colnames(taxa_melt)[c(6,7)] = c("Sample", "Taxon")
head(taxa_melt)

# factor: taxa
taxa_melt$Taxon = factor(taxa_melt$Taxon,
                         levels = c("Other", 
                                    "p__Acidobacteria", "p__Actinobacteria", "Microbacteriaceae sp.",
                                    "p__Bacteroidetes", "p__Chloroflexi", 
                                    "p__Firmicutes", "Bacillus", "p__Planctomycetes", 
                                    "p__Proteobacteria", "Pantoea", "Gluconobacter", 
                                    "Enterobacteriaceae.sp.", "Kosakonia", "Rahnella", "Pseudomonas"))

# factor: type
taxa_melt$Sample_type = factor(taxa_melt$Sample_type,
                               levels = c("Native", "Enrich", "Digest"))

# factor: pH
taxa_melt$pH = as.factor(taxa_melt$pH)
taxa_melt$pH = factor(taxa_melt$pH,
                      levels = c("Native", "7", "5"))


# plot EP
EP = ggplot(taxa_melt[grep("EP", taxa_melt$Apple_type),],
            aes(x=label, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  scale_fill_manual(values=c(
    "lightgray",
    "gold", 
    "brown4", "brown1",  
    "#666633", "#CCFF33", 
    "#0099FF", "#33FFFF", "#CCFFFF", # blues
    "#330066", "#9900CC", "#CC00CC", "#FF33FF", "#FF66FF", "#FF99FF", "#FFCCFF" # earth-tones
  )) +
  facet_grid(~Sample_type, scales="free", space="free") +
  theme(legend.position = "none",
    strip.text = element_blank(),
    legend.text = element_text(size=12),
    axis.text = element_text(size=12),
    axis.text.x = element_text(size=9, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=12, color = "black"))

# plot GD
GD = ggplot(taxa_melt[grep("GD", taxa_melt$Apple_type),],
            aes(x=label, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  scale_fill_manual(values=c(
    "lightgray",
    "gold", 
    "brown4", "brown1",  
    "#666633", "#CCFF33", 
    "#0099FF", "#33FFFF", "#CCFFFF", # blues
    "#330066", "#9900CC", "#CC00CC", "#FF33FF", "#FF66FF", "#FF99FF", "#FFCCFF" # earth-tones
  )) +
  facet_grid(~Sample_type, scales="free", space="free") +
  theme(legend.position = "none",
    strip.text = element_blank(),
    legend.text = element_text(size=12),
    axis.text = element_text(size=12),
    axis.text.x = element_text(size=9, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=12, color = "black"))

# plot all
grid.arrange(EP, GD)


###############################
### qPCR absolute abundance ###
###############################


## plot
ggplot(meta, 
       aes(x = factor(Time_Point,
                      levels = c("Native", "Enrich", "Digest")), 
           y = log_copies_16Sgene_ml,
           fill = Time_Point)) +
  geom_boxplot(outlier.size = 0, alpha=0.7, outlier.color = "white") +
  geom_jitter(size = 1, alpha=0.5, pch=meta$alpha_shape, width = 0.15) +
  theme_bw() +
  ylab("Log copies 16S rRNA gene") +
  ylim(0, 11.2) +
  scale_fill_manual(values = c("yellow", "red","blue")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14))

## stats

# native vs. enrich vs. digest
summary(aov(log_copies_16Sgene_ml~Time_Point, meta)) # Sample type (native, enrich, digest): SIGNIFICANT
TukeyHSD(aov(log_copies_16Sgene_ml~Time_Point, meta))

# enrich: effect of apple type and pH
summary(aov(log_copies_16Sgene_ml~Enrich_pH*Apple_type, meta[grep("Enrich", meta$Time_Point),])) # pH: SIGNIFICANT
TukeyHSD(aov(log_copies_16Sgene_ml~Enrich_pH*Apple_type, meta[grep("Enrich", meta$Time_Point),])) # apple type: not sig

# digest: effect of apple type and pH
summary(aov(log_copies_16Sgene_ml~Enrich_pH*Apple_type, meta[grep("Digest", meta$Time_Point),])) # Apple type & pH: not sig, though interactions are p=0.095
TukeyHSD(aov(log_copies_16Sgene_ml~Enrich_pH*Apple_type, meta[grep("Digest", meta$Time_Point),]))

# digest: effect of digest matrix
summary(aov(log_copies_16Sgene_ml~Digest_matrix, meta[grep("Digest", meta$Time_Point),])) # Digest matrix: SIGNIFICANT

