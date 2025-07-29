# Libraries
library(vcfR)
library(tidyverse)
library(adegenet)
library(poppr)
library(ggsci)
library(RColorBrewer)

my_path <- "/path/to/data/"

############ Loading data and creating R objects #################
#load Rdata
load("pop_gen_objects.Rdata")

#Import VCF
vcf <- read.vcfR(paste(my_path, "pop.vcf.gz", sep=""), verbose = FALSE)

# convert from vcfR to genind/genlight objects
Rsol_gind <- vcfR2genind(vcf)
Rsol_glight <- vcfR2genlight(vcf)

#Save genlight and genind file in R object
save(Rsol_gind, Rsol_glight, file = "pop_gen_objects.Rdata")

#importing metadata
Rsol_metadata <- read_csv(file = "pop_metadata_str.csv", col_names = T)
str(Rsol_metadata)
head(Rsol_metadata)

Rsol_metadata$ID <- as.factor(Rsol_metadata$ID)
Rsol_metadata$Host <- as.factor(Rsol_metadata$Host)
Rsol_metadata$Location <- as.factor(Rsol_metadata$Location)
Rsol_metadata$Year <- as.factor(Rsol_metadata$Year)

############ Preparation for Poppr #################

# Defining ploidy
#Rhizoctonia is a dikaryon so we are doing ploidy = 2
ploids <- 2

#Incorporating metadata
#genind object
strata(Rsol_gind) <- data.frame(Rsol_metadata)
#strata(Rsol_gind) <- subset(Rsol_metadata[,], Rsol_metadata$ID %in% rownames(Rsol_gind@tab))
pop(Rsol_gind) <- strata(Rsol_gind)$Location

#genlight object
pop(Rsol_glight) <- strata(Rsol_gind)$Location
ploidy(Rsol_glight) <- ploids

#Add strata to genlight object
strata(Rsol_glight) <- data.frame(Rsol_metadata)
Rsol_glight

##Analyze population by strata --- Host within location
setPop(Rsol_glight) <- ~Location/Host
setPop(Rsol_gind) <- ~Location
Rsol_gind@pop

############### PCA analysis ###############

# PCA
set.seed(999)
pca1 <- glPca(Rsol_glight, nf=10, parallel = T, useC = F, n.cores = 20)  # where nf is the number of principal components to keep, adjust as needed

save(pca1, file = "newPCA.RData")
load("newPCA.RData")

## Percent variance explained by eigenvalues
eigenvalues <- pca1$eig # eigenvalues
var_explained <- eigenvalues/sum(eigenvalues)  # variance explained by each PC
# keep only two decimal places
var_explained_round <- round(var_explained, digits = 3)

PCs <- pca1$scores # PC scores
PCs <- as.data.frame(PCs) # as dataframe

#Adding metadata
PCs$pop <- pop(Rsol_glight)
PCs$sample <- rownames(PCs)
PCs$host <- strata(Rsol_gind)$Host
PCs$location <- strata(Rsol_gind)$Location

### Create plot ###
p <- ggplot(PCs, aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_point(aes(fill=host, shape=pop), size = 4.5, alpha = 0.45, colour="black", stroke = 0.69) +
  scale_shape_manual(values = 21:24, name = "Location") +
  guides(fill = guide_legend(override.aes = list(shape = 21), title = "Host")) +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right") +
  labs(x = paste("PC1 (", var_explained_round[1]*100, "%)", sep = ""),
       y = paste("PC2 (", var_explained_round[2]*100, "%)", sep = ""))

p

### add elipse to PCA
p <- p + stat_ellipse(level = 0.95, position = "identity", inherit.aes = TRUE)
p

ggsave(plot=p, "PCA_ellipses.png", dpi = 600, units = "in", height = 5, width = 9)


############ Multilocus genotypes ##############
setPop(Rsol_gind) <- ~Location

# Gene clone
Rsol_gclone <- as.genclone(Rsol_gind)
Rsol_gclone
Rsol_gclone@mlg

# SNP clone
Rsol_snpclone <- as.snpclone(Rsol_glight, parallel = T, n.cores = 20)

#Save genclone and snpclone files in R object
save(Rsol_gclone, Rsol_snpclone, file = "pop_clone_objects.Rdata")

## load clone objects
load("pop_clone_objects.Rdata")

## Determine the multilocus genotypes in the pop #####
Rsol_mlg <- mlg(Rsol_gind) ### equal to the number of individuals. Too many unique combinations

# Choosing a threshold ----------------------------------------------------

## Collect information of the thresholds. We can set threshold = 1
# because we know that this will capture the maximum possible distance:
(thresholds <- mlg.filter(Rsol_snpclone, distance = "bitwise.dist", stats = "THRESHOLDS",
                          threshold = 1))

# We can use these thresholds to find an appropriate cutoff
(pcut <- cutoff_predictor(thresholds))

#Default MLG: all alleles must match to make a unique multilocus genotype
#Threshold value equals pcut
(mlg_filter <- mlg.filter(Rsol_snpclone, threshold = 0.03, distance = "bitwise.dist", stats = "ALL"))
mlg_sclone <- mlg.table(Rsol_snpclone, strata = ~Location/Host) #Same number of MLGs as individuals
mlg_vector <- mlg.vector(Rsol_snpclone, reset = FALSE) # Each individual belongs to a single MLG

# See which individuals belong to each MLG
mlg_ind <- mlg.id(Rsol_snpclone) # all are unique

#because the MLGs are unique we will identify MLLs
#Apply a threshold to define MLLs 
#Filtered (“contracted”)
Rsol_snpclone #A threshold of 0.03 also gives 99 contracted multilocus genotypes
mlg_ind2 <- mlg.id(Rsol_snpclone) # See which ind belong to each MLG

mll(Rsol_snpclone) <- "contracted" 


##### ===============================================
# This part of the script is not really necessary
#The threshold estimated above works.
#The following steps serve as a prove of concept

#We can utilize genetic distance, which will allow us to collapse multilocus genotypes that are under a specific distance threshold. 

#Calculate raw genetic distance with bitwise.dist(), this distance is in %, see ?bitwise.dist()
#Fraction of different alleles in percentage (this is actually proportion 0 to 1):
gc.dist <- bitwise.dist(Rsol_snpclone)
gc.dist
hist(gc.dist, breaks = 100000)
median(gc.dist)

#Number of allelic differences. percent = FALSE will return the distance represented as integers from 1 to n where n is the number of loci
gc.dist.nl <- bitwise.dist(Rsol_snpclone, percent = F)
gc.dist.nl
hist(gc.dist.nl, breaks = 10000)
median(gc.dist.nl)
max(gc.dist.nl)

#Threshold based on distance calculated with bitwise.dist (dissimilarity distance)

#The most familiar name might be the Hamming distance, or the number of differences between two strings.
#Should the distance be represented from 0 to 1? Default set to TRUE. 
#FALSE will return the distance represented as integers from 1 to n where n is the number of loci. 
#This option has no effect if euclidean = TRUE
#If the user supplies a genind or genclone object, prevosti.dist() will be used for calculation.
#So here it is Prevosti's distance
#Set again to the original (default) 
mll(Rsol_snpclone) <- "original"
mll(Rsol_snpclone) # original
Rsol_snpclone

#Choosing an algorithm and a threshold to represent the minimum genetic distance at which two individuals would be considered from different clonal lineages.

Rsol_snpclone.filtered <- filter_stats(Rsol_snpclone, distance = bitwise.dist, plot = TRUE)
Rsol_snpclone.filtered

# One method described in the literature of choosing a threshold is to look 
#for an initial, small peak in the histogram of pairwise genetic distances and set the threshold 
#to be between that peak and the larger peak `(Arnaud-Haond et al. 2007, @bailleul2016rclone).
#This initial peak likely represents clones differentiated by a small set of random mutations. 
# 
#Closer look, to identify the initial peak that likely represents clones differentiated by a small set of random mutations
hist(gc.dist, breaks = 10000, xlim= c(0, 0.1)) 
hist(gc.dist, breaks = 10000, xlim= c(0, 0.01)) #It looks like the very first peak is below 0.01

#hist(gc.dist, breaks = 1000000, xlim= c(0, 0.004))
#hist(gc.dist, breaks = 1000000, xlim= c(0, 0.001))

#-- Run from here again --#
#Threshold based on 0.01
mlg.filter(Rsol_snpclone, distance = gc.dist, algorithm = "a") <- 0.01
Rsol_snpclone #A threshold of 0.01 gives 99 contracted multilocus genotypes
#$SIZES
#[1]  1  1  0  2  1  1  2  1  1  1  1  1  1  1  1  0  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  2  0  0  3  1  1  1  1  1  1  1  1  1
#[45]  1  2  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  8  0  0  0  0  0  5  0  0  0  1  2  1  1  1  1
#[89]  0  0  1  4  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0  2  0  5  0  0  0  1  0  0  0  2  0  1  1  1  2  0  1  1  1
#[133]  1  0  1 18  0  1  1  1  0  2  1  1  1

mlg_ind <- mlg.id(Rsol_snpclone)

#Repeat with previously estimated threshold
Rsol_clones <- mlg.filter(Rsol_snpclone, threshold = 0.03, distance = "bitwise.dist", stats = "ALL")
Rsol_snpclone #A threshold of 0.03 also gives 99 contracted multilocus genotypes
mlg_ind2 <- mlg.id(Rsol_snpclone) 

mlg.table(Rsol_snpclone)

####
mlg_table <- t(mlg.table(Rsol_snpclone))
mlg_table <- data.frame(mlg_table) %>% rownames_to_column(var = "MLG") %>% 
  mutate_if(is.integer, as.double)
mlg_table$Sums <- rowSums(mlg_table[,2:11])
mlg_table[mlg_table == 0] <- NA

mlg_table_sort <- mlg_table[order(-mlg_table$Sums),]
order <- mlg_table_sort$MLG
#mlg_table <- filter(mlg_table, Sums > 1) %>% select(!contains("Sums")) %>% pivot_longer(!MLG, names_to = "subpop", values_to = "count")
mlg_table_final <- arrange_all(mlg_table) %>% pivot_longer(!MLG, names_to = "subpop", values_to = "count") %>%
  filter(subpop != "Sums") 


(ranges <- ggplot(mlg_table_final, aes(x = subpop, y = MLG, group = MLG)) + 
    #geom_line(aes(color = MLG), size = 1, linetype = 1) + 
    geom_point(aes(color = MLG, size = count*3), pch = 21, fill = "white", stroke = 2) +
    geom_text(aes(label = count), size = 2.5, fontface = "bold") + 
    #scale_color_manual(values = myPal) + 
    #guides(color = guide_legend(ncol = 3)) +
    ylab("Multilocus Genotype") +
    xlab("Subpopulation") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 7, color = "black"),
           text = element_text(family = "Helvetica"),
           legend.position = "none",
           axis.line = element_line(colour = "black")) +
    #    theme(panel.grid.major.y = element_line(size = 0.5, color = "grey")) +
    #    theme(panel.grid.major.x = element_line(size = 1, color = "grey")) 
    scale_y_discrete(limits = rev(order)))

ggsave("mlg_range.pdf", device = "pdf", plot = ranges,
       width = unit(5, "in"), height = unit(12, "in"))


(mlg_counts <- ggplot(mlg_table, aes(x = MLG, y = Sums, fill = MLG)) +
  geom_bar(stat = "identity", colour = "grey50", linewidth = 0.3) +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,20), position = "bottom", guide = guide_axis(position = "top")) +
  guides(fill = guide_legend(ncol = 2, title = NULL)) +
  geom_text(aes(label = Sums), size = 2.5, hjust = 0, fontface = "bold") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        #element_text(angle = -45, vjust = 1, hjust = 0, size = 5),
        #legend.position = c(1, 0), legend.direction = "vertical",
        legend.position = "none",
        #legend.text = element_text(size = rel(0.65)),
        #legend.justification = c(1,0),
        #legend.background = element_rect(color = "black"),
        #legend.margin = grid::unit(0, "pt"),
        #legend.key.size = grid::unit(15, "pt"),
        text = element_text(family = "Helvetica"),
        axis.title.y = element_blank()) +
  scale_x_discrete(limits = rev(order)) +
  ylab("Count") +
  coord_flip())

ggsave("mlg_dist.svg", plot = set_panel_size(mlg_counts, width = unit(3, "in"),
                                             height = unit(12, "in")),
                                             width = 6, height = 15)

####

Rsol_clones$MLGS
Rsol_clones$SIZES

mll(Rsol_snpclone) <- "contracted"

(mlls <- mll(Rsol_snpclone, type = NULL))
(mlls <- mll(Rsol_gclone, type = NULL)) # This generates the multilocus lineages

#(mll_gc <- mll(Rsol_snpclone, type = NULL, mlg_filter$MLGS))
##mlg_gclone <- mlg.table(mll_gclone, strata = ~Location/Host)

#Find all multilocus genotypes that cross populations 
mlg_gc.cp <- mlg.crosspop(Rsol_snpclone, strata = ~Location/Host, df = T)
mlg_gc.cp2 <- mlg.crosspop(Rsol_snpclone, strata = ~Location, df = T)
mlg_gc.cp3 <- mlg.crosspop(Rsol_snpclone, strata = ~Host, df = T)

ggsave("Rsol_mlg_Loc-Host_new.pdf")

mlg_snpclone_loc <- mlg.table(Rsol_snpclone, strata = ~Location)
mlg_snpclone_host <- mlg.table(Rsol_snpclone, strata = ~Host)

#run AMOVA
amova_of_Rsol_snpclone <- poppr.amova(Rsol_snpclone,
                                      #hier = ~Location,
                                      hier = ~Location/Host,
                                      clonecorrect = T,
                                      within = T,
                                      filter = FALSE,
                                      threshold = 0.03,
                                      algorithm = "farthest_neighbor",
                                      threads = 12,
                                      missing = "loci",
                                      method = "ade4",
                                      #nperm = 0
)

randtest_snpclone <- randtest(amova_of_Rsol_snpclone, 999)
plot(randtest_snpclone)

############ ==========================================================
rowSums(mlg_gclone != 0)
#AR_Corn AR_Grass  AR_Rice  AR_Sorg   AR_Soy  CU_Bean  LA_Rice   LA_Soy  TX_Rice   TX_Soy 
#1        1       33        1       25        1       29       38       15        1 

#AR_Bermuda grass          AR_Corn          AR_Rice       AR_Sorghum       AR_Soybean 
#1                1               26                1                6 
#Cuba_Common bean          LA_Rice       LA_Soybean          TX_Rice 
#1                3               34                1 


##### Calculate diversity statistics for MLG s#######
diversity_stats(mlg_gc)
diversity_stats(mlg_gc_loc)
diversity_stats(mlg_gc_host)


##### Calculate Gst for the population #########
##Calculate index by subpopulation = location

gi_AR <- popsub(Rsol_gi, sublist = "AR")
gi_LA <- popsub(Rsol_gi, sublist = "LA")

setPop(Rsol_gi) <- ~Location/Host
gi_AR_soy <- popsub(Rsol_gi, sublist = "AR_Soybean")
gi_AR_rice <- popsub(Rsol_gi, sublist = "AR_Rice")
gi_LA_soy <- popsub(Rsol_gi, sublist = "LA_Soybean")
gi_LA_rice <- popsub(Rsol_gi, sublist = "LA_Rice")

#Median pairwise genetic distance within location
#AR
AR_dist <- bitwise.dist(gi_AR)
#ARdist
median(AR_dist)
#[1] 0.1597418

#LA
LA_dist <- bitwise.dist(gi_LA)
#LAdist
median(LA_dist)
#[1] 0.2529625

#AR_soybean
AR_soy_dist <- bitwise.dist(gi_AR_soy)
median(AR_rice_dist)
#1] 0.1729424

#AR_rice
AR_rice_dist <- bitwise.dist(gi_AR_rice)
median(AR_rice_dist)
#[1] 0.1280537

#AR_soybean
AR_soy_dist <- bitwise.dist(gi_AR_soy)
median(AR_soy_dist)
#[1] 0.1729424

#AR_rice
AR_rice_dist <- bitwise.dist(gi_AR_rice)
median(AR_rice_dist)
#[1] 0.1280537

#LA_soybean
LA_soy_dist <- bitwise.dist(gi_LA_soy)
median(LA_soy_dist)
#[1] 0.2478683

#LA_rice
LA_rice_dist <- bitwise.dist(gi_LA_rice)
median(LA_rice_dist)
#[1] 0.2458453

#Expected heterozygosity (Hexp) with Hs function (adegenet)
#Hs for each Country
setPop(Rsol_gi) <- ~Location
Hs.location <- Hs(Rsol_gi)
Hs.location
#       AR        LA      Cuba        TX 
#0.1983948 0.3399681 0.1497989 0.2846491

setPop(Rsol_gi) <- ~Location/Host
Hs.location.host <- Hs(Rsol_gi)
Hs.location.host
#AR_Soybean       LA_Soybean Cuba_Common bean          AR_Rice 
#0.19419880       0.32261516       0.14979890       0.18293476 
#AR_Bermuda grass          LA_Rice          AR_Corn          TX_Rice 
#0.08428795       0.29424904       0.09012410       0.28464912 
#AR_Sorghum 
#0.10330970


#Test difference in expected heterozygosity (Gene diversity)
#Pairwise AR-LA
AR_LA_Ht <- Hs.test(Rsol_gi[pop="AR"], Rsol_gi[pop="LA"], n.sim=499)
AR_LA_Ht
plot(AR_LA_Ht)

#Global GST
#Clades
diff.AR <- diff_stats(gi_AR) # this function calculates overall Nei's Gst, Hedrick's Gst and  of the dataset
diff.LA <- diff_stats(gi_LA)

#Pairwise genetic differentiation using mmod package and genind object
library(mmod)
#vignette("mmod-demo", package="mmod")

setPop(Rsol_gi) <- ~Location
popNames(Rsol_gi)
pair.diff.location  <- pairwise_Gst_Nei(Rsol_gi, linearized = FALSE) # Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- Gst)

setPop(Rsol_gi) <- ~Host
popNames(Rsol_gi)
pair.diff.host  <- pairwise_Gst_Nei(Rsol_gi, linearized = FALSE)

setPop(Rsol_gi) <- ~Location/Host
popNames(Rsol_gi)
pair.diff.location.host  <- pairwise_Gst_Nei(Rsol_gi, linearized = FALSE)


######## Genetic differentiation ############
Rsol_gi_cc <- clonecorrect(Rsol_gi, strata = ~Location, combine = T)
Rsol_gi_cc

Rsol_gc_filtered <- filter_stats(Rsol_gc, distance = bitwise.dist, plot=T)
Rsol_gc_filtered

print(farthest_thresh <- cutoff_predictor(Rsol_gc_filtered$farthest$THRESHOLDS))
##[1] 0.04191202
print(average_thresh <- cutoff_predictor(Rsol_gc_filtered$average$THRESHOLDS))
##[1] 0.03224383
print(nearest_thresh <- cutoff_predictor(Rsol_gc_filtered$nearest$THRESHOLDS))
##[1] 0.02409391

##### Diversity statistics with poppr ######

Rsol_stats <- poppr(Rsol_gc, strata = ~Location/Host, missing = "zero", legend = T, parallel = T, n.cores = 20)



############ Minimum spanning network ###################

# calculate and create MSN object
#set population to location
setPop(Rsol_gi) <- ~ Location
#Calculate dissimilarity distances
Rsol_gi_dist <- bitwise.dist(Rsol_gi, percent = TRUE, mat = FALSE, threads = 20L)

#Estimate  minimum spanning network
min_span_net <- poppr.msn(Rsol_gi, Rsol_gi_dist, showplot = FALSE, include.ties = TRUE)

# calculate edge weight cutoff, emphasizing clone
library(igraph)
min_span_net_dists <- E(min_span_net$graph)$weight

# a cutoff of the 90th percentile was selected for edges that should be plotted
cutoff_dist <- quantile(min_span_net_dists, 0.99)

# plot MSN
set.seed(999)
plot_msn <- plot_poppr_msn(Rsol_gi,
                           min_span_net,
                           inds = "NONE",
                           mlg = FALSE,
                           gadj = 50,
                           gweight = 1,
                           palette = spectral,
                           beforecut=TRUE,
                           cutoff = cutoff_dist,
                           quantiles = FALSE,
                           pop.leg = TRUE,
                           size.leg = TRUE,
                           scale.leg = TRUE,
                           nodelab = 30,
                           layfun = igraph::layout_with_gem,
                           nodescale = 50
)






diversity_stats(Rsol_gl)
Rsol.mlg <- mlg.filter()

### DartR for exporting vcf to structure

library(dartR)
gl2structure(Rsol_gl, ploidy = 2, exportMarkerNames = T,
             outfile = "Rsol_str.str")

gl2faststructure(Rsol_gl, outfile = "Rsol_fstr.str", outpath = ".", verbose = 5)


Rsol_Ia <- bitwise.ia(Rsol_gl, missing_match = F, differences_only = F, threads = 10)
#[1] -0.7979109


library(hierfstat)
Rsol_hf <- genind2hierfstat(Rsol_gi)
Rsol_dos <- fstat2dos(Rsol_hf[,-1])
#
Rsol_fs <- hierfstat::fs.dosage(Rsol_dos, pop = Rsol_hf[,1])
#
# AR_Soybean LA_Soybean Cuba_Common bean AR_Rice
# Fis    -0.2731    -0.0064              NaN -0.0574
# Fst     0.3980     0.0011              NaN  0.4344
# AR_Bermuda grass LA_Rice AR_Corn TX_Rice AR_Sorghum
# Fis              NaN -0.4436     NaN -0.9991        NaN
# Fst              NaN  0.0675     NaN  0.1330        NaN
# All
# Fis -0.3559
# Fst  0.2068

test_pairwise_fst <- genet.dist(Rsol_hf, method = "WC84")

#Saving objects
save(Rsol_hf, Rsol_dos, file = "pop_gen_hierfstat.Rdata")

#run AMOVA
amova_of_genid2snp2 <- poppr.amova(Rsol_snpcl,
                                   hier = ~Location/Host,
                                   clonecorrect = F,
                                   within = T,
                                   filter = FALSE,
                                   threshold = 0.04,
                                   algorithm = "farthest_neighbor",
                                   threads = 12,
                                   missing = "loci",
                                   method = "ade4",
                                   #nperm = 0
)
randtest
set.seed(123)
signif_genid2snp <- ade4::randtest(amova_of_genid2snp, nrepet = 999)

set.seed(123)
signif_genid2snp2 <- ade4::randtest(amova_of_genid2snp2, nrepet = 999)
