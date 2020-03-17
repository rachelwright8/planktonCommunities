setwd("~/Dropbox/andrea_plankton/clean_for_upload/2dada_deseq/")

# DOWNLOADING (installing) packages ----
# YOU ONLY HAVE TO execute just once when first using a script:
# source("https://bioconductor.org/biocLite.R")
# biocLite("labdsv")
# biocLite("dada2")
# biocLite('phyloseq')
# biocLite('ShortRead')
# install.packages("tidyverse")
# install.packages("vegan")
# install.packages("MCMC.OTU")
# install.packages("cluster")
# install.packages("pheatmap")

# Load the libraries you need (DO THIS EVERY TIME YOU OPEN A SCRIPT)
library(labdsv)  # this loads MASS, which conflicts with the "select" function in dplyr(tidyverse). Load tidyverse last.
library(dada2) # for dada2
library(ShortRead) # for dada2
library(phyloseq) # for phyloseq
library(tidyverse) # for data wrangling and ggplot2
library(vegan) # for adonis
library(MCMC.OTU) # for MCMC.oTU
library(cluster)
library(pheatmap) # for pretty heatmaps
library(DESeq2)

# Making sample information table with raw data-----
sample_info <- read.delim("orig_raw_data_total.mapping.txt")
head(sample_info)
summary(sample_info$Site)

# Make samples 37, 38, 39 mid instead of late (labeling error- confirmed with lab notebook)
sample_info[36, 6] = "Mid"  
sample_info[37, 6] = "Mid"
sample_info[38, 6] = "Mid"
sample_info[75, 6] = "Mid"  
sample_info[76, 6] = "Mid"
sample_info[77, 6] = "Mid"
# remove samples 1-3 for both forward and reverse reads because of sampling depth inconsistency (confirmed with field notebook)
sample_info <- sample_info[-c(1,2,39,40,41), ]

summary(sample_info$Site)
# BastimentosN BastimentosS    Cristobal     DragoMar   PopaIsland  PuntaDonato  PuntaLaurel    STRIPoint 
# 6            6            6            6            6            6            6           30
summary(sample_info$Time)
# Early  Late   Mid 
# 12    12    48

# Make a vector called `inshore_sites` that lists all of the inshore sites
inshore_sites <- c("PuntaDonato", "STRIPoint", "Cristobal", "PuntaLaurel")

# Make a new column called `siteType` that (as a factor) enters the text "inshore" if the site name is contained in the vector `inshore_sites` and enters the text "offshore" if it isn't
sample_info$siteType <- as.factor(ifelse(sample_info$Site %in% inshore_sites, "inshore","offshore"))
summary(sample_info)

# Rename the column called `Number` to `tech_rep` and make it a factor
names(sample_info)
colnames(sample_info)[8] <- "techRep"

# Get rid of columns you don't need. Only keep SampleID, Site, Time, techRep, and siteType
names(sample_info)
sam_info <- sample_info %>% 
            dplyr::select(SampleID, Site, Time, techRep, siteType)
head(sam_info)

# Set path to unzipped, renamed fastq files (sequencing samples) 
path <- "../Plankton_data/"
fns <- list.files(path)
fns

# remove samples 1-3 for both forward and reverse reads because of sampling depth inconsistency (confirmed with field notebook)
fns <- fns[-c(2,3,4,5,6,7,8,9,10,11)]

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files

# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq; OTHERWISE MODIFY
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) #the last number will select the field for renaming
head(sample.names)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Visualize raw data 

# First, lets look at quality profile of R1 reads. Plot the first and last 4 samples.
# This function plots a visual summary of the distribution of quality scores as a function of sequence position for the input fastq file.
plotQualityProfile(fnFs[c(1:4)])
plotQualityProfile(fnFs[c(71:72)])
# Where do the base call qualities get lower than ~30? -----> ~250 bp in forward reads

# Then look at quality profile of R2 reads
plotQualityProfile(fnRs[c(1:4)])
plotQualityProfile(fnRs[c(71:72)])
# Where do the base call qualities get lower than ~30? -----> ~200 bp in reverse reads
# The reverse reads are significantly worse quality, especially at the end, common in Illumina sequencing.
# This isn’t too worrisome, DADA2 incorporates quality information into its error model which makes the algorithm more robust, 
# but trimming as the average qualities crash is still a good idea as long as our reads will still overlap. 

# Recommend trimming where quality profile crashes ~30

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter and Trim (this takes awhile) 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
              truncLen=c(250,200),
              maxN=0, # DADA does not allow Ns
              maxEE=c(1,1), # allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
              truncQ=2, # truncate reads at the first instance of a quality score less than or equal to 2
              trimLeft=c(24,19), #N nucleotides to remove from the start of each read to remove sequencing primers
              rm.phix=TRUE, # remove reads matching phiX genome
              matchIDs=TRUE, # enforce matching between id-line sequence identifiers of F and R reads
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
tail(out)
class(out)

# How many reads did we lose?
summary(out)

out_stats <- as.data.frame(out) %>% mutate(perc_reads_remaining = reads.out/reads.in*100)
mean(out_stats$perc_reads_remaining) # we only lost 20% of the reads
sum(out_stats)

# Save the out file 
# save(sam_info, out, filtFs, filtRs, sample.names, file="outData.RData")

# For DADA2 alaysis-------
load("outData.RData")

setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
# Forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
# Reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate reads

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer Sequence Variants 
# dada takes a long time

setDadaOpt(BAND_SIZE=32)

# DADA analysis on forward and reverse reads 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[70]]
dadaRs[[70]]
head(dadaFs)

# Merge paired reads. Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample

head(mergers[[70]])
summary((mergers[[70]]))

# We now have a data.frame for each sample with the merged $sequence, its $abundance, and the 
# indices of the merged $forward and $reverse denoised sequences. Paired reads that did not 
# exactly overlap were removed by mergePairs.

# Construct sequence table
# a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab))), xlab="BP", ylab="abundance", main="Histogram of sequence lengths")

# 365-386 based on histogram of sequence lengths and where there were more (at a faily conservative length)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(365 ,386)] 
table(nchar(getSequences(seqtab2)))
dim(seqtab2)

 
# Remove chimeras 

# The core dada method removes substitution and indel errors, but chimeras remain. 
# Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
# than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
# a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)


# The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
# but can be substantial. Here chimeras make up about 36% of the inferred sequence variants (138-89 = 49 => 49/138), 
# BUT those variants account for only about 0.5% of the total sequence reads
# Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)

# Track Read Stats -----

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track) 

write.csv(track,file="final_sequence.csv",row.names=TRUE,quote=FALSE)

# Checking how much was lost over all through out the analysis
track <- read.csv("final_sequence.csv")
head(track)

tracklost <- as.data.frame(track) %>% 
  mutate(remaining_afterfilter = nonchim/input*100)
mean(tracklost$remaining_afterfilter) 
# lost 62% of reads

# How much was lost at each step all looks good 
tracklostind <- as.data.frame(track) %>% 
 mutate(remaining_filtered = filtered/input*100) %>% # lost 20% of reads
 mutate(remaining_denoised = denoised/filtered*100) %>% # lost 0% of reads
 mutate(remaining_merged = merged/denoised*100) %>% # lost 49% of reads
 mutate(remaining_tabled = tabled/merged*100) %>% # lost <1% of reads
 mutate(remaining_nonchim = nonchim/tabled*100) # lost <1% of reads
mean(tracklostind$remaining_merged) 

head(tracklost)
write.csv(tracklost,file="track_lost.csv",row.names=T,quote=F)

# Assign Taxonomy
# It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to classify sequence variants taxonomically. 
# The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, and outputs 
# the taxonomic assignments with at least minBoot bootstrap confidence.

# Silva downloaded from database and then matched 
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa", 
                       minBoot = 50,
                       multithread = TRUE,
                       tryRC = TRUE,
                       outputBootstraps = FALSE)
summary(taxa)

# Rownames of the sample variable table (sam_info) and the OTU table (seqtab.nochim) don't match
rownames(seqtab.nochim)
rownames(sam_info) <- sam_info$SampleID
rownames(sam_info)

# Make them match
rownames(seqtab.nochim) <- sub("-",".",rownames(seqtab.nochim))
rownames(seqtab.nochim)

identical(sort(rownames(seqtab.nochim)),sort(rownames(sam_info)))
# they match now!

# START PHYLOSEQ (load) ------
#save(sam_info, seqtab.nochim, taxa, file = "dada2_output.Rdata")

load("dada2_output.Rdata")

# Make OTU - sequence - taxa table for later
colnames(seqtab.nochim)[c(1:3)]
colnames(taxa)[c(1:10)]
rownames(taxa)[c(1)]

seqtab.trans <- as.data.frame(t(seqtab.nochim))
head(seqtab.trans)

# create short OTU IDs
seqtab.trans$ids <- paste0("OTU", seq(1, length(colnames(seqtab.nochim))))
ids <-paste0("OTU", seq(1, length(colnames(seqtab.nochim))))
head(seqtab.trans)

# merge with taxa
otu_taxa_seq <- merge(seqtab.trans, taxa, by = 0)
# coorelate otu with taxa (removing seq info)
otu_taxa <- select(otu_taxa_seq, "ids", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# make sure it worked
head(otu_taxa)
write.csv(otu_taxa, file="taxa.csv")
# relevel sam_info so that the order is Early -> Mid -> Late, instead of alphabetical (Early, Late, Mid)
head(sam_info)
levels(sam_info$Time)
sam_info$Time <- factor(sam_info$Time, levels = c("Early", "Mid", "Late"))
levels(sam_info$Time)

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sam_info), 
               tax_table(taxa))
ps

# Let's look at what we have for phyloseq
otu_table(ps)[1:5, 1:5]

# yikes, nasty column names. Change that.
colnames(seqtab.nochim) <- ids
head(seqtab.nochim)[2,]

# Replace taxa names in the phyloseq object
taxa_names(ps) <- ids

# Try again...
otu_table(ps)[1:5, 1:5]
# Much better! Rows = samples. Columns = OTUs. Abundances (counts) fill the cells.

# What are the sample varaibles in our 'ps' object?
sample_variables(ps)
# "SampleID" "Site"     "Time"     "techRep"  "siteType"

# taxonomic ranks are "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"  
rank_names(ps)

# How many samples do we have?
nsamples(ps)
# 72, correct

# What does the taxonomic table look like?
tax_table(ps)[1:5, 1:6]
# Each OTU is associated with six levels of taxonomy (KPCOFG --- no species)

# Separate data and conds objects into two groups: STRI only (AM+midday+PM) and All sites (Midday only)
head(sam_info)

# Subset the `ps` object into STRI only (all timepoints) and all site (midday only)
psSTRI <- subset_samples(ps, Site=="STRIPoint")
psMid <- subset_samples(ps, Time=="Mid")

# Save
save(sam_info, seqtab.nochim, taxa, ps, psSTRI, psMid, otu_taxa, otu_taxa_seq, file="startHere4vegan.Rdata")

# Count the copepods ---------
load("startHere4vegan.Rdata")
head(otu_taxa_seq)

copepod_order <- c("Calanoida", "Canuelloida", "Cyclopoida", "Gelyelloida", "Harpacticoida", "Misophrioida", "Monstrilloida", "Mormonilloida", "Platycopioida", "Siphonostomatoida")
copepods <- otu_taxa_seq %>% filter(Order %in% copepod_order) %>% select(-c(Row.names,ids,Kingdom,Phylum,Class,Family,Genus))
summary(copepods$Order)
copepod_counts <- copepods %>% select(-Order) %>% sum()
copepod_counts #1,281,071 counts for copepods
total_counts <- otu_taxa_seq %>% select(-c(Row.names,ids,Kingdom,Phylum,Class,Family,Order,Genus)) %>% sum()
total_counts #3,675,476 total counts
(copepod_counts/total_counts)*100
# 34.9 percent copepod counts

just_calanoida <- otu_taxa_seq %>% filter(Order %in% c("Calanoida")) %>% select(-c(Row.names,ids,Kingdom,Phylum,Class,Family,Genus))
just_calanoida <- just_calanoida %>% select(-Order) %>% sum()
just_calanoida #1,259,402
((just_calanoida)/(total_counts))*100
# 34.265

most_order <- otu_taxa_seq %>% filter(!is.na(Order)) %>% select(-c(Row.names,ids,Kingdom,Phylum,Class,Family,Genus)) %>%
  gather(sample, count, P04.1:P39.2, factor_key=TRUE) %>% select(-sample) %>%
  group_by(Order) %>% summarize(sum = sum(count)) %>%
  arrange(sum)
head(most_order)
tail(most_order)

# how many Bacillariophytina?
just_Bacillariophytina <- otu_taxa_seq %>% filter(Order %in% c("Bacillariophytina")) %>% select(-c(Row.names,ids,Kingdom,Phylum,Class,Family,Genus))
just_Bacillariophytina <- just_Bacillariophytina %>% select(-Order) %>% sum()
((just_Bacillariophytina)/(total_counts))*100

# make pretty for posting
head(otu_taxa_seq)
otu_table <- otu_taxa_seq %>% select(-Row.names) %>% mutate(classification = paste(Kingdom,Phylum,Class,Order,Family,Genus,sep="_")) %>%
  select(ids, classification,P04.1:P39.2)
head(otu_table)
write.table(otu_table, file="counts_with_classification.txt",sep="\t",row.names=F,quote=F)

# START HERE FOR VEGAN ---------

# Load in data
load("startHere4vegan.Rdata")

write.table(otu_taxa,file="ASV_Taxonomic_assignment.txt",quote=F,sep="\t")

# Sum technical replicates

# Tidy up "sam_info
sam_info <- sam_info %>% 
  mutate(sitecode = ifelse(Site=="STRIPoint", "ST_in", 
                       ifelse(Site=="BastimentosN", "BN_off", 
                              ifelse(Site=="BastimentosS", "BS_off", 
                                     ifelse(Site=="Cristobal", "CR_in", 
                                            ifelse(Site=="DragoMar", "DM_off", 
                                                   ifelse(Site=="PopaIsland","PI_off",
                                                          ifelse(Site=="PuntaLaurel","PL_in", "PD_in"))))))), 
         tow = substr(SampleID,2,3),
         samID = as.factor(paste(Time,tow,sitecode,sep="_")))

head(sam_info)
summary(sam_info)

summary(sam_info$Site)
summary(sam_info$Time)
# Early   Mid  Late 
# 12    48    12 

nrow(sam_info[grep("Mid_[0-9][0-9]_ST_in",sam_info$samID),]) #6 -- 3 sets of tech reps, matches pcoa
sam_info[grep("Mid_[0-9][0-9]_ST_in",sam_info$samID),]
nrow(sam_info[grep("Early_[0-9][0-9]_ST_in",sam_info$samID),]) #12
nrow(sam_info[grep("Late_[0-9][0-9]_ST_in",sam_info$samID),]) #12

old2new <- sam_info %>% 
  select(SampleID, samID) %>%
  column_to_rownames("SampleID")
head(old2new)

# sum tech reps
alldat0 <- as.data.frame(seqtab.nochim)
head(alldat0)[c(1)]

alldat1 <- merge(alldat0,old2new,by=0)
head(alldat1[c(1,length(alldat1))])
alldat1[grep("Late",alldat1$samID),][c(1,length(alldat1))]

alldat <- alldat1 %>% 
  group_by(samID) %>% # group by sample ID (will group technical replicates)
  mutate_all(function(x) as.numeric(as.character(x))) %>% # convert to numeric so "sum" will work
  summarise_all(funs(sum)) %>% # add up the counts
  select(-Row.names) # it makes a column called "Row.names", remove it
head(alldat)[1:10]

# manual spot check to make sure that worked
sum(alldat1[alldat1$samID=="Early_26_ST_in",][5])==alldat[alldat$samID=="Early_26_ST_in",][5]
sum(alldat1[alldat1$samID=="Early_20_ST_in",][10])==alldat[alldat$samID=="Early_20_ST_in",][10]

# subset alldat for samples in "mid" only(run to make goods only Mid time pt)----
summary(sam_info)

Mid_samples <- sam_info %>%
  filter(Time=="Mid")
summary(Mid_samples)
head(Mid_samples)
nrow(Mid_samples)
# 3 tows * 8 sites * 2 library replicates = 48

alldat.Mid <- alldat %>%
  filter(samID %in% Mid_samples$samID) %>%
  rename(sample = samID) # MCMC.OTU needs this column to be named "sample"
head(alldat.Mid[c(1:3)])

# subset alldat for samples in "STRI" only(run to make goods only STRI)----
STRI_samples <- sam_info %>%
  filter(Site=="STRIPoint")
summary(STRI_samples)

alldat.STRI <- alldat %>%
  filter(samID %in% STRI_samples$samID) %>%
  rename(sample = samID) # MCMC.OTU needs this column to be named "sample"
alldat.STRI$sample
# [1] Early_19_ST_in Early_20_ST_in Early_21_ST_in Early_25_ST_in Early_26_ST_in Early_27_ST_in Late_22_ST_in 
# [8] Late_23_ST_in  Late_24_ST_in  Late_31_ST_in  Late_32_ST_in  Late_33_ST_in  Mid_37_ST_in   Mid_38_ST_in
# [15] Mid_39_ST_in

head(alldat.STRI[c(1:3)])
table(sapply(alldat.STRI,is.numeric))

# purging under-sequenced samples; 
# change what is being purged in alldat. depending on site or time 

goods.STRI <- purgeOutliers(alldat.STRI,
                       count.columns = c(2:ncol(alldat.STRI)),
                       sampleZcut = (-2.5),
                       otu.cut = 0.00001,
                       zero.cut = 0.01)
summary(goods.STRI)[,1:6]
goods.STRI$sample
# [1] Early_19_ST_in Early_20_ST_in Early_21_ST_in Early_25_ST_in Early_26_ST_in Early_27_ST_in Late_22_ST_in # [8] Late_23_ST_in  Late_24_ST_in  Late_31_ST_in  Late_32_ST_in  Late_33_ST_in  Mid_37_ST_in   Mid_38_ST_in
# [15] Mid_39_ST_in   

goods.Mid <- purgeOutliers(alldat.Mid,
                            count.columns = c(2:ncol(alldat.Mid)),
                            otu.cut = 0.00001,
                            zero.cut = 0.01)
summary(goods.Mid)[,1:6]
# "samples with counts below z-score -2.5 :"
# [1] "Mid_07_BS_off"

# creating a log-transfromed normalized dataset for PCoA:
goods.log.STRI <- logLin(data = goods.STRI,
                    count.columns = 2:length(names(goods.STRI)))
summary(goods.log.STRI)[,1:6]

goods.log.Mid <- logLin(data = goods.Mid,
                         count.columns = 2:length(names(goods.Mid)))
summary(goods.log.Mid)[,1:6]

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
goods.dist.STRI <- vegdist(goods.log.STRI, method = "manhattan")
goods.pcoa.STRI <- pcoa(goods.dist.STRI)

goods.dist.Mid <- vegdist(goods.log.Mid, method = "manhattan")
goods.pcoa.Mid <- pcoa(goods.dist.Mid)

# make conditions
head(sam_info)
conditions.STRI <- sam_info %>%
  select(Site,Time,siteType,sitecode,samID) %>%
  filter(samID %in% goods.STRI$sample) %>%
  distinct()
head(conditions.STRI)
dim(conditions.STRI)
dim(goods.STRI)

conditions.Mid <- sam_info %>%
  select(Site,Time,siteType,sitecode,samID) %>%
  filter(samID %in% goods.Mid$sample) %>%
  distinct()
head(conditions.Mid)
dim(conditions.Mid)
dim(goods.Mid)
goods.Mid[1:2,1:2]

# does sample order match?
table(conditions.STRI$samID == goods.STRI$sample)
# Reorder the samples conds table to match the order in the counts matrix
conditions.STRI <- conditions.STRI[match(goods.STRI$sample, conditions.STRI$samID),]
table(conditions.STRI$samID == goods.STRI$sample)
# fixed!
table(conditions.Mid$samID == goods.Mid$sample)
# good

# Analysis for PCoA p-value----

# set vectors for analysis 
siteType <- as.vector(conditions.Mid$siteType)
Site <- as.vector(conditions.Mid$Site)
time <- as.factor(conditions.STRI$Time)
# siteShape <- ifelse(site=="Cristobal", 15, ifelse(site=="DragoMar", 1, ifelse(site=="PopaIsland", 17, ifelse(site=="PuntaDonato", 19, 5))))
# siteTypeColor <- ifelse(siteType=="inshore","salmon", "royalblue4")
# timeColor <- ifelse(time=="Early","orange",ifelse(time=="Mid","green","blue"))

# create the df of sites 
metaData.Mid <- data.frame(cbind(siteType, Site))
head(metaData.Mid)

# check order one more time...
cbind(as.character(goods.Mid$sample), siteType)

# analysis by site
# p = 0.001
set.seed(1)
adonis.Midsite <- adonis(goods.log.Mid~Site,metaData.Mid,distance="manhattan")
adonis.Midsite$aov.tab$`Pr(>F)`[1]

# analysis by siteType
# new = 0.015
set.seed(1)
adonis.MidsiteType <- adonis(goods.log.Mid~siteType,metaData.Mid,distance="manhattan")
adonis.MidsiteType$aov.tab$`Pr(>F)`[1]

# create df for stri
time <- as.data.frame(time)

# analysis for stri
# new = 0.001
set.seed(1)
adonis.STRI <- adonis(goods.log.STRI~time,time,distance="manhattan")
adonis.STRI$aov.tab$`Pr(>F)`[1]

# Save
save(goods.pcoa.STRI, goods.pcoa.Mid, conditions.Mid, conditions.STRI, adonis.Midsite, 
       adonis.MidsiteType, adonis.STRI, goods.pcoa.STRI, goods.dist.Mid, goods.dist.STRI,
       goods.log.STRI, goods.log.Mid, goods.STRI, goods.Mid, goods.log.Mid, goods.log.STRI, otu_taxa, 
     file="startgraphs.Rdata")

# START HERE FOR GRAPHS ---------

# Load in data
load("startgraphs.Rdata")

# plotting PCoA----
scores.STRI <- goods.pcoa.STRI$vectors
scores.Mid <- goods.pcoa.Mid$vectors
margin <- 0.01

# play around with these numbers
xaxis <- 1
yaxis <- 2
# PCoA for mid by site type
quartz()
plot(scores.Mid[,xaxis], scores.Mid[,2],type="n",
     xlim=c(min(scores.Mid[,xaxis])-margin,max(scores.Mid[,xaxis])+margin),
     ylim=c(min(scores.Mid[,2])-margin,max(scores.Mid[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("Axis", xaxis,"(", round(goods.pcoa.Mid$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("Axis", yaxis,"(", round(goods.pcoa.Mid$values$Relative_eig[yaxis]*100,1),"%)",sep="")) + 
  ordihull(scores.Mid,conditions.Mid$Site,label=F, draw = "polygon", col = c("royalblue4", "royalblue4", "salmon","royalblue4", "royalblue4", "salmon","salmon", "salmon", alpha = 255))
  points(scores.Mid[conditions.Mid$Site=="PuntaDonato",xaxis],scores.Mid[conditions.Mid$Site=="PuntaDonato",yaxis], col="salmon", pch=19) +
  points(scores.Mid[conditions.Mid$Site=="STRIPoint",xaxis],scores.Mid[conditions.Mid$Site=="STRIPoint",yaxis], col="salmon", pch=17) +  
  points(scores.Mid[conditions.Mid$Site=="Cristobal",xaxis],scores.Mid[conditions.Mid$Site=="Cristobal",yaxis], col="salmon", pch=15) +
  points(scores.Mid[conditions.Mid$Site=="PuntaLaurel",xaxis],scores.Mid[conditions.Mid$Site=="PuntaLaurel",yaxis], col="salmon", pch=18) + 
  points(scores.Mid[conditions.Mid$Site=="DragoMar",xaxis],scores.Mid[conditions.Mid$Site=="DragoMar",yaxis], col="royalblue4", pch=19) +
  points(scores.Mid[conditions.Mid$Site=="BastimentosN",xaxis],scores.Mid[conditions.Mid$Site=="BastimentosN",yaxis], col="royalblue4", pch=17) +
  points(scores.Mid[conditions.Mid$Site=="BastimentosS",xaxis],scores.Mid[conditions.Mid$Site=="BastimentosS",yaxis], col="royalblue4", pch=15) +
  points(scores.Mid[conditions.Mid$Site=="PopaIsland",xaxis],scores.Mid[conditions.Mid$Site=="PopaIsland",yaxis], col="royalblue4", pch=18) 
 # legend of sites 
  legend(80,120, 
         c("PD_in","ST_in","CR_in","PL_in"),
           pch=c(19,17,15,18), 
         col=c("salmon", "salmon","salmon","salmon"), 
         cex=1, bty = "n")
  legend(-100,120,
         c("DM_off","BN_off","BS_off","PI_off"),
         pch=c(19,17,15,18),
         col=c("royalblue4","royalblue4","royalblue4","royalblue4"), cex=1, bty = 'n')
 #insert p value 
   legend("bottomright", inset=.02, paste("p = ",adonis.Midsite$aov.tab$`Pr(>F)`[1], sep=" "), cex=1, bty='n')  
  
  
# PCoA midday hull by inshore/offshore
# PCoA for mid by site type
quartz()
plot(scores.Mid[,xaxis], scores.Mid[,2],type="n",
       xlim=c(min(scores.Mid[,xaxis])-margin,max(scores.Mid[,xaxis])+margin),
       ylim=c(min(scores.Mid[,2])-margin,max(scores.Mid[,2])+margin),
       mgp=c(2.3,1,0),
       xlab=paste("Axis", xaxis,"(", round(goods.pcoa.Mid$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
       ylab=paste("Axis", yaxis,"(", round(goods.pcoa.Mid$values$Relative_eig[yaxis]*100,1),"%)",sep=""),) +
    ordihull(scores.Mid,conditions.Mid$siteType,label=F, draw = "polygon", col = c( "salmon","royalblue4", alpha = 255))
  points(scores.Mid[conditions.Mid$Site=="PuntaDonato",xaxis],scores.Mid[conditions.Mid$Site=="PuntaDonato",yaxis], col="salmon", pch=19) +
    points(scores.Mid[conditions.Mid$Site=="STRIPoint",xaxis],scores.Mid[conditions.Mid$Site=="STRIPoint",yaxis], col="salmon", pch=19) +  
    points(scores.Mid[conditions.Mid$Site=="Cristobal",xaxis],scores.Mid[conditions.Mid$Site=="Cristobal",yaxis], col="salmon", pch=19) +
    points(scores.Mid[conditions.Mid$Site=="PuntaLaurel",xaxis],scores.Mid[conditions.Mid$Site=="PuntaLaurel",yaxis], col="salmon", pch=19) + 
    points(scores.Mid[conditions.Mid$Site=="DragoMar",xaxis],scores.Mid[conditions.Mid$Site=="DragoMar",yaxis], col="royalblue4", pch=19) +
    points(scores.Mid[conditions.Mid$Site=="BastimentosN",xaxis],scores.Mid[conditions.Mid$Site=="BastimentosN",yaxis], col="royalblue4", pch=19) +
    points(scores.Mid[conditions.Mid$Site=="BastimentosS",xaxis],scores.Mid[conditions.Mid$Site=="BastimentosS",yaxis], col="royalblue4", pch=19) +
    points(scores.Mid[conditions.Mid$Site=="PopaIsland",xaxis],scores.Mid[conditions.Mid$Site=="PopaIsland",yaxis], col="royalblue4", pch=19)
  legend(60,120, 
         c("Inshore","Offshore"),
         pch=c(19,19), 
         col=c("salmon","royalblue4"), cex=1, bty = "n")
  legend("topleft", inset=.02, paste("p = ",adonis.MidsiteType$aov.tab$`Pr(>F)`[1], sep=" "), cex=1, bty='n')  

# PCOA STRIPoint by time of day
quartz()
plot(scores.STRI[,xaxis], scores.STRI[,2],type="n",
       xlim=c(min(scores.STRI[,xaxis])-margin,max(scores.STRI[,xaxis])+margin),
       ylim=c(min(scores.STRI[,2])-margin,max(scores.STRI[,2])+margin),
       mgp=c(2.3,1,0),
       xlab=paste("Axis", xaxis,"(", round(goods.pcoa.STRI$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
       ylab=paste("Axis", yaxis,"(", round(goods.pcoa.STRI$values$Relative_eig[yaxis]*100,1),"%)",sep=""),las=2)
      ordihull(scores.STRI,conditions.STRI$Time,label=F, 
               draw = "polygon", col = c("orange", "salmon","green", alpha = 255))
  points(scores.STRI[conditions.STRI$Time=="Early",xaxis],scores.STRI[conditions.STRI$Time=="Early",yaxis], col="orange", pch=17) +
  points(scores.STRI[conditions.STRI$Time=="Mid",xaxis],scores.STRI[conditions.STRI$Time=="Mid",yaxis], col="salmon", pch=17) +
  points(scores.STRI[conditions.STRI$Time=="Late",xaxis],scores.STRI[conditions.STRI$Time=="Late",yaxis], col="green", pch=17) 
  legend(100,100, c("Early", "Mid", "Late"), pch = c(17), col=c("orange","salmon","green"), bty="n", cex=1)
  legend(80,-50, paste("p = ",adonis.STRI$aov.tab$`Pr(>F)`[1], sep=" "), cex=1, bty='n')  
  
  
# analysis of multivariate homogeneity of group dispersions (variances) -----

# centroid distance for siteType
groups_siteType <- factor(conditions.Mid$siteType)
disp_siteType <- betadisper(goods.dist.Mid, groups_siteType, type = "centroid")
disp_siteType
# Average distance to centroid:
# inshore offshore 
# 164.9    149.4

# significance of centroid distance for siteType
anovadisp_siteType <- anova(disp_siteType)
# Response: Distances
# Df Sum Sq Mean Sq F value Pr(>F)
# Groups     1  230.1  230.11    1.01 0.3263
# Residuals 21 4784.6  227.84

# significance between each of the distances for siteType 
tukey_siteType <- TukeyHSD(disp_siteType, which = "group", ordered = FALSE, conf.level = 0.95)
# $group
# diff       lwr      upr     p adj
# offshore-inshore -15.43778 -30.48454 -0.391012 0.0448268

# centroid distance for site
groups_site <- factor(conditions.Mid$Site)
disp_site <- betadisper(goods.dist.Mid, groups_site, type = "centroid")
disp_site
# Average distance to centroid:
# BastimentosN BastimentosS    Cristobal     DragoMar   PopaIsland  PuntaDonato  PuntaLaurel    STRIPoint 
# 111.05        94.79       136.20       147.03       108.44       141.32       137.44       115.59 

# significance of centroid distance for site
anova_site <- anova(disp_site)
anova_site
# Response: Distances
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Groups     7 6591.4  941.62  4.9804 0.004415 **
#   Residuals 15 2836.0  189.07

# significance between each of the distances for site
tukey_site <- TukeyHSD(disp_site, which = "group", ordered = FALSE, conf.level = 0.95)
tukey_site

# centroid distance for time of day 
groups_time <- factor(conditions.STRI$Time)
disp_time <- betadisper(goods.dist.STRI, groups_time, type = "centroid")
disp_time
# Average distance to centroid:
#   Early   Mid  Late 
# 143.1 125.0 164.5

# significance of centroid distance for time of day
anova_time <- anova(disp_time)
anova_time
# Response: Distances
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Groups     2 3374.0 1686.99  8.0995 0.005939 **
# Residuals 12 2499.4  208.28         

# significance between each of the distances for time of day 
tukey_time <- TukeyHSD(disp_time, which = "group", ordered = FALSE, conf.level = 0.95)
tukey_time

# $group
# diff         lwr       upr     p adj
# Mid-Early  -18.14895 -45.3744794  9.076573 0.2180725
# Late-Early  21.41350  -0.8160528 43.643046 0.0593529
# Late-Mid    39.56245  12.3369230 66.787976 0.0057639

# all three plots
quartz()
par(mfrow=c(1,3))
boxplot(disp_time, ylab = "Distance to Centroid",
        col = c("gold","orange","blue"))
legend("bottomright", paste("anova p =", round(anova_time$`Pr(>F)`[1],2), sep=" "), bty="n")
boxplot(disp_site, ylab = "Distance to Centroid",
        col = c("royalblue4","royalblue4","salmon","royalblue4","royalblue4","salmon","salmon","salmon"),
        main = "Midday Only- Dispersion by Site")
legend("bottomright", paste("anova p =", round(anova_site$`Pr(>F)`[1],2), sep=" "), bty="n")
boxplot(disp_siteType, ylab = "Distance to Centroid", col = c("salmon","royalblue4"),
        main = "Midday Only- Dispersion by Reef Zone")
legend("bottomright", paste("p.adj =", round(tukey_siteType$group[4],2), sep=" "), bty="n")


# Shannon diversity -----
# shannon by time   
sample <- as.vector(goods.STRI$sample)
time <- as.vector(conditions.STRI$Time)
shannonSTRI <- diversity(goods.log.STRI, "shannon")
df.shannon.STRI <- data.frame(sample,time,shannonSTRI)
df.shannon.STRI$time <- factor(df.shannon.STRI$time, levels = c("Early", "Mid", "Late"))

# stats for shannon by time p= 0.489
aovshannonstri <- aov(shannonSTRI~time, data=df.shannon.STRI)
summary(aovshannonstri)

# check for normal data distribution of shannon time
par(mfrow=c(2,2))
plot(aovshannonstri)

# shannon by site 
sample <- as.vector(goods.Mid$sample)
Site <- as.vector(conditions.Mid$sitecode)
shannonMid <- diversity(goods.log.Mid, "shannon")
df.shannon.Mid <- data.frame(sample,Site,shannonMid)
levels(df.shannon.Mid$Site) 
df.shannon.Mid$Site <- factor(df.shannon.Mid$Site, levels(df.shannon.Mid$Site)[c(5,8,3,7,1,2,4,6)])

# stats for shannon by site
aovshannonsite <- aov(shannonMid~Site, data=df.shannon.Mid)
summary(aovshannonsite)

# check for normal data distribution of shannon site
par(mfrow=c(2,2))
plot(aovshannonsite)

#mean for each Site 
aggregate(. ~ Site, df.shannon.Mid[-1], mean)

# shannon by siteType 
sample <- as.vector(goods.Mid$sample)
siteType <- as.vector(conditions.Mid$siteType)
shannonMid <- diversity(goods.log.Mid, "shannon")
df.shannon.Mid <- data.frame(sample,Site,shannonMid)
levels(df.shannon.Mid$Site) 
df.shannon.Mid$Site <- factor(df.shannon.Mid$Site, levels(df.shannon.Mid$Site)[c(5,8,3,7,1,2,4,6)])

# stats for shannon by siteType
aovsiteType <- aov(shannonMid~siteType, data=df.shannon.Mid)
summary(aovsiteType)
# Df   Sum Sq   Mean Sq F value Pr(>F)
# siteType     1 0.000102 1.020e-04   1.683  0.209
# Residuals   21 0.001273 6.062e-05 

# check for normal data distribution of shannon siteType
par(mfrow=c(2,2))
plot(aovsiteType)

# Simpson-------
# simpson by time 
sample <- as.vector(goods.STRI$sample)
time <- as.vector(conditions.STRI$Time)
simpsonSTRI <- diversity(goods.log.STRI, "simpson")
df.simpson.STRI <- data.frame(sample,time,simpsonSTRI,3)
df.simpson.STRI$time <- factor(df.simpson.STRI$time, levels = c("Early", "Mid", "Late"))

# stats for simpson by time
aovsimpsontime <- aov(simpsonSTRI~time, data=df.simpson.STRI)
summary(aovsimpsontime)
# Df    Sum Sq   Mean Sq F value Pr(>F)
# time         2 4.910e-09 2.455e-09   0.669   0.53
# Residuals   12 4.403e-08 3.669e-09

# check for normal data distribution of simpson time
par(mfrow=c(2,2))
plot(aovsimpsontime)

# Plot Shannon and Simpson
quartz()
par(mfrow=c(1,2))
boxplot(shannonSTRI~time,data=df.shannon.STRI, col = c("orange","salmon","green"), 
        main = "Shannon's Index of Diversity (H)", cex.lab = 1.5, cex.axis = 1.5)
legend("bottomright", paste("p =", round(summary(aovshannonstri)[[1]][["Pr(>F)"]][1],3), 
                            sep = " "), cex=1.5, bty='n')
boxplot(simpsonSTRI~time,data=df.simpson.STRI, col = c("orange","salmon","green"),las=1,
        main = "Simpson's Index of Diversity (1-D)", cex.lab = 1.5, cex.axis = 1.5)
legend("bottomright", paste("p =", round(summary(aovsimpsontime)[[1]][["Pr(>F)"]][1],3)), 
       cex=1.5, bty='n')

# simpson by site
sample <- as.vector(goods.Mid$sample)
Site <- as.vector(conditions.Mid$sitecode)
simpsonMid <- diversity(goods.log.Mid, "simpson")
df.simpson.Mid <- data.frame(sample,Site,simpsonMid)
levels(df.simpson.Mid$Site) 
df.simpson.Mid$Site <- factor(df.simpson.Mid$Site, levels(df.simpson.Mid$Site)[c(5,8,3,7,1,2,4,6)])
# stats for simpson by site
aovsimpsonsite <- aov(simpsonMid~Site, data=df.simpson.Mid)
summary(aovsimpsonsite)

# Plot Shannon and Simpson
quartz()
par(mfrow=c(1,2))
boxplot(shannonMid~Site,data=df.shannon.Mid, col =c("salmon","salmon","salmon","salmon","royalblue4", "royalblue4", "royalblue4", "royalblue4"), las=2,
        main = "Shannon's Index of Diversity (H)", cex.lab = 1.5, cex.axis = 1.5)
legend("bottomleft", inset=.02, paste("p =", round(summary(aovshannonsite)[[1]][["Pr(>F)"]][1],3)), 
       cex=1.5, bty='n')
boxplot(simpsonMid~Site,data=df.simpson.Mid,col =c("salmon","salmon","salmon","salmon","royalblue4", "royalblue4", "royalblue4", "royalblue4"), las=2,
        main = "Simpson's Index of Diversity (1-D)", cex.lab = 1.5, cex.axis = 1.5)
legend("bottomleft", paste("p = ", round(summary(aovsimpsonsite)[[1]][["Pr(>F)"]][1],3)), cex=1.5, bty='n')

# #mean for each Site 
aggregate(. ~ Site, df.simpson.Mid[-1], mean)

# simpson by siteType
sample <- as.vector(goods.Mid$sample)
Site <- as.vector(conditions.Mid$sitecode)
simpsonMid <- diversity(goods.log.Mid, "simpson")
df.simpson.Mid <- data.frame(sample,Site,simpsonMid)
levels(df.simpson.Mid$Site) 
df.simpson.Mid$Site <- factor(df.simpson.Mid$Site, levels(df.simpson.Mid$Site)[c(5,8,3,7,1,2,4,6)])

# stats for simpson by siteType p= 0.19 
aovsimpsonsiteType <- aov(simpsonMid~siteType, data=df.simpson.Mid)
summary(aovsimpsonsiteType)

# check for normal data distribution of simpson siteType
par(mfrow=c(2,2))
plot(aovsimpsonsiteType)
dev.off()

# Plot Shannon and Simpson
quartz()
par(mfrow=c(1,2))
boxplot(shannonMid~siteType, data=df.shannon.Mid, col =c("salmon", "royalblue4"),las=1,
        main = "Shannon's Index of Diversity (H)", cex.lab = 1.5, cex.axis = 1.5)
legend("bottomleft", paste("p =", round(summary(aovsiteType)[[1]][["Pr(>F)"]][1],3)), 
       cex=1.5, bty='n')
boxplot(simpsonMid~siteType,data=df.simpson.Mid,col =c("salmon","royalblue4"), las=1,
        main = "Simpson's Index of Diversity (1-D)", cex.lab = 1.5, cex.axis = 1.5)
legend("bottomleft", paste("p = ", round(summary(aovsimpsonsiteType)[[1]][["Pr(>F)"]][1],3)), cex=1.5, bty='n')


# MCMC OTU analysis on site type midday only -----
  
# reformat data for mcmc.otu
goods.mcmc.STRI <- merge(goods.STRI,conditions.STRI, by.x = "sample", by.y="samID")
names(goods.mcmc.STRI)[c(1:3)]
names(goods.mcmc.STRI)[c((length(goods.mcmc.STRI)-5):length(goods.mcmc.STRI))]

goods.mcmc.Mid <- merge(goods.Mid,conditions.Mid, by.x = "sample", by.y="samID")
names(goods.mcmc.Mid)[c(1:3)]
names(goods.mcmc.Mid)[c((length(goods.mcmc.Mid)-5):length(goods.mcmc.Mid))]

# stacking the data table
names(goods.mcmc.STRI[c(497:501)])
gs.STRI <- otuStack(goods.mcmc.STRI,
    count.columns=c(2:(ncol(goods.mcmc.STRI)-4)),
    condition.columns = c(1,498:501))
head(gs.STRI)

names(goods.mcmc.Mid[c(475:483)])
gs.Mid <- otuStack(goods.mcmc.Mid,
                    count.columns=c(2:(ncol(goods.mcmc.Mid)-4)),
                    condition.columns = c(1,480:483))
head(gs.Mid)
write.csv(gs.STRI, file="STRI_counts.csv")

# fitting the model
set.seed(1)
mm.STRI <- mcmc.otu(
  fixed = "Time",
  data = gs.STRI
)

summary(mm.STRI)

set.seed(1)
mm.Mid <- mcmc.otu(
  fixed = "siteType",
  data = gs.Mid
)

summary(mm.Mid)

# selecting the OTUs that were modeled reliably
acpass.STRI <- otuByAutocorr(mm.STRI,gs.STRI)
head(acpass.STRI)

acpass.Mid <- otuByAutocorr(mm.Mid,gs.Mid)
head(acpass.Mid)

# calculating effect sizes and p-values:
ss.STRI <- OTUsummary(mm.STRI,gs.STRI,summ.plot=FALSE)
ss.Mid <- OTUsummary(mm.Mid,gs.Mid,summ.plot=FALSE)

# correcting for mutliple comparisons (FDR)
ss.STRI <- padjustOTU(ss.STRI)
ss.Mid <- padjustOTU(ss.Mid)

# getting significatly changing OTUs (FDR<0.05)
sigs.STRI <- signifOTU(ss.STRI, p.cutoff = 0.005)
length(sigs.STRI) #16 when cutoff is 0.005
sigs.STRI

sigs.Mid <- signifOTU(ss.Mid, p.cutoff = 0.05)
length(sigs.Mid) #12 when cutoff is 0.05
sigs.Mid

# Filter stack for significant OTUs
sig16_STRI <- gs.STRI %>%
  filter(otu %in% sigs.STRI)

sig16_Mid <- gs.Mid %>%
  filter(otu %in% sigs.Mid)

# Add taxonomy
gs_order_STRI <- merge(sig16_STRI, otu_taxa, by.x=2, by.y=1)
head(gs_order_STRI)
gs_order_Mid <- merge(sig16_Mid, otu_taxa, by.x=2, by.y=1)
head(gs_order_Mid)

# SAVE/LOAD MCMC
save(sam_info, goods.STRI, goods.Mid, goods.log.STRI, goods.log.Mid,
     mm.STRI, mm.Mid,
     sigs.STRI, sigs.Mid, ss.STRI, ss.Mid,
     gs.STRI, gs.Mid,
     gs_order_STRI, gs_order_Mid,
     sig16_Mid, sig16_STRI,
     file="mcmcanalysis.Rdata")

load("mcmcanalysis.Rdata")

# Make heatmap -------------------------
library(pheatmap)

# Midday only 
# reformat for making the heatmap
counts_Mid <- gs_order_Mid %>%
  mutate(taxa = paste(Kingdom,Phylum,Class,Order,Family,Genus,sep="_")) %>% # make full taxa name
  mutate(taxa_otu = paste(otu,taxa,sep="_")) %>% # make unique taxa names
  select(otu,count,sample,taxa_otu) %>% # only keep columns you need
  spread(sample,count) %>% # format for heatmap
  filter(otu %in% sig16_Mid$otu) %>%# filter for significant
  column_to_rownames("taxa_otu") %>% # make the unique taxa the rownames
  select(-"otu") # get rid of otu column

head(counts_Mid)
nrow(counts_Mid)

# change order for plotting
names(counts_Mid)
order <- c(grep("in",names(counts_Mid)),
           grep("off",names(counts_Mid)))
counts_Mid <- counts_Mid[c(order)]
head(counts_Mid)


# Make color scale
heat.colors <- colorRampPalette(rev(c("red","white","gold")),bias=1)(100)

# Plot heatmap
pheatmap(as.matrix(log(counts_Mid+0.001)),
         scale = "row", cex=0.9,
         color = heat.colors, border_color=NA,
         clustering_distance_rows="correlation", cluster_cols=F)

# STRI only 
# reformat for making the heatmap
counts_STRI <- gs_order_STRI %>%
  mutate(taxa = paste(Kingdom,Phylum,Class,Order,Family,Genus,sep="_")) %>% # make full taxa name
  mutate(taxa_otu = paste(otu,taxa,sep="_")) %>% # make unique taxa names
  select(otu,count,sample,taxa_otu) %>% # only keep columns you need
  spread(sample,count) %>% # format for heatmap
  filter(otu %in% sig16_STRI$otu) %>%# filter for significant
  column_to_rownames("taxa_otu") %>% # make the unique taxa the rownames
  select(-"otu") # get rid of otu column

head(counts_STRI)
nrow(counts_STRI)

# change order for plotting
names(counts_STRI)
order <- c(grep("Early",names(counts_STRI)),
           grep("Mid",names(counts_STRI)),
           grep("Late",names(counts_STRI)))
counts_STRI <- counts_STRI[c(order)]

# Make color scale
heat.colors <- colorRampPalette(rev(c("red","white","gold")),bias=1)(100)

# Plot heatmap
pheatmap(as.matrix(counts_STRI),
         scale = "row", cex=0.9,
         color = heat.colors, border_color=NA,
         clustering_distance_rows="correlation", cluster_cols=F)

# DESeq2 by site (midday only)--------
load("mcmcanalysis.Rdata")
load("startgraphs.Rdata")

head(gs.Mid)
tail(gs.Mid)

dat <- gs.Mid
dat1 <- gs.Mid %>%
  dplyr::select(-Site, -Time, -siteType, -sitecode) %>%
  spread(sample, count) %>%
  filter(!otu == "summ")
dim(dat1)
# 478 24
head(dat1[1:10,1:10])

# replace OTU names with taxonomy names
head(otu_taxa)
taxanames <- otu_taxa %>% mutate(fullID = paste(ids, Phylum, Class, Order, Family, Genus)) %>% select(ids, fullID)
head(taxanames)
dat1 <- merge(dat1, taxanames, by.x="otu",by.y="ids")
dim(dat1)
#478 25
head(dat1[1:10,c(1:10,25)])

row.names(dat1) <- dat1$fullID
dat1$otu <- NULL
dat1$fullID <- NULL
head(dat1) 
dim(dat1)
# 478 23

countData <- dat1
head(countData)
length(countData[,1])
#478

totalCounts <- colSums(countData)
mean(totalCounts) #98385.83
min(totalCounts) #63926
max(totalCounts) #132301

# make conditions data for DESeq
head(dat1)
names(dat1)
site <- rep("site",length(dat1))
site[grep("in",names(dat1))] <- "in"
site[grep("off",names(dat1))] <- "off"
site
colData <- data.frame(site)

# perform DESeq2
deds <- DESeqDataSetFromMatrix(countData=countData, 
                               colData=colData, 
                               design=~site)
dds <- DESeq(deds)

#log2 fold change (MLE): site off vs in 
res <- results(dds)
head(res)
table(res$padj<0.1)
# TRUE = 18

as.data.frame(res) %>% filter(padj<0.1)

plotMA(res, main="offshore vs Inshore")  
write.table(res, file="InOff_DE.txt", quote=F, sep="\t")

# rlog transformation for plotting
rlog <- rlogTransformation(dds, blind=TRUE) 
rld <- assay(rlog)
head(rld)
colnames(rld) <- paste(colnames(countData))
head(rld)
length(rld[,1]) #478
valIn <- cbind(res$pvalue, res$padj)
colnames(valIn)=c("pval", "padj")
length(valIn[,1])
table(complete.cases(valIn))
rldpvals <- as.data.frame(cbind(rld,valIn))
head(rldpvals)

# make heatmap
rldpvals$ids <- row.names(rldpvals)
cutoff <- 0.1 # FDR cutoff
conds <- rldpvals %>% filter(padj < cutoff) %>% select(-pval,-padj) %>% column_to_rownames("ids")
head(conds)
length(conds[,1])

means <- apply(conds,1,mean) # means of rows
head(means)
explc <- conds-means # subtracting them
head(explc)

ann <- data.frame(site=colData$site)
rownames(ann) <- names(explc)
Var1 <- c("royalblue4", "salmon")
names(Var1) <- c("off", "in")
anno_colors <- list(site = Var1)
paletteLength <- 100
myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)

myBreaks <- c(seq(min(explc), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(explc)/paletteLength, max(explc), length.out=floor(paletteLength/2)))

quartz()
pheatmap(as.matrix(explc),annotation_col=ann,annotation_colors=anno_colors,
         cex=1,color=myColor, breaks=myBreaks,border_color=NA,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation", show_rownames=T)

# DESeq2 by time of day  --------
dat <- gs.STRI
dat1 <- gs.STRI %>%
  dplyr::select(-Site, -Time, -siteType, -sitecode) %>%
  spread(sample, count) %>%
  filter(!otu == "summ")
dim(dat1)
# 496 14
head(dat1[1:10,1:10])

# replace OTU names with taxonomy names
head(otu_taxa)
taxanames <- otu_taxa %>% mutate(fullID = paste(ids, Phylum, Class, Order, Family, Genus)) %>% select(ids, fullID)
head(taxanames)
dat1 <- merge(dat1, taxanames, by.x="otu",by.y="ids")
dim(dat1)
#496 15
head(dat1[1:10,c(1:10,15)])

row.names(dat1) <- dat1$fullID
dat1$otu <- NULL
dat1$fullID <- NULL
head(dat1) 
dim(dat1)

countData <- dat1
head(countData)
length(countData[,1])
#496

totalCounts <- colSums(countData)
mean(totalCounts) #128503.5
min(totalCounts) #65959
max(totalCounts) #219946

# make conditions data for DESeq
head(dat1)
time <- rep("time",length(dat1))
time[grep("Early",names(dat1))] <- "Early"
time[grep("Mid",names(dat1))] <- "Mid"
time[grep("Late",names(dat1))] <- "Late"
time
colData <- data.frame(time)
summary(colData)

# check order
names(countData)
colData

# perform DESeq2
deds <- DESeqDataSetFromMatrix(countData=countData, 
                               colData=colData, 
                               design=~time)
dds <- DESeq(deds)

# results by contract
resEvM <- results(dds, contrast = c("time","Early","Mid"))
head(resEvM)
table(resEvM$padj<0.1)
# TRUE = 38

resEvL <- results(dds, contrast = c("time","Early","Late"))
head(resEvL)
table(resEvL$padj<0.05)
# TRUE = 2

resMvL <- results(dds, contrast = c("time","Mid","Late"))
head(resMvL)
table(resMvL$padj<0.05)
# TRUE = 31

plotMA(resEvM, main="Early vs. Mid")  
plotMA(resEvL, main="Early vs. Late")
plotMA(resMvL, main="Mid vs. Late")

write.table(resEvM, file="EarlyvsMid_DE.txt", quote=F, sep="\t")
write.table(resEvL, file="EarlyvsLate_DE.txt", quote=F, sep="\t")
write.table(resMvL, file="MidvsLate_DE.txt", quote=F, sep="\t")

# rlog transformation for plotting
rlog <- rlogTransformation(dds, blind=TRUE) 
rld <- assay(rlog)
head(rld)
colnames(rld) <- paste(colnames(countData))
head(rld)
length(rld[,1]) #496

valEvM <- cbind(resEvM$pvalue, resEvM$padj)
colnames(valEvM)=c("pval", "padj")
length(valEvM[,1])
table(complete.cases(valEvM))
rldpvalsEvM <- as.data.frame(cbind(rld,valEvM))
head(rldpvalsEvM)
write.csv(rldpvalsEvM, "RLDandPVALS_EvM.csv", quote=F)

valEvL <- cbind(resEvL$pvalue, resEvL$padj)
colnames(valEvL)=c("pval", "padj")
length(valEvL[,1])
table(complete.cases(valEvL))
rldpvalsEvL <- as.data.frame(cbind(rld,valEvL))
head(rldpvalsEvL)
write.csv(rldpvalsEvL, "RLDandPVALS_EvL.csv", quote=F)

valMvL <- cbind(resMvL$pvalue, resMvL$padj)
colnames(valMvL)=c("pval", "padj")
length(valMvL[,1])
table(complete.cases(valMvL))
rldpvalsMvL <- as.data.frame(cbind(rld,valMvL))
head(rldpvalsMvL)
write.csv(rldpvalsMvL, "RLDandPVALS_MvL.csv", quote=F)

# Heatmaps -------

# Mid vs. Late
rldpvals <- read.csv(file="RLDandPVALS_MvL.csv", row.names=1)
head(rldpvals)
p.val <- 0.05 # FDR cutoff
conds <- rldpvals[rldpvals$padj<=p.val & !is.na(rldpvals$padj),]
head(conds)
length(conds[,1])
# 28

rld <- conds %>% select(c(grep("Mid",names(conds)),grep("Late",names(conds))))
head(rld)
means=apply(rld,1,mean) # means of rows
head(means)
explc=rld-means # subtracting them
head(explc)

ann <- data.frame(cond = c(rep("Mid",3),rep("Late",6)))
rownames(ann) <- names(explc)
Var1        <- c("red", "green")
names(Var1) <-c("Mid", "Late")
anno_colors <- list(cond = Var1)

quartz()
pheatmap(as.matrix(explc),annotation_col=ann,annotation_colors=anno_colors,cex=1,
         color=myColor, breaks=myBreaks,border_color=NA,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation", show_rownames=T)

# Early vs. Mid
rldpvals <- read.csv(file="RLDandPVALS_EvM.csv", row.names=1)
head(rldpvals)
p.val <- 0.05 # FDR cutoff
conds <- rldpvals[rldpvals$padj<=p.val & !is.na(rldpvals$padj),]
head(conds)
length(conds[,1])
# 19

rld <- conds %>% select(c(grep("Early",names(conds)),grep("Mid",names(conds))))
head(rld)
means=apply(rld,1,mean) # means of rows
head(means)

explc=rld-means # subtracting them
head(explc)

ann <- data.frame(cond = c(rep("Early",6),rep("Mid",3)))
rownames(ann) <- names(explc)
Var1        <- c("orange", "red")
names(Var1) <-c("Early", "Mid")
anno_colors <- list(cond = Var1)

quartz()
pheatmap(as.matrix(explc),annotation_col=ann,annotation_colors=anno_colors,cex=1,
         color=myColor, breaks=myBreaks,border_color=NA,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation", show_rownames=T)

# Make big table for supplementary data
write.table(countData, file="ASV_taxa_counts.txt",sep="\t",quote=F)
write.table(res, file="DESeq2_results_inshore_v_offshore.txt", sep="\t", quote=F)
write.table(resEvL, file="DESeq2_results_early_v_late.txt", sep="\t", quote=F)
write.table(resEvM, file="DESeq2_results_early_v_mid.txt", sep="\t", quote=F)
write.table(resMvL, file="DESeq2_results_mid_v_late.txt", sep="\t", quote=F)

