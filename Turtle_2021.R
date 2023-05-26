######TURTLE 2021 THESIS & PUBLICATION DATA
#Set up working directory
#Can use "Session" > "Set Working Directory" > "Choose Directory" > Select file you want to use
setwd("~/Grad School/Thesis/Data")

####Load Packages
library(decontam) #Contamination cleaning
library(phyloseq)
library(ggplot2)
library(qiime2R)
library(base)
library(microbiome)
library(pgirmess) #Kruskal-Wallis Test
library(iNEXT) #Hill Number Analysis
library(RVAideMemoire) #PermANOVA Test
library(vegan)
library(venn) #Venn Diagrams
library(lme4) #GLMM
library(emmeans) #GLMM pairwise test
library(report)


#Load data, before loading make sure .tsv does not have words/phrases above ASV ID row, remove # from ASV
data <- read.delim("~/Grad School/Thesis/Data/feature-table_2021.tsv", row.names=1)
#Flip columns to rows
t.dat <- as.data.frame(t(data))
#Set t.dat as dat file
dat <- t.dat

####Convert Raw Data to Relative Abundance
dat.ra<-decostand(dat, method = "total")
##Export relative abundance table 
write.csv(dat.ra, "Turtle_2021_RelativeAbundance_Raw.csv")
##Export relative abundance table transposed (needed for merging of taxonomy)
dat.rat <- as.data.frame(t(dat.ra))
write.csv(dat.rat, "Turtle_2021_RelativeAbundance_RawT.csv")

#Import taxonomy table using the taxonomy file from QIIME (edit in excel, separate KPODCFGS and remove confidence)
taxonomy = read.table(file= "Taxonomy.txt", header = TRUE, sep ="\t", row.names = 1)
#Import metadata
metadata = read.delim("~/Grad School/Thesis/Data/Contamination_Metadata.txt", row.names=1)
#Load in QIIME tree
phy_tree = qza_to_phyloseq(tree="rooted-tree.qza")
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.UF = otu_table(as.matrix(dat), taxa_are_rows=FALSE)
tax.UF = tax_table(as.matrix(taxonomy))
meta.UF = sample_data(metadata)
#Check ASV name consistency
taxa_names(tax.UF)
taxa_names(ASV.UF)
taxa_names(phy_tree)
#Merge Files
physeq = phyloseq(ASV.UF,tax.UF,meta.UF,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Check raw counts
physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16646 taxa and 245 samples ]
# sample_data() Sample Data:       [ 245 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 16646 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 16646 tips and 16617 internal nodes ]


#########################DEEP CLEAN
#Identify Contaminants - Prevalence
#Threshold 0.5
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control"
contamdf.prev05 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
# FALSE  TRUE 
# 16533   113 
###Remove contaminants
which(contamdf.prev05$contaminant)
noncontam.prev <- prune_taxa(!contamdf.prev05$contaminant, physeq)
noncontam.prev
# otu_table()   OTU Table:         [ 16533 taxa and 245 samples ]
# sample_data() Sample Data:       [ 245 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 16533 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 16533 tips and 16506 internal nodes ]

#Export, Remove Controls
write.table(t(otu_table(noncontam.prev)), "NonContaminated_Prev.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


#Identify Contaminants - Frequency, using Prevelance Table
###Threshold=0.5
contamdf.freq05 <- isContaminant(noncontam.prev, method="frequency", conc="QIIME_Reads", threshold=0.5)
table(contamdf.freq05$contaminant)
# FALSE  TRUE 
# 13447  3086 
which(contamdf.freq05$contaminant)
###Remove contaminants
noncontam.freq <- prune_taxa(!contamdf.freq05$contaminant, physeq_clean)
noncontam.freq
# otu_table()   OTU Table:         [ 13447 taxa and 243 samples ]
# sample_data() Sample Data:       [ 243 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 13447 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 13447 tips and 13435 internal nodes ]

#Export and Use as Data Frame
write.table(t(otu_table(noncontam.freq)), "NonContaminated_Freq.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


#Import Data
data <- read.delim("~/Grad School/Thesis/Data/NonContaminated_Freq.txt", row.names=1)
#Flip columns to rows
t.dat <- as.data.frame(t(data))
#Set t.dat as dat
dat <- t.dat
###Remove ASVs that occur <0.001 ---> increases the number of ASVs - includes more "microdiversity" 
dat.pa<-decostand(dat, method ="pa") 
dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
dat.001per<-dat[,dat.otus.001per]

###Remove ASVs that occur <0.005
dat.otus.005per<-which(colSums(dat.pa) > (0.005*nrow(dat.pa)))
dat.005per<-dat[,dat.otus.005per]
#3,109 taxa

####Convert Clean Data to Relative Abundance
dat.ra<-decostand(dat.001per, method = "total")
#Export relative abundance table 
write.csv(dat.ra, "Turtle_2021_RelativeAbundance_05.csv")
#Export relative abundance table transposed (needed for merging of taxonomy)
dat.rat <- as.data.frame(t(dat.ra))
write.csv(dat.rat, "Turtle_2021_RelativeAbundance_05T.csv")


######################CONSERVATIVE CLEANING
#Identify Contaminants - Prevalence
###Threshold 0.1
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold = 0.1)
table(contamdf.prev$contaminant)
# FALSE  TRUE 
# 16625    21
###Remove contaminants
noncontam.prev <- prune_taxa(!contamdf.prev$contaminant, physeq)
noncontam.prev
# otu_table()   OTU Table:         [ 13447 taxa and 243 samples ]
# sample_data() Sample Data:       [ 243 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 13447 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 13447 tips and 13435 internal nodes ]

#Export, Remove Controls
write.table(t(otu_table(noncontam.prev)), "NonContaminated_Prev1.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

#Identify Contaminants - Frequency, using Prevelance Table
###Threshold=0.1
contamdf.freq01 <- isContaminant(noncontam.prev, method="frequency", conc="QIIME_Reads", threshold=0.1)
table(contamdf.freq01$contaminant)
# FALSE  TRUE 
# 16518  107 
###Remove contaminants
noncontam.freq01 <- prune_taxa(!contamdf.freq01$contaminant, physeq_clean)
noncontam.freq01
# otu_table()   OTU Table:         [ 16518 taxa and 243 samples ]
# sample_data() Sample Data:       [ 243 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 16518 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 16518 tips and 16489 internal nodes ]

#Export and Use as Data Frame
write.table(t(otu_table(noncontam.freq01)), "NonContaminated_Freq01.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

#Import Data
Data <- read.delim("~/Grad School/Thesis/Data/NonContaminated_Freq01.txt", row.names=1)
#Flip columns to rows
t.dat <- as.data.frame(t(Data))
#Set t.dat as dat
dat <- t.dat

###Remove ASVs that occur <0.001 ---> increases the number of ASVs - includes more "microdiversity" 
dat.pa<-decostand(dat, method ="pa") 
dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
dat.001per<-dat[,dat.otus.001per]
#16,516 taxa

#Export Clean Count Table
write.csv(dat.001per, "Turtle_2021_Count_01.csv")
#Export count table transposed (needed for merging of taxonomy)
dat.001per.t <- as.data.frame(t(dat.001per))
write.csv(dat.001per.t, "Turtle_2021_Count_01T.csv")

####Convert Clean Data to Relative Abundance
dat.ra<-decostand(dat.001per, method = "total")
#Export relative abundance table 
write.csv(dat.ra, "Turtle_2021_RelativeAbundance_01.csv")
#Export relative abundance table transposed (needed for merging of taxonomy)
dat.rat <- as.data.frame(t(dat.ra))
write.csv(dat.rat, "Turtle_2021_RelativeAbundance_01T.csv")
write.table(dat.rat, "Turtle_2021_RelativeAbundance_01T.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)





#################Merge Tables for Primer
#Import Data
dat <- read.delim("~/Grad School/Thesis/Data/Turtle_2021_Count_01T.txt", row.names=1)
#Import Taxonomy
taxonomy = read.table(file= "Taxonomy.txt", header = TRUE, sep ="\t", row.names = 1)

#Match Taxonomy in data and taxonomy files
common.rownames <- intersect(rownames(dat), rownames(taxonomy))
#Remove Taxonomy that got removed via Cleaning process 
dat <- dat[common.rownames,]
Taxonomy <- taxonomy[common.rownames,]
##Check that all rows match
all.equal(rownames(dat),rownames(Taxonomy))
#Merge Taxonomy and Data
Primer <- merge(dat, Taxonomy, by = 0)  
#Export as csv
write.csv(Primer, "Primer.csv")
####Merge Metadata in New Primer File using Excel (transpose metadata file and sort both files by column name)
#Saved as "Primer_Turtle2021_01Clean.xlsx"



###################ANALYSIS

###Look at Raw Data Stats
#Import Data: before loading make sure .tsv does not have words/phrases above ASV ID row, remove # from ASV
dat_raw <- read.delim("~/Grad School/Thesis/Data/feature-table_2021.tsv", row.names=1)
#Flip columns to rows
t.dat_raw <- as.data.frame(t(dat_raw))
#Set t.dat as dat
dat_raw <- t.dat_raw
#Import Taxonomy Table using the taxonomy file from QIIME (edit in excel, separate KPODCFGS and remove confidence)
taxonomy = read.table(file= "Taxonomy.txt", header = TRUE, sep ="\t", row.names = 1)
#Import metadata for cleaning
metadata = read.delim("~/Grad School/Thesis/Data/Contamination_Metadata.txt", row.names=1)
#Load in QIIME tree
phy_tree = qza_to_phyloseq(tree="rooted-tree.qza")
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.UF = otu_table(as.matrix(dat_raw), taxa_are_rows=FALSE)
tax.UF = tax_table(as.matrix(taxonomy))
meta.UF = sample_data(metadata)
#Check ASV name consistency
taxa_names(tax.UF)
taxa_names(ASV.UF)
taxa_names(phy_tree)
#Merge Files
physeq = phyloseq(ASV.UF,tax.UF,meta.UF,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Phyloseq
#######READ STATISTICS
#Check number of microbes observed in each sample
sample_sums(physeq)
#Basic stats for reads of samples
mean(sample_sums(physeq))
##Mean = 94,013.58
min(sample_sums(physeq))
##Min= 36,264
max(sample_sums(physeq))
##Max = 206,332
sd(sample_sums(physeq))
##SD = 27,693.09
#View ASV count
physeq
##16646 taxa and 245 samples


###Look at Conservative Clean Data Stats
#Import Data
dat <- read.delim("~/Grad School/Thesis/Data/NonContaminated_Freq01.txt", row.names=1)
#Flip columns to rows
t.dat <- as.data.frame(t(dat))
#Set t.dat as dat
dat <- t.dat
###Remove ASVs that occur <0.001 ---> increases the number of ASVs - includes more "microdiversity" 
dat.pa<-decostand(dat, method ="pa") 
dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
dat.001per<-dat[,dat.otus.001per]
#16,516 taxa

#Convert to relative abundance and standardize by sample total
dat.ra<-decostand(dat.001per, method = "total")

#Import Taxonomy Table using the taxonomy file from QIIME (edit in excel, separate KPODCFGS and remove confidence)
taxonomy = read.table(file= "Taxonomy.txt", header = TRUE, sep ="\t", row.names = 1)
#Import metadata for cleaning
metadata = read.delim("~/Grad School/Thesis/Data/Turtle_Metadata_2021.txt", row.names=1)
#Load in QIIME tree
phy_tree = qza_to_phyloseq(tree="rooted-tree.qza")
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.UF = otu_table(as.matrix(dat.001per), taxa_are_rows=FALSE)
tax.UF = tax_table(as.matrix(taxonomy))
meta.UF = sample_data(metadata)
#Check ASV name consistency
taxa_names(tax.UF)
taxa_names(ASV.UF)
taxa_names(phy_tree)
#Merge Files
physeq_clean = phyloseq(ASV.UF,tax.UF,meta.UF,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Phyloseq
#######READ STATISTICS
#Check number of microbes observed in each sample
sample_sums(physeq_clean)
#Basic stats for reads of samples
mean(sample_sums(physeq_clean))
##Mean = 92982.84
min(sample_sums(physeq_clean))
##Min= 35598
max(sample_sums(physeq_clean))
##Max = 206332
sd(sample_sums(physeq_clean))
##SD = 27573.97
#View ASV count
physeq_clean
##16516 taxa and 243 samples


##########Alpha Diversity - vegan package
###Diversity by Sample
#Species Richness:
S <- specnumber(dat.001per)
#No. individuals:
N <- rowSums(dat.001per)
#Margalef's Species Richness:
d = (S-1)/log(N)
#Shannon-Weiner Diversity (ln):
H <- diversity(dat.001per, index="shannon")
#Pielou's Evenness:
J <- H/log(N)
#Simpson's Diversity (1-D):
D <- diversity(dat.001per, index="simpson")
#Simpson's Diversity (1/D):
inv.D <- diversity(dat.001per, index="inv")

#Combine data together into a single new dataframe, export as CSV
diversity.by.sample <- cbind(S, N, d, H, J, D, inv.D)
write.csv(diversity.by.sample, "diversity.by.sample.csv")

###Diversity by Sample Type
#Merge Count and Environmental Data (only merge with factor wanted)
MergedData <- merge(dat, Turtle_Metadata_2021[1], by=0, all=T)
#Aggregate Count by Factor (Sample Type)
SumSampleType <- aggregate(. ~ Sample.Type, MergedData, sum)
#Make Sample Type row labels
SumSampleType.Label <- SumSampleType[,-1]
rownames(SumSampleType.Label) <- SumSampleType[,1]

#Species Richness
(S <- specnumber(SumSampleType.Label))
#No. individuals
(N <- rowSums(SumSampleType.Label))
#Margalef's Species Richness
(d = (S-1)/log(N))
#Shannon-Weiner Diversity (ln)
(H <- diversity(SumSampleType.Label, index="shannon"))
#Pielou's Evenness
(J <- H/log(N))
#Simpson's Diversity (1-D)
(D <- diversity(SumSampleType.Label, index="simpson"))
#Simpson's Diversity (1/D)
(inv.D <- diversity(SumSampleType.Label, index="invsimpson"))

#Bind all diversity indices into a new data frame
diversity.by.sampletype <- cbind(S, N, d, H, J, D, inv.D)
#Save the file
write.csv(diversity.by.sampletype, "diversity.by.sampletype.csv")


###Diversity by Species 
########CC
#Merge Count and Environmental Data (only merge with factor wanted)
MergedData <- merge(dat.001per, metadata[3], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Merge Sample Types
MergedData <- merge(MergedData.Label, metadata[1], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Remove CM species
datf <- droplevels(MergedData[!MergedData$Species == 'CM',])
datf.label <- datf[,-1]
rownames(datf.label) <- datf[,1]
CCdata <- subset(datf.label, select = -Species)
#Aggregate Count by Factor (Sample Type)
CCSumSampleType <- aggregate(. ~ Sample.Type, CCdata, sum)
#Make Sample Type row labels
CCSumSampleType.Label <- CCSumSampleType[,-1]
rownames(CCSumSampleType.Label) <- CCSumSampleType[,1]

#Species Richness
(S <- specnumber(CCSumSampleType.Label))
#No. individuals
(N <- rowSums(CCSumSampleType.Label))
#Margalef's Species Richness
(d = (S-1)/log(N))
#Shannon-Weiner Diversity (ln)
(H <- diversity(CCSumSampleType.Label, index="shannon"))
#Pielou's Evenness
(J <- H/log(N))
#Simpson's Diversity (1-D)
(D <- diversity(CCSumSampleType.Label, index="simpson"))
#Simpson's Diversity (1/D)
(inv.D <- diversity(CCSumSampleType.Label, index="invsimpson"))

#Bind all diversity indices into a new data frame
diversity.by.CC <- cbind(S, N, d, H, J, D, inv.D)
#Save the file
write.csv(diversity.by.CC, "diversity.by.CC.csv")

#######CM
#Merge Count and Environmental Data (only merge with factor wanted)
MergedData <- merge(dat, metadata[3], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Merge Sample Types
MergedData <- merge(MergedData.Label, metadata[1], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Remove CM species
datf <- droplevels(MergedData[!MergedData$Species == 'CC',])
datf.label <- datf[,-1]
rownames(datf.label) <- datf[,1]
CMdata <- subset(datf.label, select = -Species)
#Aggregate Count by Factor (Sample Type)
CMSumSampleType <- aggregate(. ~ Sample.Type, CMdata, sum)
#Make Sample Type row labels
CMSumSampleType.Label <- CMSumSampleType[,-1]
rownames(CMSumSampleType.Label) <- CMSumSampleType[,1]

#Species Richness
(S <- specnumber(CMSumSampleType.Label))
#No. individuals
(N <- rowSums(CMSumSampleType.Label))
#Margalef's Species Richness
(d = (S-1)/log(N))
#Shannon-Weiner Diversity (ln)
(H <- diversity(CMSumSampleType.Label, index="shannon"))
#Pielou's Evenness
(J <- H/log(N))
#Simpson's Diversity (1-D)
(D <- diversity(CMSumSampleType.Label, index="simpson"))
#Simpson's Diversity (1/D)
(inv.D <- diversity(CMSumSampleType.Label, index="invsimpson"))

#Bind all diversity indices into a new data frame
diversity.by.CM <- cbind(S, N, d, H, J, D, inv.D)
#Save the file
write.csv(diversity.by.CM, "diversity.by.CM.csv")


###Diversity by Beach (CC only)
########FT
#Merge Count and Environmental Data (only merge with factor wanted)
MergedData <- merge(dat, metadata[3], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Merge Sample Types
MergedData <- merge(MergedData.Label, metadata[1], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Merge Beach
MergedData <- merge(MergedData.Label, metadata[2], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Remove CM species
datf <- droplevels(MergedData[!MergedData$Species == 'CM',])
datf.label <- datf[,-1]
rownames(datf.label) <- datf[,1]
CCdata <- subset(datf.label, select = -Species)
#Remove Hillsboro beach
datf <- droplevels(CCdata[!CCdata$Beach == 'Hillsboro',])
FTdata <- subset(datf, select = -Beach)

#Aggregate Count by Factor (Sample Type)
FTSumSampleType <- aggregate(. ~ Sample.Type, FTdata, sum)
#Make Sample Type row labels
FTSumSampleType.Label <- FTSumSampleType[,-1]
rownames(FTSumSampleType.Label) <- FTSumSampleType[,1]

#Species Richness
(S <- specnumber(FTSumSampleType.Label))
#No. individuals
(N <- rowSums(FTSumSampleType.Label))
#Margalef's Species Richness
(d = (S-1)/log(N))
#Shannon-Weiner Diversity (ln)
(H <- diversity(FTSumSampleType.Label, index="shannon"))
#Pielou's Evenness
(J <- H/log(N))
#Simpson's Diversity (1-D)
(D <- diversity(FTSumSampleType.Label, index="simpson"))
#Simpson's Diversity (1/D)
(inv.D <- diversity(FTSumSampleType.Label, index="invsimpson"))

#Bind all diversity indices into a new data frame
diversity.by.FT <- cbind(S, N, d, H, J, D, inv.D)
#Save the file
write.csv(diversity.by.FT, "diversity.by.FT.csv")

########Hillsboro
#Merge Count and Environmental Data (only merge with factor wanted)
MergedData <- merge(dat, metadata[3], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Merge Sample Types
MergedData <- merge(MergedData.Label, metadata[1], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Merge Beach
MergedData <- merge(MergedData.Label, metadata[2], by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Remove CM species
datf <- droplevels(MergedData[!MergedData$Species == 'CM',])
datf.label <- datf[,-1]
rownames(datf.label) <- datf[,1]
CCdata <- subset(datf.label, select = -Species)
#Remove Hillsboro beach
datf <- droplevels(CCdata[!CCdata$Beach == 'Fort Lauderdale',])
Hdata <- subset(datf, select = -Beach)

#Aggregate Count by Factor (Sample Type)
HSumSampleType <- aggregate(. ~ Sample.Type, Hdata, sum)
#Make Sample Type row labels
HSumSampleType.Label <- HSumSampleType[,-1]
rownames(HSumSampleType.Label) <- HSumSampleType[,1]

#Species Richness
(S <- specnumber(HSumSampleType.Label))
#No. individuals
(N <- rowSums(HSumSampleType.Label))
#Margalef's Species Richness
(d = (S-1)/log(N))
#Shannon-Weiner Diversity (ln)
(H <- diversity(HSumSampleType.Label, index="shannon"))
#Pielou's Evenness
(J <- H/log(N))
#Simpson's Diversity (1-D)
(D <- diversity(HSumSampleType.Label, index="simpson"))
#Simpson's Diversity (1/D)
(inv.D <- diversity(HSumSampleType.Label, index="invsimpson"))

#Bind all diversity indices into a new data frame
diversity.by.H <- cbind(S, N, d, H, J, D, inv.D)
#Save the file
write.csv(diversity.by.H, "diversity.by.H.csv")

###PLOT Diversity Tables
AlphaDiversity_ALL <- read.delim("~/Grad School/Thesis/Data/AlphaDiversity_ALL.txt", row.names=1)
AlphaDiversity_CC <- read.delim("~/Grad School/Thesis/Data/AlphaDiversity_CC.txt", row.names=1)
AlphaDiversity_CM <- read.delim("~/Grad School/Thesis/Data/AlphaDiversity_CM.txt", row.names=1)
AlphaDiversity_F <- read.delim("~/Grad School/Thesis/Data/AlphaDiversity_F.txt", row.names=1)
AlphaDiversity_H <- read.delim("~/Grad School/Thesis/Data/AlphaDiversity_H.txt", row.names=1)
AlphaDiversity_Species <- read.delim("~/Grad School/Thesis/Data/AlphaDiversity_Species.txt", row.names=1)
AlphaDiversity_Beaches <- read.delim("~/Grad School/Thesis/Data/AlphaDiversity_Beaches.txt", row.names=1)

#ALL SAMPLES
group.colors <- c("Cloaca"="blue", "Control Sand" = "pink", "Hatched Egg" = "red", "Nest Sand" = "green", "Unhatched Egg" = "turquoise")
ggplot(AlphaDiversity_ALL) +
  geom_boxplot(aes(x = Group, y = S, fill = Group)) + ylab("Total Species (S)") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_ALL) +
  geom_boxplot(aes(x = Group, y = N, fill = Group)) + ylab("Total Individuals (N)") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_ALL) +
  geom_boxplot(aes(x = Group, y = d, fill = Group)) + ylab("Margalef's Species Richness (d)")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_ALL) +
  geom_boxplot(aes(x = Group, y = J, fill = Group)) + ylab("Pielou's Evenness (J')")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_ALL) +
  geom_boxplot(aes(x = Group, y = D, fill = Group)) + ylab("Inverse Simpson's Diversity (1- ??)")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_ALL) +
  geom_boxplot(aes(x = Group, y = H, fill = Group)) + ylab("Shannon Diversity (H')")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=group.colors)

#CC SAMPLES
group.colors <- c("Cloaca"="blue", "Control Sand" = "pink", "Hatched Egg" = "red", "Nest Sand" = "green", "Unhatched Egg" = "turquoise")
windows(15,10)
ggplot(AlphaDiversity_CC, aes(x = Group, y = S, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors) +
ggplot(AlphaDiversity_CC, aes(x = Group, y = N, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CC, aes(x = Group, y = d, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CC, aes(x = Group, y = J, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CC, aes(x = Group, y = D, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CC, aes(x = Group, y = H, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)

#CM SAMPLES
group.colors <- c("Cloaca"="blue", "Control Sand" = "pink", "Hatched Egg" = "red", "Nest Sand" = "green", "Unhatched Egg" = "turquoise")
windows(15,10)
ggplot(AlphaDiversity_CM, aes(x = Group, y = S, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CM, aes(x = Group, y = N, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CM, aes(x = Group, y = d, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CM, aes(x = Group, y = J, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CM, aes(x = Group, y = D, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_CM, aes(x = Group, y = H, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)

#FT SAMPLES
group.colors <- c("Cloaca"="blue", "Control Sand" = "pink", "Hatched Egg" = "red", "Nest Sand" = "green", "Unhatched Egg" = "turquoise")
windows(15,10)
ggplot(AlphaDiversity_F, aes(x = Group, y = S, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_F, aes(x = Group, y = N, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_F, aes(x = Group, y = d, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_F, aes(x = Group, y = J, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_F, aes(x = Group, y = D, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_F, aes(x = Group, y = H, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)

#Hillsboro SAMPLES
group.colors <- c("Cloaca"="blue", "Control Sand" = "pink", "Hatched Egg" = "red", "Nest Sand" = "green", "Unhatched Egg" = "turquoise")
windows(15,10)
ggplot(AlphaDiversity_H, aes(x = Group, y = S, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_H, aes(x = Group, y = N, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_H, aes(x = Group, y = d, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_H, aes(x = Group, y = J, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_H, aes(x = Group, y = D, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)
ggplot(AlphaDiversity_H, aes(x = Group, y = H, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.colors)


#Both Species Alpha Combined
species.colors <- c("CC"="dark orange", "CM" = "forestgreen")
ggplot(AlphaDiversity_Species) +
  geom_boxplot(aes(x = Group, y = S, fill = Species)) + ylab("Total Species (S)") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=species.colors)
ggplot(AlphaDiversity_Species) +
  geom_boxplot(aes(x = Group, y = N, fill = Species)) + ylab("Total Individuals (N)") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=species.colors)
ggplot(AlphaDiversity_Species) +
  geom_boxplot(aes(x = Group, y = d, fill = Species)) + ylab("Margalef's Species Richness (d)")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=species.colors)
ggplot(AlphaDiversity_Species) +
  geom_boxplot(aes(x = Group, y = J, fill = Species)) + ylab("Pielou's Evenness (J')")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=species.colors)
ggplot(AlphaDiversity_Species) +
  geom_boxplot(aes(x = Group, y = D, fill = Species)) + ylab("Inverse Simpson's Diversity (1- ??)")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=species.colors)
ggplot(AlphaDiversity_Species) +
  geom_boxplot(aes(x = Group, y = H, fill = Species)) + ylab("Shannon Diversity (H')")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=species.colors)

#Both Beaches Alpha Combined
beaches.colors <- c("FT"="red", "H" = "blue")
ggplot(AlphaDiversity_Beaches) +
  geom_boxplot(aes(x = Group, y = S, fill = Beach)) + ylab("Total Species (S)") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=beaches.colors)
ggplot(AlphaDiversity_Beaches) +
  geom_boxplot(aes(x = Group, y = N, fill = Beach)) + ylab("Total Individuals (N)") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=beaches.colors)
ggplot(AlphaDiversity_Beaches) +
  geom_boxplot(aes(x = Group, y = d, fill = Beach)) + ylab("Margalef's Species Richness (d)")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=beaches.colors)
ggplot(AlphaDiversity_Beaches) +
  geom_boxplot(aes(x = Group, y = J, fill = Beach)) + ylab("Pielou's Evenness (J')")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=beaches.colors)
ggplot(AlphaDiversity_Beaches) +
  geom_boxplot(aes(x = Group, y = D, fill = Beach)) + ylab("Inverse Simpson's Diversity (1- ??)")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=beaches.colors)
ggplot(AlphaDiversity_Beaches) +
  geom_boxplot(aes(x = Group, y = H, fill = Beach)) + ylab("Shannon Diversity (H')")  + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_manual(values=beaches.colors)


#############Final Alpha Diversity Plot (Figure 1)
###Create file with alpha diversity separated by sample and include two columns identifying beach and turtle species
AlphaDiversity_AllSplit <- read.delim("~/Grad School/Thesis/Data/AlphaDiversity_AllSplit.txt", row.names=1)

#ALL SAMPLES
ggplot(AlphaDiversity_AllSplit) +
  geom_boxplot(aes(x = Group, y = S, fill = Group)) + ylab("Total Species (S)") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_grey(start=0.4, end=1.0) +
  facet_grid(. ~ Sample.Type) +
  theme(strip.text.x = element_text(size = 17))
ggplot(AlphaDiversity_AllSplit) +
  geom_boxplot(aes(x = Group, y = d, fill = Group)) + ylab("Margalef's Species Richness (d)") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_grey(start=0.4, end=1.0) +
  facet_grid(. ~ Sample.Type) +
  theme(strip.text.x = element_text(size = 17))
ggplot(AlphaDiversity_AllSplit) +
  geom_boxplot(aes(x = Group, y = D, fill = Group)) + ylab("Inverse Simpson's Diversity") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_grey(start=0.4, end=1.0) +
  facet_grid(. ~ Sample.Type) +
  theme(strip.text.x = element_text(size = 17))
ggplot(AlphaDiversity_AllSplit) +
  geom_boxplot(aes(x = Group, y = H, fill = Group)) + ylab("Shannon Diversity") + xlab("") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_fill_grey(start=0.4, end=1.0) +
  facet_grid(. ~ Sample.Type) +
  theme(strip.text.x = element_text(size = 17))



###ALL SAMPLES Statistical Significance
#Shapiro Test = Normality
shapiro.test(AlphaDiversity_ALL$S)
#W = 0.61204, p-value < 2.2e-16
shapiro.test(AlphaDiversity_ALL$N)
#W = 0.94261, p-value = 3.63e-08
shapiro.test(AlphaDiversity_ALL$d)
#W = 0.60183, p-value < 2.2e-16
shapiro.test(AlphaDiversity_ALL$H)
#W = 0.89851, p-value = 9.789e-12
shapiro.test(AlphaDiversity_ALL$J)
#W = 0.88216, p-value = 8.548e-13
shapiro.test(AlphaDiversity_ALL$D)
#W = 0.55447, p-value < 2.2e-16

##p<0.05 = not normally distributed >> transformations did not work = non-parametric (KRUSKAL-WALLIS)

########Kruskal Wallis: Nonparametric Data (not normal)
kruskal.test(AlphaDiversity_ALL$S ~ AlphaDiversity_ALL$Group)
#Kruskal-Wallis chi-squared = 144.44, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_ALL$S, AlphaDiversity_ALL$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_ALL$N ~ AlphaDiversity_ALL$Group)
#Kruskal-Wallis chi-squared = 24.508, df = 4, p-value = 6.318e-05
pairwise.wilcox.test(AlphaDiversity_ALL$N, AlphaDiversity_ALL$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_ALL$d ~ AlphaDiversity_ALL$Group)
#Kruskal-Wallis chi-squared = 146.76, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_ALL$d, AlphaDiversity_ALL$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_ALL$H ~ AlphaDiversity_ALL$Group)
#Kruskal-Wallis chi-squared = 106.35, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_ALL$H, AlphaDiversity_ALL$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_ALL$J ~ AlphaDiversity_ALL$Group)
#Kruskal-Wallis chi-squared = 104.35, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_ALL$J, AlphaDiversity_ALL$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_ALL$D ~ AlphaDiversity_ALL$Group)
#Kruskal-Wallis chi-squared = 89.151, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_ALL$D, AlphaDiversity_ALL$Group, p.adjust.method = "fdr")

###CC Statistical Significance
#Shapiro Test = Normality
shapiro.test(AlphaDiversity_CC$S)
#W = 0.63831, p-value < 2.2e-16
shapiro.test(AlphaDiversity_CC$N)
#W = 0.94592, p-value = 2.451e-06
shapiro.test(AlphaDiversity_CC$d)
#W = 0.62809, p-value < 2.2e-16
shapiro.test(AlphaDiversity_CC$H)
#W = 0.92505, p-value = 5.311e-08
shapiro.test(AlphaDiversity_CC$J)
#W = 0.90945, p-value = 4.505e-09
shapiro.test(AlphaDiversity_CC$D)
#W = 0.61523, p-value < 2.2e-16

##p<0.05 = not normally distributed >> transformations did not work = non-parametric (KRUSKAL-WALLIS)

########Kruskal Wallis: Nonparametric Data (not normal)
kruskal.test(AlphaDiversity_CC$S ~ AlphaDiversity_CC$Group)
#Kruskal-Wallis chi-squared = 114.27, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_CC$S, AlphaDiversity_CC$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CC$N ~ AlphaDiversity_CC$Group)
#Kruskal-Wallis chi-squared = 22.881, df = 4, p-value = 0.0001337
pairwise.wilcox.test(AlphaDiversity_CC$N, AlphaDiversity_CC$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CC$d ~ AlphaDiversity_CC$Group)
#Kruskal-Wallis chi-squared = 117.09, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_CC$d, AlphaDiversity_CC$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CC$H ~ AlphaDiversity_CC$Group)
#Kruskal-Wallis chi-squared = 97.851, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_CC$H, AlphaDiversity_CC$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CC$J ~ AlphaDiversity_CC$Group)
#Kruskal-Wallis chi-squared = 98.017, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_CC$J, AlphaDiversity_CC$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CC$D ~ AlphaDiversity_CC$Group)
#Kruskal-Wallis chi-squared = 85.974, df = 4, p-value < 2.2e-16
pairwise.wilcox.test(AlphaDiversity_CC$D, AlphaDiversity_CC$Group, p.adjust.method = "fdr")


###CM Statistical Significance
#Shapiro Test = Normality
shapiro.test(AlphaDiversity_CM$S)
#W = 0.678, p-value = 1.785e-10
shapiro.test(AlphaDiversity_CM$N)
#W = 0.97241, p-value = 0.1684
##NO significant differences
shapiro.test(AlphaDiversity_CM$d)
#W = 0.6783, p-value = 1.807e-10
shapiro.test(AlphaDiversity_CM$H)
#W = 0.79808, p-value = 7.083e-08
shapiro.test(AlphaDiversity_CM$J)
#W = 0.79232, p-value = 5.078e-08
shapiro.test(AlphaDiversity_CM$D)
#W = 0.29993, p-value = 9.286e-16

##p<0.05 = not normally distributed >> transformations did not work = non-parametric (KRUSKAL-WALLIS)
##p>0.05 = normally distributed (ANOVA)

#ANOVA: Normal distribution
anova.N<-aov(AlphaDiversity_CM$N ~ factor(AlphaDiversity_CM$Group))
summary(anova.N)
#p=0.655
TukeyHSD(anova.N)
##No significant pairwise comparisons

########Kruskal Wallis: Nonparametric Data (not normal)
kruskal.test(AlphaDiversity_CM$S ~ AlphaDiversity_CM$Group)
#Kruskal-Wallis chi-squared = 34.342, df = 4, p-value = 6.339e-07
pairwise.wilcox.test(AlphaDiversity_CM$S, AlphaDiversity_CM$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CM$d ~ AlphaDiversity_CM$Group)
#Kruskal-Wallis chi-squared = 34.604, df = 4, p-value = 5.601e-07
pairwise.wilcox.test(AlphaDiversity_CM$d, AlphaDiversity_CM$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CM$H ~ AlphaDiversity_CM$Group)
#Kruskal-Wallis chi-squared = 19.246, df = 4, p-value = 0.000703
pairwise.wilcox.test(AlphaDiversity_CM$H, AlphaDiversity_CM$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CM$J ~ AlphaDiversity_CM$Group)
#Kruskal-Wallis chi-squared = 17.35, df = 4, p-value = 0.001653
pairwise.wilcox.test(AlphaDiversity_CM$J, AlphaDiversity_CM$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_CM$D ~ AlphaDiversity_CM$Group)
#Kruskal-Wallis chi-squared = 23.946, df = 4, p-value = 8.188e-05
pairwise.wilcox.test(AlphaDiversity_CM$D, AlphaDiversity_CM$Group, p.adjust.method = "fdr")


###FT SAMPLES Statistical Significance
#Shapiro Test = Normality
shapiro.test(AlphaDiversity_F$S)
#W = 0.57148, p-value = 1.534e-15
shapiro.test(AlphaDiversity_F$N)
#W = 0.98176, p-value = 0.1871
shapiro.test(AlphaDiversity_F$d)
#W = 0.56001, p-value = 9.577e-16
shapiro.test(AlphaDiversity_F$H)
#W = 0.8971, p-value = 1.149e-06
shapiro.test(AlphaDiversity_F$J)
#W = 0.87249, p-value = 9.935e-08
shapiro.test(AlphaDiversity_F$D)
#W = 0.59405, p-value = 3.982e-15

##p<0.05 = not normally distributed >> transformations did not work = non-parametric (KRUSKAL-WALLIS)
##p>0.05 = normally distributed (ANOVA)

#ANOVA: Normal distribution
anova.N<-aov(AlphaDiversity_F$N ~ factor(AlphaDiversity_F$Group))
summary(anova.N)
#p=0.0377
TukeyHSD(anova.N)
##No significant pairwise comparisons


########Kruskal Wallis: Nonparametric Data (not normal)
kruskal.test(AlphaDiversity_F$S ~ AlphaDiversity_F$Group)
#Kruskal-Wallis chi-squared = 57.665, df = 4, p-value = 8.971e-12
pairwise.wilcox.test(AlphaDiversity_F$S, AlphaDiversity_F$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_F$d ~ AlphaDiversity_F$Group)
#Kruskal-Wallis chi-squared = 58.936, df = 4, p-value = 4.853e-12
pairwise.wilcox.test(AlphaDiversity_F$d, AlphaDiversity_F$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_F$H ~ AlphaDiversity_F$Group)
#Kruskal-Wallis chi-squared = 44.761, df = 4, p-value = 4.457e-09
pairwise.wilcox.test(AlphaDiversity_F$H, AlphaDiversity_F$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_F$J ~ AlphaDiversity_F$Group)
#Kruskal-Wallis chi-squared = 45.961, df = 4, p-value = 2.509e-09
pairwise.wilcox.test(AlphaDiversity_F$J, AlphaDiversity_F$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_F$D ~ AlphaDiversity_F$Group)
#Kruskal-Wallis chi-squared = 38.56, df = 4, p-value = 8.587e-08
pairwise.wilcox.test(AlphaDiversity_F$D, AlphaDiversity_F$Group, p.adjust.method = "fdr")


###H Statistical Significance
#Shapiro Test = Normality
shapiro.test(AlphaDiversity_H$S)
#W = 0.71233, p-value = 2.64e-11
shapiro.test(log(AlphaDiversity_H$N))
#W = 0.98748, p-value = 0.6218
shapiro.test(AlphaDiversity_H$d)
#W = 0.7031, p-value = 1.655e-11
shapiro.test(AlphaDiversity_H$H)
#W = 0.9473, p-value = 0.002246
shapiro.test(AlphaDiversity_H$J)
#W = 0.93835, p-value = 0.0007225
shapiro.test(AlphaDiversity_H$D)
#W = 0.63515, p-value = 7.007e-13

##p<0.05 = not normally distributed >> transformations did not work = non-parametric (KRUSKAL-WALLIS)
##p>0.05 = normally distributed (ANOVA)

#ANOVA: Normal distribution
anova.N<-aov(log(AlphaDiversity_H$N) ~ factor(AlphaDiversity_H$Group))
summary(anova.N)
#p=0.000852
TukeyHSD(anova.N)
##Significant pairwise comparisons

########Kruskal Wallis: Nonparametric Data (not normal)
kruskal.test(AlphaDiversity_H$S ~ AlphaDiversity_H$Group)
#Kruskal-Wallis chi-squared = 58.803, df = 4, p-value = 5.176e-12
pairwise.wilcox.test(AlphaDiversity_H$S, AlphaDiversity_H$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_H$d ~ AlphaDiversity_H$Group)
#Kruskal-Wallis chi-squared = 59.879, df = 4, p-value = 3.075e-12
pairwise.wilcox.test(AlphaDiversity_H$d, AlphaDiversity_H$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_H$H ~ AlphaDiversity_H$Group)
#Kruskal-Wallis chi-squared = 56.982, df = 4, p-value = 1.248e-11
pairwise.wilcox.test(AlphaDiversity_H$H, AlphaDiversity_H$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_H$J ~ AlphaDiversity_H$Group)
#Kruskal-Wallis chi-squared = 56.557, df = 4, p-value = 1.532e-11
pairwise.wilcox.test(AlphaDiversity_H$J, AlphaDiversity_H$Group, p.adjust.method = "fdr")

kruskal.test(AlphaDiversity_H$D ~ AlphaDiversity_H$Group)
#Kruskal-Wallis chi-squared = 51.219, df = 4, p-value = 2.009e-10
pairwise.wilcox.test(AlphaDiversity_H$D, AlphaDiversity_H$Group, p.adjust.method = "fdr")



#########GLMM Alpha Diversity
#Create a new object for common rownames for hatched and unhatched egg data
common.rownames <- intersect(rownames(AlphaDiversity_ALL),rownames(metadata))
#Set the data file and metadata file to have only the data that includes these common names 
AD <- AlphaDiversity_ALL[common.rownames,]
metadata <- metadata[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(AD), rownames(metadata), ignore.row.order = TRUE)
#Create one file with alpha diversity metric and metadata
glmmdata_ALL = merge(AD, metadata, by=0)

###Look for covariation among metadata
Turtle_Metadata_2021_simplified <- read.delim("~/Grad School/Thesis/Data/Turtle_Metadata_2021_simplified.txt", row.names=1)
metadata_covar = subset(Turtle_Metadata_2021_simplified, select = -c(Beach, Species, R.Zone, Recent.Renourishment, Turtle.ID, Nesting.Date, Excavation.Date, LIN, LPIP, DIN, DPIP, Whole, Hatched, Roots, Washover, Relocated, Predated, Temperature.Average, pH.Average, Conductivity.Average))
#Run correlation
cor(metadata_covar)
###Look for correlations above out of the -0.6 to 0.6 range (strong correlation) for removal
# pH.side - pH.Bottom = 0.86836426
# Temperature.Side - Temperature.Bottom = 0.87331573
# Latitude-Longitude = 0.99557176
# Dune.Distance-Latitude = -0.79102842
# Dune.Distance-Longitude = -0.75737137

############Margalef's Species Richness
###Create full model
mod.1 = lmer(d ~ Incubation.Length + Clutch.Size + Hatch.Success + 
                Chamber.Depth + High.Tide.Distance + Dune.Distance + 
                Temperature.Bottom + pH.Bottom + 
                Conductivity.Side + Conductivity.Bottom + 
                Sand.Grain.Size + Sorting.Coefficient + 
                Sample.Type*Species + Sample.Type*Beach +
                (1|Nest.Number), data = glmmdata_ALL)

summary(mod.1)
drop1(mod.1)
#                           AIC
# <none>                   1934.7
# Incubation.Length      1 1939.4
# Clutch.Size            1 1932.9
# Hatch.Success          1 1936.4
# Chamber.Depth          1 1936.8
# High.Tide.Distance     1 1938.7
# Dune.Distance          1 1935.1
# Temperature.Bottom     1 1942.9
# pH.Bottom              1 1932.8
# Conductivity.Side      1 1934.9
# Conductivity.Bottom    1 1936.1
# Sand.Grain.Size        1 1932.8
# Sorting.Coefficient    1 1933.4
# Sample.Type:Species    4 2032.1
# Sample.Type:Beach      4 1998.8

#Remove terms
#Keep terms with AIC scores that increase by 2 compared to the <none>
mod.2 = lmer(d ~ Incubation.Length + High.Tide.Distance + Temperature.Bottom +  
               Sample.Type*Species + Sample.Type*Beach +
               (1|Nest.Number), data = glmmdata_ALL)
summary(mod.2)
report(mod.2)

#Validate Model (mean should be roughly at zero)
plot(mod.2)
plot(residuals(mod.2)~fitted(mod.2))
plot(residuals(mod.2)~ glmmdata_ALL$Incubation.Length)
plot(residuals(mod.2)~ glmmdata_ALL$High.Tide.Distance)
plot(residuals(mod.2)~ glmmdata_ALL$Temperature.Bottom)
plot(residuals(mod.2)~factor(glmmdata_ALL$Beach))
plot(residuals(mod.2)~factor(glmmdata_ALL$Species))
plot(residuals(mod.2)~factor(glmmdata_ALL$Sample.Type))
plot(predict(mod.2, type = 'response')~ glmmdata_ALL$d)

#Pairwise
pw.tests_1 <- pairs(emmeans(mod.2, ~ Sample.Type, by=c("Species")))
pw.tests_1
pw.tests_2 <- pairs(emmeans(mod.2, ~ Sample.Type, by=c("Species","Beach")))
pw.tests_2
pw.tests_3 <- pairs(emmeans(mod.2, ~ Species, by=c("Sample.Type")))
pw.tests_3
pw.tests_4 <- pairs(emmeans(mod.2, ~ Beach, by=c("Sample.Type", "Species")))
pw.tests_4



############Shannon Diversity
###Create full model
mod.1 = lmer(H ~ Incubation.Length + Clutch.Size + Hatch.Success + 
               Chamber.Depth + High.Tide.Distance + Dune.Distance + 
               Temperature.Bottom + pH.Bottom + 
               Conductivity.Side + Conductivity.Bottom + 
               Sand.Grain.Size + Sorting.Coefficient + 
               Sample.Type*Species + Sample.Type*Beach +
               (1|Nest.Number), data = glmmdata_ALL)

summary(mod.1)
drop1(mod.1)
#                           AIC
# <none>                   480.56
# Incubation.Length      1 478.78
# Clutch.Size            1 489.44
# Hatch.Success          1 478.66
# Chamber.Depth          1 479.90
# High.Tide.Distance     1 483.15
# Dune.Distance          1 479.43
# Temperature.Bottom     1 485.16
# pH.Bottom              1 480.33
# Conductivity.Side      1 480.36
# Conductivity.Bottom    1 479.52
# Sand.Grain.Size        1 478.62
# Sorting.Coefficient    1 478.57
# Sample.Type:Species    4 556.93
# Sample.Type:Beach      4 508.37

#Remove terms
#Keep terms with AIC scores that increase by 2 compared to the <none>
mod.2 = lmer(H ~ Clutch.Size + High.Tide.Distance + Temperature.Bottom +
               Sample.Type*Species + Sample.Type*Beach +
               (1|Nest.Number), data = glmmdata_ALL)
summary(mod.2)
report(mod.2)

#Validate Model (mean should be roughly at zero)
plot(mod.2)
plot(residuals(mod.2)~fitted(mod.2))
plot(residuals(mod.2)~ glmmdata_ALL$Clutch.Size)
plot(residuals(mod.2)~ glmmdata_ALL$High.Tide.Distance)
plot(residuals(mod.2)~ glmmdata_ALL$Temperature.Bottom)
plot(residuals(mod.2)~factor(glmmdata_ALL$Beach))
plot(residuals(mod.2)~factor(glmmdata_ALL$Species))
plot(residuals(mod.2)~factor(glmmdata_ALL$Sample.Type))
plot(predict(mod.2, type = 'response')~ glmmdata_ALL$H)

#Pairwise
pw.tests_1 <- pairs(emmeans(mod.2, ~ Sample.Type, by=c("Species")))
pw.tests_1
pw.tests_2 <- pairs(emmeans(mod.2, ~ Sample.Type, by=c("Species","Beach")))
pw.tests_2
pw.tests_3 <- pairs(emmeans(mod.2, ~ Species, by=c("Sample.Type")))
pw.tests_3
pw.tests_4 <- pairs(emmeans(mod.2, ~ Beach, by=c("Sample.Type", "Species")))
pw.tests_4




############Inverse Simpson Diversity
###Create full model
mod.1 = lmer(D ~ Incubation.Length + Clutch.Size + Hatch.Success + 
               Chamber.Depth + High.Tide.Distance + Dune.Distance + 
               Temperature.Bottom + pH.Bottom + 
               Conductivity.Side + Conductivity.Bottom + 
               Sand.Grain.Size + Sorting.Coefficient + 
               Sample.Type*Species + Sample.Type*Beach +
               (1|Nest.Number), data = glmmdata_ALL)

summary(mod.1)
drop1(mod.1)
#                           AIC
# <none>                   -387.90
# Incubation.Length      1 -388.42
# Clutch.Size            1 -377.85
# Hatch.Success          1 -387.73
# Chamber.Depth          1 -389.89
# High.Tide.Distance     1 -388.95
# Dune.Distance          1 -389.71
# Temperature.Bottom     1 -389.13
# pH.Bottom              1 -388.47
# Conductivity.Side      1 -389.89
# Conductivity.Bottom    1 -388.14
# Sand.Grain.Size        1 -389.89
# Sorting.Coefficient    1 -389.90
# Sample.Type:Species    4 -381.31
# Sample.Type:Beach      4 -392.88

#Remove terms
#Keep terms with AIC scores that increase by 2 compared to the <none>
mod.2 = lmer(D ~ Sorting.Coefficient + Sample.Type*Beach + 
               (1|Nest.Number), data = glmmdata_ALL)
summary(mod.2)
report(mod.2)

#Validate Model (mean should be roughly at zero)
plot(mod.2)
plot(residuals(mod.2)~fitted(mod.2))
plot(residuals(mod.2)~ glmmdata_ALL$Sorting.Coefficient)
plot(residuals(mod.2)~factor(glmmdata_ALL$Beach))
plot(residuals(mod.2)~factor(glmmdata_ALL$Sample.Type))
plot(predict(mod.2, type = 'response')~ glmmdata_ALL$D)

#Pairwise
pw.tests_1 <- pairs(emmeans(mod.2, ~ Sample.Type, by=c("Beach")))
pw.tests_1
pw.tests_2 <- pairs(emmeans(mod.2, ~ Beach, by=c("Sample.Type")))
pw.tests_2







##########CREATE SEPARATED DATA FILES - by Turtle Species, Beach, Sample Types
metadata = read.delim("~/Grad School/Thesis/Data/Turtle_Metadata_2021.txt", row.names=1)
##CC
#Remove CM species
metadata_cc <- droplevels(metadata[!metadata$Species == 'CM',])
#Create a new object for common rownames for CC data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_cc))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cc <- dat.001per[common.rownames,]
metadata_cc <- metadata_cc[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cc), rownames(metadata_cc), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_cc<-decostand(dat_cc, method = "total")

##CM
#Remove CC species
metadata_cm <- droplevels(metadata[!metadata$Species == 'CC',])
#Create a new object for common rownames for CM data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_cm))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cm <- dat.001per[common.rownames,]
metadata_cm <- metadata_cm[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cm), rownames(metadata_cm), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_cm<-decostand(dat_cm, method = "total")

##FT
#Remove H Beach
metadata_FT <- droplevels(metadata_cc[!metadata_cc$Beach == 'Hillsboro',])
#Create a new object for common rownames for CC data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_FT))
#Set the data file and metadata file to have only the data that includes these common names 
dat_FT <- dat.001per[common.rownames,]
metadata_FT <- metadata_cc[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_FT), rownames(metadata_FT), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_FT<-decostand(dat_FT, method = "total")

##H
#Remove FT Beach
metadata_H <- droplevels(metadata_cc[!metadata_cc$Beach == 'Fort Lauderdale',])
#Create a new object for common rownames for CC data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_H))
#Set the data file and metadata file to have only the data that includes these common names 
dat_H <- dat.001per[common.rownames,]
metadata_H <- metadata_cc[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_H), rownames(metadata_H), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_H<-decostand(dat_H, method = "total")

##Cloaca
#Separate both files by sample type
metadata_cloaca <- metadata %>%
  filter(Sample.Type =='Cloaca')
#Create a new object for common rownames for data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_cloaca))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cloaca <- dat.001per[common.rownames,]
metadata_cloaca <- metadata_cloaca[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cloaca), rownames(metadata_cloaca), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_cloaca<-decostand(dat_cloaca, method = "total")

##Hatched Eggs
#Separate both files by sample type
metadata_HatchedOnly <- metadata %>%
  filter(Sample.Type =='Hatched Egg')
#Create a new object for common rownames for data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_HatchedOnly))
#Set the data file and metadata file to have only the data that includes these common names 
dat_HatchedOnly <- dat.001per[common.rownames,]
metadata_HatchedOnly <- metadata_HatchedOnly[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_HatchedOnly), rownames(metadata_HatchedOnly), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_HatchedOnly<-decostand(dat_HatchedOnly, method = "total")

##Unhatched Eggs
#Separate both files by sample type
metadata_UnhatchedOnly <- metadata %>%
  filter(Sample.Type =='Unhatched Egg')
#Create a new object for common rownames for data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_UnhatchedOnly))
#Set the data file and metadata file to have only the data that includes these common names 
dat_UnhatchedOnly <- dat.001per[common.rownames,]
metadata_UnhatchedOnly <- metadata_UnhatchedOnly[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_UnhatchedOnly), rownames(metadata_UnhatchedOnly), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_UnhatchedOnly<-decostand(dat_UnhatchedOnly, method = "total")

##Control Sand
#Separate both files by sample type
metadata_sand <- metadata_cc %>%
  filter(Sample.Type =='Control Sand')
#Create a new object for common rownames for data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_sand))
#Set the data file and metadata file to have only the data that includes these common names 
dat_sand <- dat.001per[common.rownames,]
metadata_sand <- metadata_sand[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_sand), rownames(metadata_sand), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_sand<-decostand(dat_sand, method = "total")

##Nest Sand
#Separate both files by sample type
metadata_nest <- metadata_cc %>%
  filter(Sample.Type =='Nest Sand')
#Create a new object for common rownames for data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_nest))
#Set the data file and metadata file to have only the data that includes these common names 
dat_nest <- dat.001per[common.rownames,]
metadata_nest <- metadata_nest[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_nest), rownames(metadata_nest), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_nest<-decostand(dat_nest, method = "total")


##Hatched and Unhatched Egg
#Remove Cloaca, Nest, Control
metadata_HU <- droplevels(metadata[!metadata$Sample.Type == c('Cloaca','Nest Sand','Control Sand'),])
#Create a new object for common rownames for all HU
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_HU <- dat.001per[common.rownames,]
metadata_HU <- metadata_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_HU), rownames(metadata_HU), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_HU<-decostand(dat_HU, method = "total")

##CC_HU
#Remove CM species
metadata_cc_HU <- droplevels(metadata_HU[!metadata_HU$Species == 'CM',])
#Create a new object for common rownames for CC data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_cc_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cc_HU <- dat.001per[common.rownames,]
metadata_cc_HU <- metadata_cc_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cc_HU), rownames(metadata_cc_HU), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_cc_HU<-decostand(dat_cc_HU, method = "total")

##CM_HU
#Remove CC species
metadata_cm_HU <- droplevels(metadata_HU[!metadata_HU$Species == 'CC',])
#Create a new object for common rownames for CM data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_cm_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cm_HU <- dat.001per[common.rownames,]
metadata_cm_HU <- metadata_cm_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cm_HU), rownames(metadata_cm_HU), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_cm_HU<-decostand(dat_cm_HU, method = "total")

##FT_HU
#Remove Hillsboro Beach
metadata_F_HU <- droplevels(metadata_cc_HU[!metadata_cc_HU$Beach == 'Hillsboro',])
#Create a new object for common rownames for FT data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_F_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_F_HU <- dat.001per[common.rownames,]
metadata_F_HU <- metadata_F_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_F_HU), rownames(metadata_F_HU), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_F_HU<-decostand(dat_F_HU, method = "total")

##H_HU
#Remove Fort Lauderdale Beach
metadata_H_HU <- droplevels(metadata_cc_HU[!metadata_cc_HU$Beach == 'Fort Lauderdale',])
#Create a new object for common rownames for Hillsboro data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_H_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_H_HU <- dat.001per[common.rownames,]
metadata_H_HU <- metadata_H_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_H_HU), rownames(metadata_H_HU), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_H_HU<-decostand(dat_H_HU, method = "total")





##########Alpha Diversity iNEXT (Hill Numbers)
#Calculate the overall totals for each species:
##CC
#Sum data
CC.totals <- colSums(dat_cc)
#Make a matrix for CC data as type "integer"
mCC <- as.matrix(CC.totals)
class(mCC)
##Separate Sample Type for CC
#Merge Sample Types
MergedData <- merge(metadata_cc[1], dat_cc, by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Aggregate Count by Factor (Sample Type)
CCSumSampleType <- aggregate(. ~ Sample.Type, MergedData.Label, sum)
#Make Sample Type row labels
CCSumSampleType.Label <- CCSumSampleType[,-1]
rownames(CCSumSampleType.Label) <- CCSumSampleType[,1]
#Transpose matrix (ASV = rows)
mCC_ST <- as.matrix(t(CCSumSampleType.Label))


##CM
#Sum data
CM.totals <- colSums(dat_cm)
#Make a matrix for CC data as type "integer"
mCM <- as.matrix(CM.totals)
##Separate Sample Type for CM
#Merge Sample Types
MergedData <- merge(metadata_cm[1], dat_cm, by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Aggregate Count by Factor (Sample Type)
CMSumSampleType <- aggregate(. ~ Sample.Type, MergedData.Label, sum)
#Make Sample Type row labels
CMSumSampleType.Label <- CMSumSampleType[,-1]
rownames(CMSumSampleType.Label) <- CMSumSampleType[,1]
#Transpose matrix (ASV = rows)
mCM_ST <- as.matrix(t(CMSumSampleType.Label))

##Overall CC vs CM
#Make a combined matrix list
CC_CM_matrix<-list(mCC,mCM)
#Compute diversity estimates of order q
hill.estimates <- iNEXT(CC_CM_matrix, q=c(0,2), datatype="abundance", nboot=100)
## Examine the results
hill.estimates
#Visualize Results
##Sample-size-based R/E curve
ggiNEXT(hill.estimates, type=1, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 
##Sample completeness curve
ggiNEXT(hill.estimates, type=2, facet.var="none", color.var="site")
##Coverage-based R/E curve 
ggiNEXT(hill.estimates, type=3, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 


##CC by Sample Type
#Compute diversity estimates of order q
hill.estimates.ST.CC <- iNEXT(mCC_ST, q=c(0,2), datatype="abundance", nboot=100)
## Examine the results
hill.estimates.ST.CC
#Visualize Results
##Sample-size-based R/E curve
ggiNEXT(hill.estimates.ST.CC, type=1, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 
##Sample completeness curve
ggiNEXT(hill.estimates.ST.CC, type=2, facet.var="none", color.var="site")
##Coverage-based R/E curve 
ggiNEXT(hill.estimates.ST.CC, type=3, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 

##CM by Sample Type
#Compute diversity estimates of order q
hill.estimates.ST.CM <- iNEXT(mCM_ST, q=c(0,2), datatype="abundance", nboot=100)
## Examine the results
hill.estimates.ST.CM
#Visualize Results
##Sample-size-based R/E curve
ggiNEXT(hill.estimates.ST.CM, type=1, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 
##Sample completeness curve
ggiNEXT(hill.estimates.ST.CM, type=2, facet.var="none", color.var="site")
##Coverage-based R/E curve 
ggiNEXT(hill.estimates.ST.CM, type=3, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 


##FT
#Sum data (for overall comparison)
FT.totals <- colSums(dat_FT)
#Make a matrix for CC data as type "integer"
mFT <- as.matrix(FT.totals)
#Merge Sample Types
MergedData <- merge(metadata_FT[1], dat_FT, by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Aggregate Count by Factor (Sample Type)
FTSumSampleType <- aggregate(. ~ Sample.Type, MergedData.Label, sum)
#Make Sample Type row labels
FTSumSampleType.Label <- FTSumSampleType[,-1]
rownames(FTSumSampleType.Label) <- FTSumSampleType[,1]
#Transpose matrix (ASV = rows)
mFT_ST <- as.matrix(t(FTSumSampleType.Label))


##H
#Sum data (for overall comparison)
H.totals <- colSums(dat_H)
#Make a matrix for CC data as type "integer"
mH <- as.matrix(H.totals)
#Merge Sample Types
MergedData <- merge(metadata_H[1], dat_H, by=0, all=T)
#Fix row names
MergedData.Label <- MergedData[,-1]
rownames(MergedData.Label) <- MergedData[,1]
#Aggregate Count by Factor (Sample Type)
HSumSampleType <- aggregate(. ~ Sample.Type, MergedData.Label, sum)
#Make Sample Type row labels
HSumSampleType.Label <- HSumSampleType[,-1]
rownames(HSumSampleType.Label) <- HSumSampleType[,1]
#Transpose matrix (ASV = rows)
mH_ST <- as.matrix(t(HSumSampleType.Label))

##Overall FT vs H
#Make a combined matrix list
FT_H_matrix <- merge(mFT, mH, by=0, all=T)
#Fix row names
FT_H_matrix.Label <- FT_H_matrix[,-1]
rownames(FT_H_matrix.Label) <- FT_H_matrix[,1]
#Fix column names
colnames(FT_H_matrix.Label) <- c("Fort Lauderdale", "Hillsboro")
#Remove rows with on zeroes
FT_H_matrix.Label.z<-FT_H_matrix.Label[rowSums(FT_H_matrix.Label[])>0,]
#Compute diversity estimates of order q
hill.estimates <- iNEXT(FT_H_matrix.Label.z, q=c(0,2), datatype="abundance", nboot=100)
## Examine the results
hill.estimates
#Visualize Results
##Sample-size-based R/E curve
ggiNEXT(hill.estimates, type=1, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 
##Sample completeness curve
ggiNEXT(hill.estimates, type=2, facet.var="none", color.var="site")
##Coverage-based R/E curve 
ggiNEXT(hill.estimates, type=3, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 


##Fort Lauderdale by Sample Type
#Compute diversity estimates of order q
hill.estimates.ST.FT <- iNEXT(mFT_ST, q=c(0,2), datatype="abundance", nboot=100)
## Examine the results
hill.estimates.ST.FT
#Visualize Results
##Sample-size-based R/E curve
ggiNEXT(hill.estimates.ST.FT, type=1, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 
##Sample completeness curve
ggiNEXT(hill.estimates.ST.FT, type=2, facet.var="none", color.var="site")
##Coverage-based R/E curve 
ggiNEXT(hill.estimates.ST.FT, type=3, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 

###Hillsboro by Sample Type
#Compute diversity estimates of order q
hill.estimates.ST.H <- iNEXT(mH_ST, q=c(0,2), datatype="abundance", nboot=100)
## Examine the results
hill.estimates.ST.H
#Visualize Results
##Sample-size-based R/E curve
ggiNEXT(hill.estimates.ST.H, type=1, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 
##Sample completeness curve
ggiNEXT(hill.estimates.ST.H, type=2, facet.var="none", color.var="site")
##Coverage-based R/E curve 
ggiNEXT(hill.estimates.ST.H, type=3, facet.var="order", color.var="site") + 
  facet_wrap(~order, scales="free") 



##########Alpha Diversity - phyloseq package
#Estimate richness
rich = estimate_richness(physeq_clean)
rich
write.csv(rich, file = "Phyloseq_Diversity_Sample.csv", row.names = TRUE)

#Mean Observed Species
mean(rich$Observed)
##Mean = 250.6091
min(rich$Observed)
##Min= 28
max(rich$Observed)
##Max = 1644
sd(rich$Observed)
##SD = 283.0561

#Check if any ASVs are not present in any samples
any(taxa_sums(physeq_clean) == 0)
##False

#Plot by Sample Type
plot_richness(physeq_clean, x="Sample.Type")+ geom_boxplot()
##Saved as "Alpha Diversity"
#Plot Simplified
plot_richness(physeq_clean, x="Sample.Type", measures = c("Observed","Chao1", "Shannon", "InvSimpson"))+ geom_boxplot()
##Saved as "Alpha Diversity_simplified"
p = plot_richness(physeq_clean, x="Sample.Type", color="Sample.Type", measures = c("Observed","Chao1", "Shannon", "InvSimpson"))
p + geom_boxplot(data=p$data, aes(x = Sample.Type, color = NULL), alpha = 0.1)
##Saved as "Alpha Diversity_colored"




########Beta Diversity
##Overall
#ANOSIM - tests whether distances between groups are greater than within groups.
ano = anosim(dat.ra, metadata$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano
#ANOSIM statistic R: -0.002302 
#Significance: 0.5269 

#Adonis - Analysis of variance using distance matrices
adonis(dat.ra~Sample.Type, data = metadata, method="bray", permutations=9999)
#R2= 0.01353, P-value= 0.869
#1.353% of the sums of squares can be explained by "Sample.Type"
#Not significant

#PerMANOVA to see what sites have the differences
pairwise.perm.manova(dat.ra,metadata$Sample.Type, nperm=9999)
#p-values greater then 0.05 then there is no significant difference between the sites
#All comparisons showed no significant differences
#               Cloaca Control Sand Hatched Egg Nest Sand
# Control Sand  0.91   -            -           -        
# Hatched Egg   0.91   0.91         -           -        
# Nest Sand     0.91   0.91         0.91        -        
# Unhatched Egg 0.91   0.91         0.91        0.91  


#Create NMDS chart 
comm.bc.mds<-metaMDS(dat.ra, distance="bray",trace=FALSE, autotransform=FALSE, k=2, trymax=100)
comm.bc.mds
##Stress = 0.1725899
stressplot(comm.bc.mds)
#Add extra space to right of plot area; change clipping to figure; used to add legend to outside the plot area
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
mds.fig<-ordiplot(comm.bc.mds, display="sites")
ordiellipse(mds.fig, metadata$Sample.Type, label = F, conf = 0.95, col = c("green4","red"))
#Adjust colors, pch=20 make it bullet points
points(mds.fig,"sites", pch=19, col= "red", select = metadata$Sample.Type == "Hatched Egg")
points(mds.fig,"sites", pch=19, col= "cyan", select = metadata$Sample.Type == "Unhatched Egg")
points(mds.fig,"sites", pch=19, col= "blue4", select = metadata$Sample.Type == "Cloaca")
points(mds.fig,"sites", pch=19, col= "green", select = metadata$Sample.Type == "Nest Sand")
points(mds.fig,"sites", pch=19, col= "pink", select = metadata$Sample.Type == "Control Sand")
#Add Stress Value
text(7,4, "Stress = 0.17", cex=0.6)
#Add legend
legend("bottomright",legend= c("Hatched","Unhatched", "Cloaca", "Nest Sand", "Control Sand"), 
       title = "Sample Type",
       col=c("red","cyan", "blue4", "green", "pink"), 
       pch=19, cex=0.55)
##Saved as "NMDS_SampleType"
#####Prefer PRIMER-e for NMDS plot generation



##########SAMPLE TYPE BY TURTLE SPECIES COMPARISON
#Two-Way ANOSIM
SampleType.Species.test <- anosim(dat.ra, metadata$Sample.Type, strata = metadata$Species,
                           permutations = 9999)
Species.SampleType.test <- anosim(dat.ra, metadata$Species, strata = metadata$Sample.Type,
                           permutations = 9999)
###Tests for differences between Sample Types within Species
SampleType.Species.test
#ANOSIM statistic R: -0.002302 
#Significance: 0.5154 

###Tests for differences between Species within Sample Types
Species.SampleType.test
#ANOSIM statistic R: 0.2641 
#Significance: 0.001 

#ANOSIM - tests whether distances between groups are greater than within groups.
ano_cc = anosim(dat.ra_cc, metadata_cc$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_cc
#ANOSIM statistic R: 0.4925  
#Significance: 1e-04 

ano_cm = anosim(dat.ra_cm, metadata_cm$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_cm
#ANOSIM statistic R: 0.1536 
#Significance: 0.0061 

#Adonis - Analysis of variance using distance matrices
adonis(dat.ra_cc~Sample.Type, data = metadata_cc, permutations = 9999, method = "bray")
#R2= 0.19036, P-value= 0.001
#19.04% of the sums of squares can be explained by "Sample.Type"
#Significant

adonis(dat.ra_cm~Sample.Type, data = metadata_cm, permutations = 9999, method = "bray")
#R2= 0.27328, P-value= 0.001
#27.33% of the sums of squares can be explained by "Sample.Type"
#Significant

#PerMANOVA to see what sites have the differences
pairwise.perm.manova(dat.ra_cc,metadata_cc$Sample.Type, nperm=9999)
#p-values greater then 0.05 then there is no significant difference between the sites
#All comparisons show significant differences
#               Cloaca Control Sand Hatched Egg Nest Sand
# Control Sand  0.001  -            -           -        
# Hatched Egg   0.001  0.001        -           -        
# Nest Sand     0.001  0.001        0.001       -        
# Unhatched Egg 0.001  0.001        0.001       0.001   

pairwise.perm.manova(dat.ra_cm,metadata_cm$Sample.Type, nperm=9999)
#p-values greater then 0.05 then there is no significant difference between the sites
#All comparisons show significant differences, except Nest vs Unhatched Eggs
#               Cloaca Control Sand Hatched Egg Nest Sand
# Control Sand  0.0025 -            -           -        
# Hatched Egg   0.0025 0.0171       -           -        
# Nest Sand     0.0050 0.0040       0.0333      -        
# Unhatched Egg 0.0025 0.0287       0.0025      0.0840



############Cloaca Comparison
#ANOSIM - tests whether distances between groups are greater than within groups.
ano_cloaca = anosim(dat_cloaca, metadata_cloaca$Species, permutations = 9999, distance = "bray", strata = NULL)
ano_cloaca
#ANOSIM statistic R: 0.5613 
#Significance: 1e-04 

#SIMPER
dat.simp<-simper(dat_cloaca, metadata_cloaca$Species, permutations = 9999)
sink("Simper_Cloaca.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

#SIMPER Cloaca Growth
dat.simp<-simper(dat_cloaca, metadata_cloaca$Cloaca.Growth, permutations = 9999)
sink("Simper_CloacaGrowth.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output





##########SAMPLE TYPE BY BEACH COMPARISON
#Two-Way ANOSIM
SampleType.Beach.test <- anosim(dat.ra_cc, metadata_cc$Sample.Type, strata = metadata_cc$Beach,
                                permutations = 9999)
Beach.SampleType.test <- anosim(dat.ra_cc, metadata_cc$Beach, strata = metadata_cc$Sample.Type,
                                permutations = 9999)
###Tests for differences between Sample Types within Beaches
SampleType.Beach.test
#ANOSIM statistic R: 0.4925 
#Significance: 1e-04 

###Tests for differences between Beaches within Sample Types
Beach.SampleType.test
#ANOSIM statistic R: 0.09194 
#Significance: 1e-04

#ANOSIM - tests whether distances between groups are greater than within groups.
ano_H = anosim(dat.ra_H, metadata_H$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_H
#ANOSIM statistic R: 0.6586 
#Significance: 1e-04 

ano_F = anosim(dat.ra_FT, metadata_FT$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_F
#ANOSIM statistic R: 0.3385 
#Significance: 1e-04 

#Adonis - Analysis of variance using distance matrices
adonis(dat.ra_H~Sample.Type, data = metadata_H, permutations = 9999, method = "bray")
#R2= 0.25184, P-value= 0.001
#25.18% of the sums of squares can be explained by "Sample.Type"
#Significant

adonis(dat.ra_FT~Sample.Type, data = metadata_FT, permutations = 9999, method = "bray")
#R2= 0.20363, P-value= 0.001
#20.36% of the sums of squares can be explained by "Sample.Type"
#Significant

#PerMANOVA to see what sites have the differences
pairwise.perm.manova(dat.ra_H,metadata_H$Sample.Type, nperm=9999)
#p-values greater then 0.05 then there is no significant difference between the sites
#All comparisons show significant differences
#               Cloaca Control Sand Hatched Egg Nest Sand
# Control Sand  0.001  -            -           -        
# Hatched Egg   0.001  0.001        -           -        
# Nest Sand     0.001  0.001        0.001       -        
# Unhatched Egg 0.001  0.001        0.001       0.001   

pairwise.perm.manova(dat.ra_FT,metadata_FT$Sample.Type, nperm=9999)
#p-values greater then 0.05 then there is no significant difference between the sites
#All comparisons show significant differences
#               Cloaca Control Sand Hatched Egg Nest Sand
# Control Sand  0.0012 -            -           -        
# Hatched Egg   0.0012 0.0012       -           -        
# Nest Sand     0.0012 0.0170       0.0100      -        
# Unhatched Egg 0.0012 0.0012       0.0012      0.0012




############Control Sand Comparsion
#ANOSIM - tests whether distances between groups are greater than within groups.
ano_sand = anosim(dat.ra_sand, metadata_sand$Beach, permutations = 9999, distance = "bray", strata = NULL)
ano_sand
#ANOSIM statistic R: 0.2738 
#Significance: 0.0118 

#SIMPER SAND
dat.simp<-simper(dat.ra_sand, metadata_sand$Beach, permutations = 9999)
sink("Simper_Sand.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output


############Nest Sand Comparsion
#ANOSIM - tests whether distances between groups are greater than within groups.
ano_nest = anosim(dat.ra_nest, metadata_nest$Beach, permutations = 9999, distance = "bray", strata = NULL)
ano_nest
#ANOSIM statistic R: 0.158
#Significance: 0.0287

##Nest Distances Correlation
#Create distance matrices
dat.bc.dist_sand<-vegdist(dat.ra_sand, method = "bray")
dat.bc.dist_nest<-vegdist(dat.ra_nest, method = "bray")
##Export and create data matrices in Excel
Sand_micro_dist <- read.delim("~/Grad School/Thesis/Data/Sand_micro_dist.txt", row.names=1)
Nest_micro_dist <- read.delim("~/Grad School/Thesis/Data/Nest_micro_dist.txt", row.names=1)
#Import Nest Distances
Nest.Distances_CC <- read.delim("~/Grad School/Thesis/Data/Nest Distances_CC.txt", row.names=1)
#Mantel Test - Correlation between Two Distance Matrices
mantel(Sand_micro_dist, Nest.Distances_CC, permutations = 9999)
mantel(Nest_micro_dist, Nest.Distances_CC, permutations = 9999)


###########HATCHED VS UNHATCHED EGG COMPARISON ALL
#ANOSIM - tests whether distances between groups are greater than within groups.
ano_HU = anosim(dat.ra_HU, metadata_HU$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_HU
#ANOSIM statistic R: 0.247  
#Significance: 1e-04 

#Adonis - Analysis of variance using distance matrices
adonis(dat.ra_HU~Sample.Type, data = metadata_HU, permutations = 9999, method = "bray")
#R2= 0.07794, P-value= 0.001
#7.794% of the sums of squares can be explained by "Sample.Type"
#Significant

#PerMANOVA to see what sites have the differences
pairwise.perm.manova(dat.ra_HU,metadata_HU$Sample.Type, nperm=9999)
#              Hatched Egg
#Unhatched Egg 0.001
###Significant differences between sample type

#SIMPER
dat.simp<-simper(dat.ra_HU, metadata_HU$Sample.Type, permutations = 9999)
sink("Simper_HU.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

#PERMDISP - see vairance within groups
dat.bc.dist<-vegdist(dat.ra_HU, method = "bray")
#Calculate multivariate dispersions
mod_all<- betadisper(dat.bc.dist, metadata_HU$Sample.Type)
mod_all
#Perform test
anova(mod_all)
#Permutation test for F
pmod_all <- permutest(mod_all, permutations = 9999, pairwise = TRUE)
pmod_all
#Tukey's Honest Significant Differences
(mod.all.HSD <- TukeyHSD(mod_all))
plot(mod.all.HSD)
##Permustats() method
pstat_all <- permustats(pmod_all)
densityplot(pstat_all, scales = list(x = list(relation = "free")))
qqmath(pstat_all, scales = list(relation = "free"))
#Plot
boxplot(mod_all)




###########HATCHED VS UNHATCHED EGG COMPARISON BY SPECIES
#ANOSIM - tests whether distances between groups are greater than within groups.
ano_cc_HU = anosim(dat.ra_cc_HU, metadata_cc_HU$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_cc_HU
#ANOSIM statistic R: 0.3449   
#Significance: 1e-04 

ano_cm_HU = anosim(dat.ra_cm_HU, metadata_cm_HU$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_cm_HU
#ANOSIM statistic R: 0.1322    
#Significance: 4e-04 

#SIMPER
dat.simp<-simper(dat.ra_cc_HU, metadata_cc_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_CC.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

dat.simp<-simper(dat.ra_cm_HU, metadata_cm_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_CM.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output


#PERMDISP - see vairance within groups
dat.bc.dist_cc<-vegdist(dat.ra_cc_HU, method = "bray")
dat.bc.dist_cm<-vegdist(dat.ra_cm_HU, method = "bray")
#Calculate multivariate dispersions
mod_cc <- betadisper(dat.bc.dist_cc, metadata_cc_HU$Sample.Type)
mod_cc
mod_cm <- betadisper(dat.bc.dist_cm, metadata_cm_HU$Sample.Type)
mod_cm
#Perform test
anova(mod_cc)
anova(mod_cm)
#Permutation test for F
pmod_cc <- permutest(mod_cc, permutations = 9999, pairwise = TRUE)
pmod_cm <- permutest(mod_cm, permutations = 9999, pairwise = TRUE)

#Tukey's Honest Significant Differences
(mod.cc.HSD <- TukeyHSD(mod_cc))
plot(mod.cc.HSD)
(mod.cm.HSD <- TukeyHSD(mod_cm))
plot(mod.cm.HSD)
##Permustats() method
pstat_cc <- permustats(pmod_cc)
densityplot(pstat_cc, scales = list(x = list(relation = "free")))
qqmath(pstat_cc, scales = list(relation = "free"))

pstat_cm <- permustats(pmod_cm)
densityplot(pstat_cm, scales = list(x = list(relation = "free")))
qqmath(pstat_cm, scales = list(relation = "free"))
#Plot
boxplot(mod_cc)
boxplot(mod_cm)





###########HATCHED VS UNHATCHED EGG COMPARISON BY BEACH (only CC)
#ANOSIM - tests whether distances between groups are greater than within groups.
ano_F_HU = anosim(dat.ra_F_HU, metadata_F_HU$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_F_HU
#ANOSIM statistic R: 0.2617
#Significance: 1e-04 

ano_H_HU = anosim(dat.ra_H_HU, metadata_H_HU$Sample.Type, permutations = 9999, distance = "bray", strata = NULL)
ano_H_HU
#ANOSIM statistic R: 0.4207   
#Significance: 1e-04 

#SIMPER
dat.simp<-simper(dat.ra_F_HU, metadata_F_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_F.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

dat.simp<-simper(dat.ra_H_HU, metadata_H_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_H.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output




##########ENVIRONMENTAL DATA CORRELATION
######ALL Hatched and Unhatched
#Make sure all environmental factors are quantitative (change yes/no to 1/0)
metadata_HU$Recent.Renourishment<-ifelse(metadata_HU$Recent.Renourishment=="Yes",1,0)
metadata_HU$Roots<-ifelse(metadata_HU$Roots=="Yes",1,0)
metadata_HU$Washover<-ifelse(metadata_HU$Washover=="Yes",1,0)

#CCA, significance of the enivornmental factors on diversity 
set.seed(55);env.cca<-cca(dat.ra_HU~Latitude+Longitude+pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+R.Zone+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_HU)
vif.cca(env.cca)
#If VIF higher than 20 remove factor with highest VIF (variance inflation factor)
#Remove R Zone
set.seed(55);env.cca<-cca(dat.ra_HU~Latitude+Longitude+pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_HU)
vif.cca(env.cca)
#Remove Latitude
set.seed(55);env.cca<-cca(dat.ra_HU~Longitude+pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_HU)
vif.cca(env.cca)

#Zero the variables
set.seed(55);lwr<- cca(dat.ra_HU~1, data=metadata_HU)
lwr
#Total Inertia = total variance in species (observdistributions
#Unconstrained Inertia = the variance explained by the environmental variables
#Using a forward selecting model, must keep our set seed 
set.seed(55);mods.all<-ordiR2step(lwr, scope = formula(env.cca))
mods.all
vif.cca(mods.all)
R2.adj.all<-RsquareAdj(mods.all)
R2.adj.all
mods.all$anova
#Add extra space to right of plot area; change clipping to figure; used to add legend to outside the plot area
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#Plot CCA
cca.p <- plot(mods.all,type = "none")
#Adjust plot point shapes and colors
points(cca.p,"sites", pch=17, col= "red", select = metadata_HU$Sample.Type == "Hatched Egg")
points(cca.p,"sites", pch=19, col= "cyan", select = metadata_HU$Sample.Type == "Unhatched Egg")
ef.all<- envfit(cca.p,metadata_HU[,c("Longitude")])
plot(ef.all)
#Add legend, click on right side of plot to make it outside the graphing areas, cex used to shrink legend
legend(locator(1),legend= as.character(paste(" ", unique(metadata_HU$Sample.Type))), col=c("red","cyan"), pch=c(17,19), cex=0.70)
###Saved as "CCA_All"


######CCA Without Location Information 
set.seed(55);env.cca<-cca(dat.ra_HU~pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_HU)
vif.cca(env.cca)
#If VIF higher than 20 remove factor with highest VIF (variance inflation factor)

#Zero the variables
set.seed(55);lwr<- cca(dat.ra_HU~1, data=metadata_HU)
lwr
#Total Inertia = total variance in species (observdistributions
#Unconstrained Inertia = the variance explained by the environmental variables
#Using a forward selecting model, must keep our set seed 
set.seed(55);mods.all<-ordiR2step(lwr, scope = formula(env.cca))
mods.all
vif.cca(mods.all)
R2.adj.all<-RsquareAdj(mods.all)
R2.adj.all
mods.all$anova
#Add extra space to right of plot area; change clipping to figure; used to add legend to outside the plot area
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#Plot CCA
cca.p <- plot(mods.all,type = "none")
#Adjust plot point shapes and colors
points(cca.p,"sites", pch=17, col= "red", select = metadata_HU$Sample.Type == "Hatched Egg")
points(cca.p,"sites", pch=19, col= "cyan", select = metadata_HU$Sample.Type == "Unhatched Egg")
ef.all<- envfit(cca.p,metadata_HU[,c("Dune.Distance", "Incubation.Length", "Temperature.Side")])
plot(ef.all)
#Add legend, click on right side of plot to make it outside the graphing areas, cex used to shrink legend
legend(locator(1),legend= as.character(paste(" ", unique(metadata_HU$Sample.Type))), col=c("red","cyan"), pch=c(17,19), cex=0.70)
###Saved as "CCA_All_noLocation"


######CCA CC Species Only
#Turn yes/no to 1/0
metadata_cc_HU$Recent.Renourishment<-ifelse(metadata_cc_HU$Recent.Renourishment=="Yes",1,0)
metadata_cc_HU$Roots<-ifelse(metadata_cc_HU$Roots=="Yes",1,0)
metadata_cc_HU$Washover<-ifelse(metadata_cc_HU$Washover=="Yes",1,0)

set.seed(55);env.cca<-cca(dat.ra_cc_HU~pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_cc_HU)
vif.cca(env.cca)
#If VIF higher than 10 remove factor with highest VIF (variance inflation factor)
#Remove Temperature.Side
set.seed(55);env.cca<-cca(dat.ra_cc_HU~pH.Side+pH.Bottom+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_cc_HU)
vif.cca(env.cca)
#Remove pH.Side
set.seed(55);env.cca<-cca(dat.ra_cc_HU~pH.Bottom+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_cc_HU)
vif.cca(env.cca)

#Zero the variables
set.seed(55);lwr<- cca(dat.ra_cc_HU~1, data=metadata_cc_HU)
lwr
#Total Inertia = total variance in species (observdistributions
#Unconstrained Inertia = the variance explained by the environmental variables
#Using a forward selecting model, must keep our set seed 
set.seed(55);mods.all<-ordiR2step(lwr, scope = formula(env.cca))
mods.all
vif.cca(mods.all)
R2.adj.all<-RsquareAdj(mods.all)
R2.adj.all
mods.all$anova
#Add extra space to right of plot area; change clipping to figure; used to add legend to outside the plot area
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#Plot CCA
cca.p <- plot(mods.all,type = "none")
#Adjust plot point shapes and colors
points(cca.p,"sites", pch=17, col= "red", select = metadata_cc_HU$Sample.Type == "Hatched Egg")
points(cca.p,"sites", pch=19, col= "cyan", select = metadata_cc_HU$Sample.Type == "Unhatched Egg")
ef.all<- envfit(cca.p,metadata_cc_HU[,c("High.Tide.Distance", "Conductivity.Bottom", "Hatch.Success", "Incubation.Length", "Clutch.Size", "Sand.Grain.Size", "pH.Bottom", "Chamber.Depth", "Temperature.Bottom", "Sorting.Coefficient")])
plot(ef.all, cex=0.7)
#Add legend, click on right side of plot to make it outside the graphing areas, cex used to shrink legend
legend(locator(1),legend= as.character(paste(" ", unique(metadata_cc_HU$Sample.Type))), col=c("red","cyan"), pch=c(17,19), cex=0.70)
###Saved as "CCA_CC"
#To see where each sample is located (to determine which nests are impacted)
cca.p[["sites"]]


######CCA CM Species Only
#Turn yes/no to 1/0
metadata_cm_HU$Recent.Renourishment<-ifelse(metadata_cm_HU$Recent.Renourishment=="Yes",1,0)
metadata_cm_HU$Roots<-ifelse(metadata_cm_HU$Roots=="Yes",1,0)
metadata_cm_HU$Washover<-ifelse(metadata_cm_HU$Washover=="Yes",1,0)

set.seed(55);env.cca<-cca(dat.ra_cm_HU~pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_cm_HU)
vif.cca(env.cca)
#If VIF higher than 10 remove factor with highest VIF (variance inflation factor)
#Remove pH.Side
set.seed(55);env.cca<-cca(dat.ra_cm_HU~pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_cm_HU)
vif.cca(env.cca)
#Remove Temperature.Side
set.seed(55);env.cca<-cca(dat.ra_cm_HU~pH.Bottom+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_cm_HU)
vif.cca(env.cca)
#Remove Sand.Grain.Size
set.seed(55);env.cca<-cca(dat.ra_cm_HU~pH.Bottom+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_cm_HU)
vif.cca(env.cca)
#Remove Incubation.Length
set.seed(55);env.cca<-cca(dat.ra_cm_HU~pH.Bottom+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sorting.Coefficient+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_cm_HU)
vif.cca(env.cca)
#Remove NAs
set.seed(55);env.cca<-cca(dat.ra_cm_HU~pH.Bottom+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sorting.Coefficient+Clutch.Size, data=metadata_cm_HU)
vif.cca(env.cca)

#Zero the variables
set.seed(55);lwr<- cca(dat.ra_cm_HU~1, data=metadata_cm_HU)
lwr
#Total Inertia = total variance in species (observdistributions
#Unconstrained Inertia = the variance explained by the environmental variables
#Using a forward selecting model, must keep our set seed 
set.seed(55);mods.all<-ordiR2step(lwr, scope = formula(env.cca))
mods.all
vif.cca(mods.all)
R2.adj.all<-RsquareAdj(mods.all)
R2.adj.all
mods.all$anova
#Add extra space to right of plot area; change clipping to figure; used to add legend to outside the plot area
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#Plot CCA
cca.p <- plot(mods.all,type = "none")
#Adjust plot point shapes and colors
points(cca.p,"sites", pch=17, col= "red", select = metadata_cm_HU$Sample.Type == "Hatched Egg")
points(cca.p,"sites", pch=19, col= "cyan", select = metadata_cm_HU$Sample.Type == "Unhatched Egg")
ef.all<- envfit(cca.p,metadata_cm_HU[,c("Temperature.Bottom", "pH.Bottom")])
plot(ef.all)
#Add legend, click on right side of plot to make it outside the graphing areas, cex used to shrink legend
legend(locator(1),legend= as.character(paste(" ", unique(metadata_cm_HU$Sample.Type))), col=c("red","cyan"), pch=c(17,19), cex=0.70)
###Saved as "CCA_CM"
#To see where each sample is located (to determine which nests are impacted)
cca.p[["sites"]]



######CCA Fort Lauderdale Beach Only
#Turn yes/no to 1/0
metadata_F_HU$Recent.Renourishment<-ifelse(metadata_F_HU$Recent.Renourishment=="Yes",1,0)
metadata_F_HU$Roots<-ifelse(metadata_F_HU$Roots=="Yes",1,0)
metadata_F_HU$Washover<-ifelse(metadata_F_HU$Washover=="Yes",1,0)

set.seed(55);env.cca<-cca(dat.ra_F_HU~pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_F_HU)
vif.cca(env.cca)
#If VIF higher than 10 remove factor with highest VIF (variance inflation factor)
#Remove pH.Bottom
set.seed(55);env.cca<-cca(dat.ra_F_HU~pH.Side+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_F_HU)
vif.cca(env.cca)
#Remove Hatch.Success
set.seed(55);env.cca<-cca(dat.ra_F_HU~pH.Side+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_F_HU)
vif.cca(env.cca)
#Remove Temperature.Side
set.seed(55);env.cca<-cca(dat.ra_F_HU~pH.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_F_HU)
vif.cca(env.cca)
#Remove Clutch.Size
set.seed(55);env.cca<-cca(dat.ra_F_HU~pH.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_F_HU)
vif.cca(env.cca)
#Remove Dune.Distance
set.seed(55);env.cca<-cca(dat.ra_F_HU~pH.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Chamber.Depth+High.Tide.Distance, data=metadata_F_HU)
vif.cca(env.cca)


#Zero the variables
set.seed(55);lwr<- cca(dat.ra_F_HU~1, data=metadata_F_HU)
lwr
#Total Inertia = total variance in species (observdistributions)
#Unconstrained Inertia = the variance explained by the environmental variables
#Using a forward selecting model, must keep our set seed 
set.seed(55);mods.all<-ordiR2step(lwr, scope = formula(env.cca))
mods.all
vif.cca(mods.all)
R2.adj.all<-RsquareAdj(mods.all)
R2.adj.all
mods.all$anova
#Add extra space to right of plot area; change clipping to figure; used to add legend to outside the plot area
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#Plot CCA
cca.p <- plot(mods.all,type = "none")
#Adjust plot point shapes and colors
points(cca.p,"sites", pch=17, col= "red", select = metadata_F_HU$Sample.Type == "Hatched Egg")
points(cca.p,"sites", pch=19, col= "cyan", select = metadata_F_HU$Sample.Type == "Unhatched Egg")
ef.all<- envfit(cca.p,metadata_F_HU[,c("High.Tide.Distance", "Conductivity.Bottom", "Sand.Grain.Size", "Incubation.Length", "pH.Side")])
plot(ef.all)
#Add legend, click on right side of plot to make it outside the graphing areas, cex used to shrink legend
legend(locator(1),legend= as.character(paste(" ", unique(metadata_F_HU$Sample.Type))), col=c("red","cyan"), pch=c(17,19), cex=0.70)
###Saved as "CCA_F"
#To see where each sample is located (to determine which nests are impacted)
cca.p[["sites"]]


######CCA Hillsboro Beach Only
#Turn yes/no to 1/0
metadata_H_HU$Recent.Renourishment<-ifelse(metadata_H_HU$Recent.Renourishment=="Yes",1,0)
metadata_H_HU$Roots<-ifelse(metadata_H_HU$Roots=="Yes",1,0)
metadata_H_HU$Washover<-ifelse(metadata_H_HU$Washover=="Yes",1,0)

set.seed(55);env.cca<-cca(dat.ra_H_HU~pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_H_HU)
vif.cca(env.cca)
#If VIF higher than 10 remove factor with highest VIF (variance inflation factor)
#Remove Temperature.Bottom
set.seed(55);env.cca<-cca(dat.ra_H_HU~pH.Side+pH.Bottom+Temperature.Side+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_H_HU)
vif.cca(env.cca)
#Remove NAs
set.seed(55);env.cca<-cca(dat.ra_H_HU~pH.Side+pH.Bottom+Temperature.Side+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length, data=metadata_H_HU)
vif.cca(env.cca)

#Zero the variables
set.seed(55);lwr<- cca(dat.ra_H_HU~1, data=metadata_H_HU)
lwr
#Total Inertia = total variance in species (observdistributions)
#Unconstrained Inertia = the variance explained by the environmental variables
#Using a forward selecting model, must keep our set seed 
set.seed(55);mods.all<-ordiR2step(lwr, scope = formula(env.cca))
mods.all
vif.cca(mods.all)
R2.adj.all<-RsquareAdj(mods.all)
R2.adj.all
mods.all$anova
#Add extra space to right of plot area; change clipping to figure; used to add legend to outside the plot area
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#Plot CCA
cca.p <- plot(mods.all,type = "none")
#Adjust plot point shapes and colors
points(cca.p,"sites", pch=17, col= "red", select = metadata_H_HU$Sample.Type == "Hatched Egg")
points(cca.p,"sites", pch=19, col= "cyan", select = metadata_H_HU$Sample.Type == "Unhatched Egg")
ef.all<- envfit(cca.p,metadata_H_HU[,c("Conductivity.Side", "Sand.Grain.Size", "pH.Side", "Temperature.Side", "pH.Bottom", "Incubation.Length", "Conductivity.Side", "Sorting.Coefficient")])
plot(ef.all)
#Add legend, click on right side of plot to make it outside the graphing areas, cex used to shrink legend
legend(locator(1),legend= as.character(paste(" ", unique(metadata_H_HU$Sample.Type))), col=c("red","cyan"), pch=c(17,19), cex=0.70)
###Saved as "CCA_H"
#To see where each sample is located (to determine which nests are impacted)
cca.p[["sites"]]



###########BEST/BIOENV ANALYSIS - Compare to CCA results
######BEST CC Species Only
metadata_cc_HU_bioenv <- subset(metadata_cc_HU, select = c(8,9,12,13,20,21,23,24,28,29,31,32,34,35,37,38))
#Normalize metadata
metadata_cc_HU_bioenv.scaled <- scale(metadata_cc_HU_bioenv)
#Run Test
best.test <- bioenv(dat.ra_cc_HU, metadata_cc_HU_bioenv.scaled, trace=TRUE, index = "bray", metric = "euclidean", method = "spearman")
best.test
summary(best.test)
#Mantel Test to estimate significance of best model
tested.dists <- bioenvdist(best.test)
spp.dist <- vegdist(dat.ra_cc_HU, method = "bray")
sig.test <- mantel(tested.dists, spp.dist, method = "spearman", permutations = 9999)
sig.test

######BEST CM Species Only
metadata_cm_HU_bioenv <- subset(metadata_cm_HU, select = c(8,9,12,13,20,21,23,24,28,29,31,32,34,35,37,38))
#Normalize metadata
metadata_cm_HU_bioenv.scaled <- scale(metadata_cm_HU_bioenv)
#Run Test
best.test <- bioenv(dat.ra_cm_HU, metadata_cm_HU_bioenv.scaled, trace=TRUE, index = "bray", metric = "euclidean", method = "spearman")
best.test
summary(best.test)
#Mantel Test to estimate significance of best model
tested.dists <- bioenvdist(best.test)
spp.dist <- vegdist(dat.ra_cm_HU, method = "bray")
sig.test <- mantel(tested.dists, spp.dist, method = "spearman", permutations = 9999)
sig.test

######BEST FT Beach Only
metadata_F_HU_bioenv <- subset(metadata_F_HU, select = c(8,9,12,13,20,21,23,24,28,29,31,32,34,35,37,38))
#Normalize metadata
metadata_F_HU_bioenv.scaled <- scale(metadata_F_HU_bioenv)
#Run Test
best.test <- bioenv(dat.ra_F_HU, metadata_F_HU_bioenv.scaled, trace=TRUE, index = "bray", metric = "euclidean", method = "spearman")
best.test
summary(best.test)
#Mantel Test to estimate significance of best model
tested.dists <- bioenvdist(best.test)
spp.dist <- vegdist(dat.ra_F_HU, method = "bray")
sig.test <- mantel(tested.dists, spp.dist, method = "spearman", permutations = 9999)
sig.test

######BEST H Beach Only
metadata_H_HU_bioenv <- subset(metadata_H_HU, select = c(8,9,12,13,20,21,23,24,28,29,31,32,34,35,37,38))
#Normalize metadata
metadata_H_HU_bioenv.scaled <- scale(metadata_H_HU_bioenv)
#Run Test
best.test <- bioenv(dat.ra_H_HU, metadata_H_HU_bioenv.scaled, trace=TRUE, index = "bray", metric = "euclidean", method = "spearman")
best.test
summary(best.test)
#Mantel Test to estimate significance of best model
tested.dists <- bioenvdist(best.test)
spp.dist <- vegdist(dat.ra_H_HU, method = "bray")
sig.test <- mantel(tested.dists, spp.dist, method = "spearman", permutations = 9999)
sig.test

######BEST Control Sand
metadata_sand_bioenv <- subset(metadata_sand, select = c(8,9,12,13,20,21,23,24,28,29,31,32,34,35,37,38))
#Normalize metadata
metadata_sand_bioenv.scaled <- scale(metadata_sand_bioenv)
#Run Test
best.test <- bioenv(dat.ra_sand, metadata_sand_bioenv.scaled, trace=TRUE, index = "bray", metric = "euclidean", method = "spearman")
best.test
summary(best.test)
#Mantel Test to estimate significance of best model
tested.dists <- bioenvdist(best.test)
spp.dist <- vegdist(dat.ra_sand, method = "bray")
sig.test <- mantel(tested.dists, spp.dist, method = "spearman", permutations = 9999)
sig.test

######BEST Nest Sand
metadata_N_bioenv <- subset(metadata_N, select = c(8,9,12,13,20,21,23,24,28,29,31,32,34,35,37,38))
#Normalize metadata
metadata_N_bioenv.scaled <- scale(metadata_N_bioenv)
#Run Test
best.test <- bioenv(dat.ra_N, metadata_N_bioenv.scaled, trace=TRUE, index = "bray", metric = "euclidean", method = "spearman")
best.test
summary(best.test)
#Mantel Test to estimate significance of best model
tested.dists <- bioenvdist(best.test)
spp.dist <- vegdist(dat.ra_N, method = "bray")
sig.test <- mantel(tested.dists, spp.dist, method = "spearman", permutations = 9999)
sig.test



###########HATCHED EGG COMPARISON
###Species Comparison
##ANOSIM - Compare Species
ano_Hatch_sp = anosim(dat.ra_HatchedOnly, metadata_HatchedOnly$Species, permutations = 9999, distance = "bray", strata = NULL)
ano_Hatch_sp
#ANOSIM statistic R: -0.1175  
#Significance: 0.9788
#####No differences between hatched eggs by Species

#SIMPER to find potential indicators
dat.simp<-simper(dat.ra_HatchedOnly, metadata_HatchedOnly$Species, permutations = 9999)
sink("Simper_HatchedEggs_Species.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

###Beach Comparison
#Remove CM species
metadata_HatchedOnly_cc <- droplevels(metadata_HatchedOnly[!metadata_HatchedOnly$Species == 'CM',])
#Create a new object for common rownames for CC data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_HatchedOnly_cc))
#Set the data file and metadata file to have only the data that includes these common names 
dat_HatchedOnly_cc <- dat.001per[common.rownames,]
metadata_HatchedOnly_cc <- metadata_HatchedOnly_cc[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_HatchedOnly_cc), rownames(metadata_HatchedOnly_cc), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_HatchedOnly_cc<-decostand(dat_HatchOnly_cc, method = "total")

##ANOSIM - Compare Beach
ano_Hatch_beach = anosim(dat.ra_HatchOnly_cc, metadata_HatchOnly_cc$Beach, permutations = 9999, distance = "bray", strata = NULL)
ano_Hatch_beach
#ANOSIM statistic R: 0.08813 
#Significance: 0.0021
#####Significant differences between hatched eggs by Beach

#SIMPER to find potential indicators
dat.simp<-simper(dat.ra_HatchOnly_cc, metadata_HatchOnly_cc$Beach, permutations = 9999)
sink("Simper_HatchedEggsCC_Beach.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output



###Nest Comparison
ano_Hatch_nestnumber = anosim(dat.ra_HatchedOnly, metadata_HatchedOnly$Nest.Number, permutations = 9999, distance = "bray", strata = NULL)
ano_Hatch_nestnumber
#ANOSIM statistic R: 0.5575 
#Significance: 1e-04 

#Adonis - Analysis of variance using distance matrices
adonis(dat.ra_HatchedOnly, data = metadata_HatchedOnly, permutations = 9999, method = "bray")
#R2= 0.12988, P-value= 0.001
#12.99% of the sums of squares can be explained significantly by "Nest.Number"

#PerMANOVA to see what sites have the differences
pairwise.perm.manova(dat.ra_HatchedOnly,metadata_HatchedOnly$Nest.Number)
#p-values greater then 0.05 then there is no significant difference between the sites
#All comparisons showed no significant differences##PerMANOVA to see what sites have the differences

#CCA, significance of the enivornmental factors on diversity 
set.seed(55);env.cca<-cca(dat.ra_HatchedOnly~Latitude+Longitude+pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+R.Zone+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_HatchOnly)
vif.cca(env.cca)
#If VIF higher than 20 remove factor with highest VIF (variance inflation factor)
#Remove R Zone
set.seed(55);env.cca<-cca(dat.ra_HatchedOnly~Latitude+Longitude+pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_HatchOnly)
vif.cca(env.cca)
#Remove Latitude
set.seed(55);env.cca<-cca(dat.ra_HatchedOnly~Longitude+pH.Side+pH.Bottom+Temperature.Side+Temperature.Bottom+Conductivity.Side+Conductivity.Bottom+Sand.Grain.Size+Sorting.Coefficient+Incubation.Length+Clutch.Size+Hatch.Success+Chamber.Depth+High.Tide.Distance+Dune.Distance, data=metadata_HatchOnly)
vif.cca(env.cca)
#Zero the variables
set.seed(55);lwr<- cca(dat.ra_HatchedOnly~1, data=metadata_HatchOnly)
lwr
#Total Inertia = total variance in species (observdistributions
#Unconstrained Inertia = the variance explained by the environmental variables
#Using a forward selecting model, must keep our set seed 
set.seed(55);mods.all<-ordiR2step(lwr, scope = formula(env.cca))
mods.all
vif.cca(mods.all)
R2.adj.all<-RsquareAdj(mods.all)
R2.adj.all
mods.all$anova


#BEST Analysis Hatched Only
metadata_HatchedOnly_bioenv <- subset(metadata_HatchedOnly, select = c(8,9,12,13,20,21,23,24,28,29,31,32,34,35,37,38))
#Normalize metadata
metadata_HatchedOnly_bioenv.scaled <- scale(metadata_HatchedOnly_bioenv)
#Run Test
best.test <- bioenv(dat.ra_HatchedOnly, metadata_HatchedOnly_bioenv.scaled, trace=TRUE, index = "bray", metric = "euclidean", method = "spearman")
best.test
summary(best.test)
#Mantel Test to estimate significance of best model
tested.dists <- bioenvdist(best.test)
spp.dist <- vegdist(dat.ra_HatchedOnly, method = "bray")
sig.test <- mantel(tested.dists, spp.dist, method = "spearman", permutations = 9999)
sig.test


###########UNHATCHED EGG COMPARISON
###Species Comparison
##ANOSIM - Compare Species
ano_Unhatched_sp = anosim(dat.ra_UnhatchedOnly, metadata_UnhatchedOnly$Species, permutations = 9999, distance = "bray", strata = NULL)
ano_Unhatched_sp
#ANOSIM statistic R: -0.06869  
#Significance: 0.8863 
#####No differences between unhatched eggs by Species

#SIMPER to find potential indicators
dat.simp<-simper(dat.ra_UnhatchedOnly, metadata_UnhatchedOnly$Species, permutations = 9999)
sink("Simper_UnhatchedEggs.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

#SIMPER for Egg Discoloration
dat.simp<-simper(dat.ra_UnhatchedOnly, metadata_UnhatchedOnly$Egg.Coloration, permutations = 9999)
sink("Simper_EggDiscoloration.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

###Beach Comparison
#Remove CM Species
metadata_UnhatchedOnly_cc <- droplevels(metadata_UnhatchedOnly[!metadata_UnhatchedOnly$Species == 'CM',])
#Create a new object for common rownames for CC data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_UnhatchedOnly_cc))
#Set the data file and metadata file to have only the data that includes these common names 
dat_UnhatchedOnly_cc <- dat.001per[common.rownames,]
metadata_UnhatchedOnly_cc <- metadata_UnhatchedOnly_cc[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_UnhatchedOnly_cc), rownames(metadata_UnhatchedOnly_cc), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_UnhatchedOnly_cc<-decostand(dat_UnhatchedOnly_cc, method = "total")

##ANOSIM - Compare Beach
ano_Unhatched_beach = anosim(dat.ra_UnhatchedOnly_cc, metadata_UnhatchedOnly_cc$Beach, permutations = 9999, distance = "bray", strata = NULL)
ano_Unhatched_beach
#ANOSIM statistic R: 0.06747
#Significance: 0.032
#####Significant differences between unhatched eggs by Beach

#SIMPER to find potential indicators
dat.simp<-simper(dat.ra_UnhatchedOnly_cc, metadata_UnhatchedOnly_cc$Beach, permutations = 9999)
sink("Simper_UnhatchedEggsCC_Beach.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output



##########Potential Pathogen List
##Export HU Relative Abundance Tables transposed for filtering in Excel
dat.rat_cc_HU <- as.data.frame(t(dat.ra_cc_HU))
write.csv(dat.rat_cc_HU, "CC_HU.csv")
dat.rat_cm_HU <- as.data.frame(t(dat.ra_cm_HU))
write.csv(dat.rat_cm_HU, "CM_HU.csv")
dat.rat_F_HU <- as.data.frame(t(dat.ra_F_HU))
write.csv(dat.rat_F_HU, "F_HU.csv")
dat.rat_H_HU <- as.data.frame(t(dat.ra_H_HU))
write.csv(dat.rat_H_HU, "H_HU.csv")




###############Top GENUS Bar Charts
#Import Taxonomy Table using the taxonomy file from QIIME (edit in excel, separate KPODCFGS and remove confidence)
taxonomy = read.table(file= "Taxonomy.txt", header = TRUE, sep ="\t", row.names = 1)
#Load in QIIME tree
phy_tree = qza_to_phyloseq(tree="rooted-tree.qza")

###ALL Samples
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.UF = otu_table(as.matrix(dat.001per), taxa_are_rows=FALSE)
tax.UF = tax_table(as.matrix(taxonomy))
meta.UF = sample_data(metadata)
#Merge Files
physeq_ALL = phyloseq(ASV.UF,tax.UF,meta.UF,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_ALL_transform <- transform(physeq_ALL, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_ALL_transform), tax_table(physeq_ALL_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_ALL_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_All"


###All Samples - CC
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.CC = otu_table(as.matrix(dat_cc), taxa_are_rows=FALSE)
tax.CC = tax_table(as.matrix(taxonomy))
meta.CC = sample_data(metadata_cc)
#Merge Files
physeq_cc = phyloseq(ASV.CC,tax.CC,meta.CC,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_cc_transform <- transform(physeq_cc, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_cc_transform), tax_table(physeq_cc_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_cc_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_CC"



###All Samples - CM
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.CM = otu_table(as.matrix(dat_cm), taxa_are_rows=FALSE)
tax.CM = tax_table(as.matrix(taxonomy))
meta.CM = sample_data(metadata_cm)
#Merge Files
physeq_cm = phyloseq(ASV.CM,tax.CM,meta.CM,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_cm_transform <- transform(physeq_cm, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_cm_transform), tax_table(physeq_cm_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_cm_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_CM"



###All Samples - Fort Lauderdale
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.F = otu_table(as.matrix(dat_FT), taxa_are_rows=FALSE)
tax.F = tax_table(as.matrix(taxonomy))
meta.F = sample_data(metadata_FT)
#Merge Files
physeq_F = phyloseq(ASV.F,tax.F,meta.F,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_F_transform <- transform(physeq_F, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_F_transform), tax_table(physeq_F_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_F_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_F"


###ALL Samples - Hillsboro
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.H = otu_table(as.matrix(dat_H), taxa_are_rows=FALSE)
tax.H = tax_table(as.matrix(taxonomy))
meta.H = sample_data(metadata_H)
#Merge Files
physeq_H = phyloseq(ASV.H,tax.H,meta.H,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_H_transform <- transform(physeq_H, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_H_transform), tax_table(physeq_H_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_H_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_H"


###Hatched vs Unhatched Only (Figure 4)
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.HU = otu_table(as.matrix(dat_HU), taxa_are_rows=FALSE)
tax.HU = tax_table(as.matrix(taxonomy))
meta.HU = sample_data(metadata_HU)
#Merge Files
physeq_HU = phyloseq(ASV.HU,tax.HU,meta.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_HU_transform <- transform(physeq_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_HU_transform), tax_table(physeq_HU_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_HU_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_HU"


###Hatched vs Unhatched Only - CC
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.CC.HU = otu_table(as.matrix(dat_cc_HU), taxa_are_rows=FALSE)
tax.CC.HU = tax_table(as.matrix(taxonomy))
meta.CC.HU = sample_data(metadata_cc_HU)
#Merge Files
physeq_cc_HU = phyloseq(ASV.CC.HU,tax.CC.HU,meta.CC.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_cc_HU_transform <- transform(physeq_cc_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_cc_HU_transform), tax_table(physeq_cc_HU_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_cc_HU_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_HU_CC"



###Hatched vs Unhatched Only - CM
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.CM.HU = otu_table(as.matrix(dat_cm_HU), taxa_are_rows=FALSE)
tax.CM.HU = tax_table(as.matrix(taxonomy))
meta.CM.HU = sample_data(metadata_cm_HU)
#Merge Files
physeq_cm_HU = phyloseq(ASV.CM.HU,tax.CM.HU,meta.CM.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_cm_HU_transform <- transform(physeq_cm_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_cm_HU_transform), tax_table(physeq_cm_HU_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_cm_HU_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_HU_CM"



###Hatched vs Unhatched Only - Fort Lauderdale
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.F.HU = otu_table(as.matrix(dat_F_HU), taxa_are_rows=FALSE)
tax.F.HU = tax_table(as.matrix(taxonomy))
meta.F.HU = sample_data(metadata_F_HU)
#Merge Files
physeq_F_HU = phyloseq(ASV.F.HU,tax.F.HU,meta.F.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_F_HU_transform <- transform(physeq_F_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_F_HU_transform), tax_table(physeq_F_HU_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_F_HU_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_HU_F"



###Hatched vs Unhatched Only - Hillsboro
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.H.HU = otu_table(as.matrix(dat_H_HU), taxa_are_rows=FALSE)
tax.H.HU = tax_table(as.matrix(taxonomy))
meta.H.HU = sample_data(metadata_H_HU)
#Merge Files
physeq_H_HU = phyloseq(ASV.H.HU,tax.H.HU,meta.H.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_H_HU_transform <- transform(physeq_H_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_H_HU_transform), tax_table(physeq_H_HU_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_H_HU_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Sample.Type", fill="Genus") + 
  labs(title = "Genus Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_HU_H"



###Cloaca Only
#Rename species column in metadata (avoid confusion with Bacterial Species)
names(metadata_cloaca)[names(metadata_cloaca) == "Species"] <- "Turtle.Species"
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.cloaca = otu_table(as.matrix(dat_cloaca), taxa_are_rows=FALSE)
tax.cloaca = tax_table(as.matrix(taxonomy))
meta.cloaca = sample_data(metadata_cloaca)
#Merge Files
physeq_cloaca = phyloseq(ASV.cloaca,tax.cloaca,meta.cloaca,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows


#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_cloaca_transform <- transform(physeq_cloaca, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 5
top10g.names = sort(tapply(taxa_sums(physeq_cloaca_transform), tax_table(physeq_cloaca_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_cloaca_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Turtle.Species", fill="Genus") + 
  labs(title = "Cloaca Genus Composition by Species")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_Cloaca"


###Sand Only
##Control Sand
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.sand = otu_table(as.matrix(dat_sand), taxa_are_rows=FALSE)
tax.sand = tax_table(as.matrix(taxonomy))
meta.sand = sample_data(metadata_sand)
#Merge Files
physeq_sand = phyloseq(ASV.sand,tax.sand,meta.sand,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_sand_transform <- transform(physeq_sand, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 10
top10g.names = sort(tapply(taxa_sums(physeq_sand_transform), tax_table(physeq_sand_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_sand_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Beach", fill="Genus") + 
  labs(title = "Control Sand Genus Composition by Beach")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_ControlSand"

#Plot New Barchart based on Renourishment/R Zone
plot_bar(top10g, x="R.Zone", fill="Genus") + 
  labs(title = "Nest Sand Genus Composition by Beach")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_CS_RZone"


##Nest Sand
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.N = otu_table(as.matrix(dat_nest), taxa_are_rows=FALSE)
tax.N = tax_table(as.matrix(taxonomy))
meta.N = sample_data(metadata_N)
#Make Nest Number a factor
meta.N <- as.factor(meta.N$Nest.Number)
#Merge Files
physeq_N = phyloseq(ASV.N,tax.N,meta.N,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_N_transform <- transform(physeq_N, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 10
top10g.names = sort(tapply(taxa_sums(physeq_N_transform), tax_table(physeq_N_transform)[, "Genus"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
top10g = subset_taxa(physeq_N_transform, Genus %in% names(top10g.names))
#Plot New Barchart
plot_bar(top10g, x="Beach", fill="Genus") + 
  labs(title = "Nest Sand Genus Composition by Beach")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_NestSand"

#Plot New Barchart by Nest Number
plot_bar(top10g, x="Nest.Number", fill="Genus") + 
  labs(title = "Nest Sand Genus Composition by Beach")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "GenusBarChart_NS_NestNumber"





###########Top PHYLA Bar Charts (Figure 3)
#Import Taxonomy Table using the taxonomy file from QIIME (edit in excel, separate KPODCFGS and remove confidence)
taxonomy = read.table(file= "Taxonomy.txt", header = TRUE, sep ="\t", row.names = 1)
#Load in QIIME tree
phy_tree=qza_to_phyloseq(tree="rooted-tree.qza")

###Hatched vs Unhatched Only - CC
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.CC.HU = otu_table(as.matrix(dat_cc_HU), taxa_are_rows=FALSE)
tax.CC.HU = tax_table(as.matrix(taxonomy))
meta.CC.HU = sample_data(metadata_cc_HU)
#Merge Files
physeq_cc_HU = phyloseq(ASV.CC.HU,tax.CC.HU,meta.CC.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_cc_HU_transform <- transform(physeq_cc_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 20
phylum_bar.names = sort(tapply(taxa_sums(physeq_cc_HU_transform), tax_table(physeq_cc_HU_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
phylum_bar = subset_taxa(physeq_cc_HU_transform, Phylum %in% names(phylum_bar.names))
#Plot New Barchart
plot_bar(phylum_bar, x="Sample.Type", fill="Phylum") + 
  labs(title = "Phylum Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "PhylumBarChart_HU_CC"


###Hatched vs Unhatched Only - CM
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.CM.HU = otu_table(as.matrix(dat_cm_HU), taxa_are_rows=FALSE)
tax.CM.HU = tax_table(as.matrix(taxonomy))
meta.CM.HU = sample_data(metadata_cm_HU)
#Merge Files
physeq_CM_HU = phyloseq(ASV.CM.HU,tax.CM.HU,meta.CM.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_CM_HU_transform <- transform(physeq_CM_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 20
phylum_bar.names = sort(tapply(taxa_sums(physeq_CM_HU_transform), tax_table(physeq_CM_HU_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
phylum_bar = subset_taxa(physeq_CM_HU_transform, Phylum %in% names(phylum_bar.names))
#Plot New Barchart
plot_bar(phylum_bar, x="Sample.Type", fill="Phylum") + 
  labs(title = "Phylum Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "PhylumBarChart_HU_CM"



###Hatched vs Unhatched Only - Fort Lauderdale
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.F.HU = otu_table(as.matrix(dat_F_HU), taxa_are_rows=FALSE)
tax.F.HU = tax_table(as.matrix(taxonomy))
meta.F.HU = sample_data(metadata_F_HU)
#Merge Files
physeq_F_HU = phyloseq(ASV.F.HU,tax.F.HU,meta.F.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_F_HU_transform <- transform(physeq_F_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 20
phylum_bar.names = sort(tapply(taxa_sums(physeq_F_HU_transform), tax_table(physeq_F_HU_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
phylum_bar = subset_taxa(physeq_F_HU_transform, Phylum %in% names(phylum_bar.names))
#Plot New Barchart
plot_bar(phylum_bar, x="Sample.Type", fill="Phylum") + 
  labs(title = "Phylum Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "PhylumBarChart_HU_F"



###Hatched vs Unhatched Only - Hillsboro
#Make a phyloseq objects which includes the dat, metadata, and taxonomy
ASV.H.HU = otu_table(as.matrix(dat_H_HU), taxa_are_rows=FALSE)
tax.H.HU = tax_table(as.matrix(taxonomy))
meta.H.HU = sample_data(metadata_H_HU)
#Merge Files
physeq_H_HU = phyloseq(ASV.H.HU,tax.H.HU,meta.H.HU,phy_tree)
#If merge failed, try changing ASV.UF taxa_as_rows

#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_H_HU_transform <- transform(physeq_H_HU, "compositional")
##Top 10 Genus Barchart
#Sort the Genus by abundance and pick the top 20
phylum_bar.names = sort(tapply(taxa_sums(physeq_H_HU_transform), tax_table(physeq_H_HU_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Classes
phylum_bar = subset_taxa(physeq_H_HU_transform, Phylum %in% names(phylum_bar.names))
#Plot New Barchart
plot_bar(phylum_bar, x="Sample.Type", fill="Phylum") + 
  labs(title = "Phylum Composition by Sample Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Saved as "PhylumBarChart_HU_H"






####SIMPER  - Phylum Level
#Make file from relative abundance data consolidated by phylum
Phylum_Data <- read.delim("~/Grad School/Thesis/Data/Phylum_Data.txt", row.names=1)
####Convert Raw Data to Relative Abundance
dat.ra_phylum<-decostand(Phylum_Data, method = "total")

###########Make Files Separated by Turtle Species and Beach
#Create a new object for common rownames for CC data
common.rownames <- intersect(rownames(dat.ra_phylum),rownames(metadata_cc_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cc_HU <- dat.ra_phylum[common.rownames,]
metadata_cc_HU <- metadata_cc_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cc_HU), rownames(metadata_cc_HU), ignore.row.order = TRUE)

#Create a new object for common rownames for CM data
common.rownames <- intersect(rownames(dat.ra_phylum),rownames(metadata_cm_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cm_HU <- dat.ra_phylum[common.rownames,]
metadata_cm_HU <- metadata_cm_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cm_HU), rownames(metadata_cm_HU), ignore.row.order = TRUE)

#Create a new object for common rownames for FT data
common.rownames <- intersect(rownames(dat.ra_phylum),rownames(metadata_F_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_F_HU <- dat.ra_phylum[common.rownames,]
metadata_F_HU <- metadata_F_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_F_HU), rownames(metadata_F_HU), ignore.row.order = TRUE)

#Create a new object for common rownames for Hillsboro data
common.rownames <- intersect(rownames(dat.ra_phylum),rownames(metadata_H_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_H_HU <- dat.ra_phylum[common.rownames,]
metadata_H_HU <- metadata_H_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_H_HU), rownames(metadata_H_HU), ignore.row.order = TRUE)


#SIMPER
###Open dat.simp to get overall dissimilarity value
dat.simp<-simper(dat_cc_HU, metadata_cc_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_CC_phylum.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

dat.simp<-simper(dat_cm_HU, metadata_cm_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_CM_phylum.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

dat.simp<-simper(dat_F_HU, metadata_F_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_F_phylum.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

dat.simp<-simper(dat_H_HU, metadata_H_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_H_phylum.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output



####SIMPER - Genus Level
#Make file with consolidated by genus
Genus_Data <- read.delim("~/Grad School/Thesis/Data/Genus_Data.txt", row.names=1)
#Flip columns to rows
t.dat <- as.data.frame(t(Genus_Data))
#Set t.dat as dat
Genus_Data <- t.dat
#Convert to relative abundance
dat.ra_genus<-decostand(Genus_Data, method = "total")

###########Make Files Separated by Turtle Species and Beach
#Create a new object for common rownames for CC data
common.rownames <- intersect(rownames(dat.ra_genus),rownames(metadata_cc_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cc_HU <- dat.ra_genus[common.rownames,]
metadata_cc_HU <- metadata_cc_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cc_HU), rownames(metadata_cc_HU), ignore.row.order = TRUE)

#Create a new object for common rownames for CM data
common.rownames <- intersect(rownames(dat.ra_genus),rownames(metadata_cm_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_cm_HU <- dat.ra_genus[common.rownames,]
metadata_cm_HU <- metadata_cm_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_cm_HU), rownames(metadata_cm_HU), ignore.row.order = TRUE)

#Create a new object for common rownames for FT data
common.rownames <- intersect(rownames(dat.ra_genus),rownames(metadata_F_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_F_HU <- dat.ra_genus[common.rownames,]
metadata_F_HU <- metadata_F_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_F_HU), rownames(metadata_F_HU), ignore.row.order = TRUE)

#Create a new object for common rownames for Hillsboro data
common.rownames <- intersect(rownames(dat.ra_genus),rownames(metadata_H_HU))
#Set the data file and metadata file to have only the data that includes these common names 
dat_H_HU <- dat.ra_genus[common.rownames,]
metadata_H_HU <- metadata_H_HU[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_H_HU), rownames(metadata_H_HU), ignore.row.order = TRUE)


###SIMPER
###Open dat.simp to get overall dissimilarity value
dat.simp<-simper(dat_cc_HU, metadata_cc_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_CC_genus.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

dat.simp<-simper(dat_cm_HU, metadata_cm_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_CM_genus.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

dat.simp<-simper(dat_F_HU, metadata_F_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_F_genus.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output

dat.simp<-simper(dat_H_HU, metadata_H_HU$Sample.Type, permutations = 9999)
sink("Simper_HU_H_genus.csv")
#Creates csv of data in folder
summary(dat.simp)
sink()
#Closes output






##########Venn Diagrams (Figure 5)
##Export Transposed Count Tables for Binary Editing in Excel
dat.t_cc <- as.data.frame(t(dat_cc))
write.csv(dat.t_cc, "CC_Binary.csv")
dat.t_cm <- as.data.frame(t(dat_cm))
write.csv(dat.t_cm, "CM_Binary.csv")
dat.t_FT <- as.data.frame(t(dat_FT))
write.csv(dat.t_FT, "F_Binary.csv")
dat.t_H <- as.data.frame(t(dat_H))
write.csv(dat.t_H, "H_Binary.csv")

##ALL DATA
Venn_Data <- read.delim("~/Grad School/Thesis/Data/SampleType_ALL_Binary.txt", row.names=1)
venn(Venn_Data, sncs=0.5, col=c("black", "red", "blue", "gold2", "darkorange2"), lwd=3)
venn(Venn_Data, sncs=0.5)

##CC DATA
Venn_Data_CC <- read.delim("~/Grad School/Thesis/Data/CC_Binary.txt", row.names=1)
venn(Venn_Data_CC, sncs=0.5, col=c("black", "red", "blue", "gold2", "darkorange2"), lwd=3)
venn(Venn_Data_CC, sncs=0.5)

##CM DATA
Venn_Data_CM <- read.delim("~/Grad School/Thesis/Data/CM_Binary.txt", row.names=1)
venn(Venn_Data_CM, sncs=0.5, col=c("black", "red", "blue", "gold2", "darkorange2"), lwd=3)
venn(Venn_Data_CM, sncs=0.5)

##Fort Lauderdale DATA
Venn_Data_F <- read.delim("~/Grad School/Thesis/Data/F_Binary.txt", row.names=1)
venn(Venn_Data_F, sncs=0.5, col=c("black", "red", "blue", "gold2", "darkorange2"), lwd=3)
venn(Venn_Data_F, sncs=0.5)

##Hillsboro DATA
Venn_Data_H <- read.delim("~/Grad School/Thesis/Data/H_Binary.txt", row.names=1)
venn(Venn_Data_H, sncs=0.5, col=c("black", "red", "blue", "gold2", "darkorange2"), lwd=3)
venn(Venn_Data_H, sncs=0.5)



###Recaptured Turtle
#Remove Hatched Eggs
metadata_RECAP <- droplevels(metadata[!metadata$Nest.Number == '127',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '236',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '242',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '282',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '323',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '324',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '395',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '392',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '463',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '695',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '696',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '1166',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '1167',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '1168',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '457',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '456',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '519',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '520',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '521',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '522',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '697',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '698',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '699',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '823',])
metadata_RECAP <- droplevels(metadata_RECAP[!metadata_RECAP$Nest.Number == '892',])

##Analyze just the hatched eggs
metadata_RECAP_hatched <- droplevels(metadata_RECAP[!metadata_RECAP$Sample.Type == 'Cloaca',])
metadata_RECAP_hatched <- droplevels(metadata_RECAP_hatched[!metadata_RECAP_hatched$Sample.Type == 'Nest Sand',])
metadata_RECAP_hatched <- droplevels(metadata_RECAP_hatched[!metadata_RECAP_hatched$Sample.Type == 'Control Sand',])
metadata_RECAP_hatched <- droplevels(metadata_RECAP_hatched[!metadata_RECAP_hatched$Sample.Type == 'Unhatched Egg',])

#Create a new object for common rownames for data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_RECAP_hatched))
#Set the data file and metadata file to have only the data that includes these common names 
dat_RECAP_hatched <- dat.001per[common.rownames,]
metadata_RECAP_hatched <- metadata_RECAP_hatched[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_RECAP_hatched), rownames(metadata_RECAP_hatched), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_RECAP_hatched<-decostand(dat_RECAP_hatched, method = "total")

#ANOSIM - Compare Beaches
ano_RECAP_hatched = anosim(dat.ra_RECAP_hatched, metadata_RECAP_hatched$Beach, permutations = 9999, distance = "bray", strata = NULL)
ano_RECAP_hatched
#ANOSIM statistic R: 0.9259  
#Significance: 0.1  

##Analyze just the Unhatched eggs
metadata_RECAP_unhatched <- droplevels(metadata_RECAP[!metadata_RECAP$Sample.Type == 'Cloaca',])
metadata_RECAP_unhatched <- droplevels(metadata_RECAP_unhatched[!metadata_RECAP_unhatched$Sample.Type == 'Nest Sand',])
metadata_RECAP_unhatched <- droplevels(metadata_RECAP_unhatched[!metadata_RECAP_unhatched$Sample.Type == 'Control Sand',])
metadata_RECAP_unhatched <- droplevels(metadata_RECAP_unhatched[!metadata_RECAP_unhatched$Sample.Type == 'Hatched Egg',])

#Create a new object for common rownames for data
common.rownames <- intersect(rownames(dat.001per),rownames(metadata_RECAP_unhatched))
#Set the data file and metadata file to have only the data that includes these common names 
dat_RECAP_unhatched <- dat.001per[common.rownames,]
metadata_RECAP_unhatched <- metadata_RECAP_unhatched[common.rownames,]
#Make sure all the row names are the same (equal) following the code
all.equal(rownames(dat_RECAP_unhatched), rownames(metadata_RECAP_unhatched), ignore.row.order = TRUE)
#Convert to relative abundance
dat.ra_RECAP_unhatched<-decostand(dat_RECAP_unhatched, method = "total")

#ANOSIM - Compare Beaches
ano_RECAP_unhatched = anosim(dat.ra_RECAP_unhatched, metadata_RECAP_unhatched$Beach, permutations = 9999, distance = "bray", strata = NULL)
ano_RECAP_unhatched
#ANOSIM statistic R: 0.7778 
#Significance: 0.1  




###########SourceTracker (Figure 6)
metadata <- read.delim("~/Grad School/Thesis/Data/Turtle_Metadata_2021_SourceTracker.txt", row.names=1)

#Extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(dat.001per))
otus <- dat.001per[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
#Double-check that the mapping file and otu table had overlapping samples
all.equal(rownames(otus),rownames(metadata))

#Extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env

#Download latest SourceTracker release from: https://github.com/danknights/sourcetracker/releases
#Load SourceTracker package
setwd("C:/Users/colmc/Downloads/sourcetracker-1.0.1/sourcetracker-1.0.1")
source('src/SourceTracker.r')

#Note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

#Train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])

#Estimate source proportions in test data
Results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)
#This can be exported as a CSV to make averaged pie charts

#Plot results
plot(Results, type='pie')



###########SourceTracker: cloaca, control sand, nest sand
metadata <- read.delim("~/Grad School/Thesis/Data/Turtle_Metadata_2021_SourceTracker2.txt", row.names=1)
#Extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(dat.001per))
otus <- dat.001per[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
#Double-check that the mapping file and otu table had overlapping samples
all.equal(rownames(otus),rownames(metadata))

#Extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env
#Download latest SourceTracker release from: https://github.com/danknights/sourcetracker/releases
#Load SourceTracker package
setwd("~/Grad School/Thesis/Data/sourcetracker-1.0.1")
source('src/SourceTracker.r')
#Note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001
#Train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])
#Estimate source proportions in test data
Results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)
#This can be exported as a CSV to make averaged pie charts






###########F/B Ratio
#Export previously made phylum relative abundance table 
#Remove all phyla except Firmicutes and Bacteriodetes
#Remove Nest 395.CC.F from analysis since there were no unhatched eggs to compare
write.csv(dat.ra_phylum, "Phylum_RAData.csv")

###Paired t-test (hatched v unhatched) by Nest Number
FB_Ratio_Nest_Pairedt <- read.delim("~/Grad School/Thesis/Data/FB_Ratio_Nest_Pairedt.txt", row.names=1)
attach(FB_Ratio_Nest_Pairedt)

#Test for normality
shapiro.test(Hatched.Egg-Unhatched.Egg)
#W = 0.7902, p-value = 0.0001205
###Data is not normally distributed >> use non-parametric test

#Non-parametric test: Two-tailed Wilcoxon test
wilcox.test (Hatched.Egg, Unhatched.Egg, paired=T)
#V = 70, p-value = 0.006131

#Non-parametric test: One-tailed Wilcoxon test
#Test if unhatched eggs have greater ratio than hatched eggs
wilcox.test (Hatched.Egg, Unhatched.Egg, paired=T,alternative='less')
#V = 70, p-value = 0.003065

#Check Reverse: hatched eggs have greater ratio than unhatched eggs
wilcox.test (Hatched.Egg, Unhatched.Egg, paired=T,alternative='greater')
#V = 70, p-value = 0.9972

detach(FB_Ratio_Nest_Pairedt)



