# Microbiome Analysis in R
# WSB Probiotic Microbiome

setwd("~/Probiotics/WSB 2024/data")

# view metadata as needed
metadata <- read.delim("~/Probiotics/WSB 2024/data/metadata.tsv")
View(metadata)

# install and load qiime2R
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

library(qiime2R)

#Import qiime2 objects as phyloseq object

# phyloseq
physeq <- qza_to_phyloseq(features="asv-table-taxa-sample-filt.qza",
                          tree="tree.qza",
                          taxonomy="gg-taxonomy.qza",
                          metadata= "metadata.tsv")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)

# Rarefy data
hist(sample_sums(physeq), main="Distribution of Sample Read Counts", xlab="Total Reads per Sample")

set.seed(6500) # For reproducibility
rarefied_physeq <- rarefy_even_depth(physeq, sample.size = 22180, replace = FALSE, trimOTUs = TRUE)


### Alpha Diversity ----------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")

install.packages("ggpubr")

library(microbiome)
library(ggpubr)
library(ggplot2)

# LCW_LF @ 11 and 18 dph

LCW_LF<-subset_samples(rarefied_physeq, sample_type != "ART")
LCW_LF_11_18<-subset_samples(LCW_LF, days_post_hatch %in% c("11", "18"))
sample_data(LCW_LF_11_18)

plot_richness(LCW_LF_11_18, color = "sample_type", x = "sample_type", 
              measures = c("Observed", "Chao1", "Shannon", "Simpson"))
# make table with alpha diversity metrics
LCW_LF_11_18_alpha_diversity<- estimate_richness(LCW_LF_11_18)
LCW_LF_11_18_alpha_diversity<- cbind(sample_data(LCW_LF_11_18), LCW_LF_11_18_alpha_diversity)
LCW_LF_11_18_evenness<- evenness(LCW_LF_11_18, index = "pielou", zeroes=TRUE, detection=0)
LCW_LF_11_18_alpha_diversity<- cbind(LCW_LF_11_18_alpha_diversity, LCW_LF_11_18_evenness)

# plot
a<-ggplot(LCW_LF_11_18_alpha_diversity, aes(x=sample_type, y=pielou))+
  geom_boxplot()
b<-ggplot(LCW_LF_11_18_alpha_diversity, aes(x=sample_type, y=Observed))+
  geom_boxplot()
c<-ggplot(LCW_LF_11_18_alpha_diversity, aes(x=sample_type, y=Shannon))+
  geom_boxplot()
ggarrange(a, b, c,
          ncol = 3, nrow = 1)
a<-ggplot(LCW_LF_11_18_alpha_diversity, aes(x=sample_type, y=pielou))+
  geom_boxplot()+
  facet_wrap(~days_post_hatch)
b<-ggplot(LCW_LF_11_18_alpha_diversity, aes(x=sample_type, y=Observed))+
  geom_boxplot()+
  facet_wrap(~days_post_hatch)
c<-ggplot(LCW_LF_11_18_alpha_diversity, aes(x=sample_type, y=Shannon))+
  geom_boxplot()+
  facet_wrap(~days_post_hatch)
ggarrange(a, b, c,
          ncol = 3, nrow = 1)

# Statistics
kruskal.test(Observed ~ sample_type, LCW_LF_11_18_alpha_diversity)
kruskal.test(pielou ~ sample_type, LCW_LF_11_18_alpha_diversity)
kruskal.test(Shannon~ sample_type, LCW_LF_11_18_alpha_diversity)


# LCW
LCW<-subset_samples(rarefied_physeq, sample_type == "LCW")
sample_data(LCW)
# make table with alpha diversity metrics
LCW_alpha_diversity<- estimate_richness(LCW)
LCW_alpha_diversity<- cbind(sample_data(LCW), LCW_alpha_diversity)
LCW_evenness<- evenness(LCW, index = "pielou", zeroes=TRUE, detection=0)
LCW_alpha_diversity<- cbind(LCW_alpha_diversity, LCW_evenness)

# plot
a<-ggplot(LCW_alpha_diversity, aes(x=treatment, y=pielou))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~days_post_hatch)
b<-ggplot(LCW_alpha_diversity, aes(x=treatment, y=Observed))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~days_post_hatch)
c<-ggplot(LCW_alpha_diversity, aes(x=treatment, y=Shannon))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~days_post_hatch)
ggarrange(a, b, c,
          ncol = 3, nrow = 1)

d<-ggplot(LCW_alpha_diversity, aes(x=treatment, y=pielou))+
  geom_boxplot()
e<-ggplot(LCW_alpha_diversity, aes(x=treatment, y=Observed))+
  geom_boxplot()
f<-ggplot(LCW_alpha_diversity, aes(x=treatment, y=Shannon))+
  geom_boxplot()
ggarrange(d, e, f,
          ncol = 1, nrow = 3)

#Statistics
kruskal.test(Observed ~ treatment, LCW_alpha_diversity)
kruskal.test(pielou ~ treatment, LCW_alpha_diversity)
kruskal.test(Shannon ~ treatment, LCW_alpha_diversity)

# LF
LF<-subset_samples(rarefied_physeq, sample_type == "LF")
sample_data(LF)
# make table with alpha diversity metrics
LF_alpha_diversity<- estimate_richness(LF)
LF_alpha_diversity<- cbind(sample_data(LF), LF_alpha_diversity)
LF_evenness<- evenness(LF, index = "pielou", zeroes=TRUE, detection=0)
LF_alpha_diversity<- cbind(LF_alpha_diversity, LF_evenness)

# mean values per days post hatch
library(dplyr)
LF_alpha_diversity %>%
  group_by(days_post_hatch) %>%
  summarise(mean_value = mean(Observed))

LF_alpha_diversity %>%
  group_by(days_post_hatch) %>%
  summarise(mean_value = mean(Shannon))

# plot
a<-ggplot(LF_alpha_diversity, aes(x=treatment, y=pielou))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~days_post_hatch)
b<-ggplot(LF_alpha_diversity, aes(x=treatment, y=Observed))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~days_post_hatch)
c<-ggplot(LF_alpha_diversity, aes(x=treatment, y=Shannon))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~days_post_hatch)
ggarrange(a, b, c,
          ncol = 3, nrow = 1)

d<-ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=pielou))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~treatment)
e<-ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=Observed))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~treatment)
f<-ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=Shannon))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~treatment)
ggarrange(d, e, f,
          ncol = 3, nrow = 1)

g<-ggplot(LF_alpha_diversity, aes(x=treatment, y=pielou))+
  geom_boxplot()
h<-ggplot(LF_alpha_diversity, aes(x=treatment, y=Observed))+
  geom_boxplot()
i<-ggplot(LF_alpha_diversity, aes(x=treatment, y=Shannon))+
  geom_boxplot()
ggarrange(g, h, i,
          ncol = 1, nrow = 3)

j<-ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=pielou))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)
k<-ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=Observed))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)
l<-ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=Shannon))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)
ggarrange(j, k, l,
          ncol = 3, nrow = 1)

# plot all treatments over time
ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=Shannon, color=treatment)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=1) +
  stat_summary(geom="line", fun.data=mean_se, linewidth=1) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days Post Hatch") +
  ylab("Shannon Diversity") +
  theme_bw()+
  scale_color_viridis_d(name="Treatment")
ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=pielou, color=treatment)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=1) +
  stat_summary(geom="line", fun.data=mean_se, linewidth=1) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days Post Hatch") +
  ylab("Evenness") +
  theme_bw()+
  scale_color_viridis_d(name="Treatment")
ggplot(LF_alpha_diversity, aes(x=days_post_hatch, y=Observed, color=treatment)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=1) +
  stat_summary(geom="line", fun.data=mean_se, linewidth=1) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days Post Hatch") +
  ylab("Richness") +
  theme_bw()+
  scale_color_viridis_d(name="Treatment")

#Statistics
kruskal.test(Observed ~ treatment, LF_alpha_diversity)
kruskal.test(pielou ~ treatment, LF_alpha_diversity)
kruskal.test(Shannon ~ treatment, LF_alpha_diversity)

library(rstatix)
kruskal.test(Observed ~ days_post_hatch, LF_alpha_diversity)
dunn_test(LF_alpha_diversity, Observed ~ days_post_hatch, p.adjust.method = "BH")
kruskal.test(pielou ~ days_post_hatch, LF_alpha_diversity)
dunn_test(LF_alpha_diversity, pielou ~ days_post_hatch, p.adjust.method = "BH")
kruskal.test(Shannon ~ days_post_hatch, LF_alpha_diversity)
dunn_test(LF_alpha_diversity, Shannon ~ days_post_hatch, p.adjust.method = "BH")

# Linear model for days post hatch
# a linear model is used because time is treated as a fixed effect and diversity as response
# https://libguides.princeton.edu/R-linear_regression

dph_lm_pielou<-lm(pielou ~ days_post_hatch, data = LF_alpha_diversity)
summary(dph_lm_pielou)
# robust regression
library(lmtest)
library(sandwich)
dph_lm_pielou$robse <- vcovHC(dph_lm_pielou, type="HC1")
coeftest(dph_lm_pielou,dph_lm_pielou$robse)
#normality check
shapiro.test(LF_alpha_diversity$pielou) #normal
library(car)
qqPlot(dph_lm_pielou) #2 outliers
ncvTest(dph_lm_pielou) #heteroscedastic
residualPlots(dph_lm_pielou) #bad 
# test normality

dph_lm_shannon<-lm(Shannon ~ days_post_hatch, data = LF_alpha_diversity)
summary(dph_lm_shannon)
dph_lm_shannon$robse <- vcovHC(dph_lm_shannon, type="HC1")
coeftest(dph_lm_shannon, dph_lm_pielou$robse)
shapiro.test(LF_alpha_diversity$Shannon) #normal
qqPlot(dph_lm_shannon) #2 outliers
ncvTest(dph_lm_shannon) #heteroscedastic
residualPlots(dph_lm_shannon) # curved - not good


# test for homogeneity of variance
bartlett.test(pielou ~ days_post_hatch, data = LF_alpha_diversity) #not equal variance
bartlett.test(Shannon ~ days_post_hatch, data = LF_alpha_diversity) #not equal variance

cor.test(x=LF_alpha_diversity$days_post_hatch, y=LF_alpha_diversity$Observed, method = 'spearman')
cor.test(x=LF_alpha_diversity$days_post_hatch, y=LF_alpha_diversity$pielou, method = 'spearman')
cor.test(x=LF_alpha_diversity$days_post_hatch, y=LF_alpha_diversity$Shannon, method = 'spearman')


# Taxa bar plots --------
# in order to merge samples within treatments, better to have a numerical column for treatment
# make this is metadata3
physeq <- qza_to_phyloseq(features="asv-table-taxa-sample-filt.qza",
                          tree="tree.qza",
                          taxonomy="gg-taxonomy.qza",
                          metadata= "metadata3.tsv")
set.seed(6500) # For reproducibility
rarefied_physeq <- rarefy_even_depth(physeq, sample.size = 22180, replace = FALSE, trimOTUs = TRUE)

# first, need to merge samples across treatment groups
# https://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html
rarefied_physeq_merged = merge_samples(rarefied_physeq, "treatment")

# repair merged values
# now using treatment_numeric
# order is 2, 1, 3, 4, 5
sample_data(rarefied_physeq_merged)$treatment <- c("Control", "Probiotic Artemia", "Probiotic fed via Artemia", "Probiotic in water", "Rotifer addition")

# transform to %
rarefied_physeq_merged = transform_sample_counts(rarefied_physeq_merged, function(x) 100 * x/sum(x))

# select only top 20 most abundant OTUs
top20otus = names(sort(taxa_sums(rarefied_physeq_merged), TRUE)[1:20])
taxtab20 = cbind(tax_table(rarefied_physeq_merged), family20 = NA)
taxtab20[top20otus, "family20"] <- as(tax_table(rarefied_physeq_merged)[top20otus, "Family"], 
                                      "character")
tax_table(rarefied_physeq_merged) <- tax_table(taxtab20)

#create abundance plot
title = "Taxa Barplot Across treatments"
plot_bar(rarefied_physeq_merged, "treatment", fill = "family20", title = title) + coord_flip()

###### create barplots per sample type ##########
# dont forget to use metadata3 phyloseq

# LCW 
LCW<-subset_samples(rarefied_physeq, sample_type == "LCW")
# merge LCW samples across treatment group
LCW_merged = merge_samples(LCW, "treatment")
# repair merged phyloseq
sample_data(LCW_merged)$treatment <- c("Control", "Probiotic fed via Artemia", 
                                       "Probiotic in water", "Rotifer addition")
# transform to %
LCW_merged = transform_sample_counts(LCW_merged, function(x) 100 * x/sum(x))
# select top 20
top20otus = names(sort(taxa_sums(LCW_merged), TRUE)[1:20])
taxtab20 = cbind(tax_table(LCW_merged), family20 = NA)
taxtab20[top20otus, "family20"] <- as(tax_table(LCW_merged)[top20otus, "Family"], 
                                      "character")
tax_table(LCW_merged) <- tax_table(taxtab20)
# plot
title = "Taxa Barplot: LCW Across treatments"
plot_bar(LCW_merged, "treatment", fill = "family20", title = title) + coord_flip()

# agglomerate at genus instead
# select top 20
taxtabgenus20 = cbind(tax_table(LCW_merged), genus20 = NA)
taxtabgenus20[top20otus, "genus20"] <- as(tax_table(LCW_merged)[top20otus, "Genus"], 
                                      "character")
tax_table(LCW_merged) <- tax_table(taxtabgenus20)
# plot
title = "Taxa Barplot: LCW Across treatments"
plot_bar(LCW_merged, "treatment", fill = "genus20", title = title) + coord_flip()

# LF
LF<-subset_samples(rarefied_physeq, sample_type == "LF")
# merge LCW samples across treatment group
LF_merged = merge_samples(LF, "treatment")
# repair merged phyloseq
sample_data(LF_merged)$treatment <- c("Control", "Probiotic fed via Artemia", 
                                       "Probiotic in water", "Rotifer addition")
# transform to %
LF_merged = transform_sample_counts(LF_merged, function(x) 100 * x/sum(x))
# select top 20
top20otus = names(sort(taxa_sums(LF_merged), TRUE)[1:20])
taxtab20 = cbind(tax_table(LF_merged), family20 = NA)
taxtab20[top20otus, "family20"] <- as(tax_table(LF_merged)[top20otus, "Family"], 
                                      "character")
tax_table(LF_merged) <- tax_table(taxtab20)
# plot
title = "Taxa Barplot: LF Across treatments"
plot_bar(LF_merged, "treatment", fill = "family20", title = title) + coord_flip()

# agglomerate at genus instead
# select top 20
taxtabgenus20 = cbind(tax_table(LF_merged), genus20 = NA)
taxtabgenus20[top20otus, "genus20"] <- as(tax_table(LF_merged)[top20otus, "Genus"], 
                                          "character")
tax_table(LF_merged) <- tax_table(taxtabgenus20)
# plot
title = "Taxa Barplot: LF Across treatments"
plot_bar(LF_merged, "treatment", fill = "genus20", title = title) + coord_flip()



# Differential Abundance of LF over time -----

# first need to generate tse from data
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mia")

library(mia)


# Venn Diagram -----

# Option 1
# https://microbiome.github.io/tutorials/core_venn.html 
install.packages("eulerr")
install.packages("microbiome")
devtools::install_github('microsud/microbiomeutilities')
library(eulerr)
library(microbiome)
library(microbiomeutilities)

# how many samples per sample type?
table(meta(rarefied_physeq)$sample_type, useNA = "always")

# convert to relative abundances
pseq.rel <- microbiome::transform(rarefied_physeq, "compositional")

# make a list
sample_types <- unique(as.character(meta(pseq.rel)$sample_type))
print(sample_types)

# Write a for loop to go through each of the sample types one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in sample_types){ # for each variable n in sample_type
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, sample_type == n) # Choose sample from sample type by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Sample type.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

# Specify colors and plot venn
# supplying colors in the order they appear in list_core
mycols <- c(ART="#d6e2e9", LCW="#cbf3f0", LF="#fcf5c7")
plot(venn(list_core),
     fills = mycols)

print(list_core)

# Option 2
#https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html#5-beta-analysis
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("MicrobiotaProcess")

# subset to only 11 and 18 dph

lcw_lf_art_11_18<-subset_samples(rarefied_physeq, days_post_hatch %in% c("11", "18"))
View(sample_data(lcw_lf_art_11_18))

# subset to only probiotic conditions

lcw_lf_art_11_18_pro<-subset_samples(lcw_lf_art_11_18, treatment %in% c("MIC", "P1", "P2"))
View(sample_data(lcw_lf_art_11_18_pro))

# create venn

library(MicrobiotaProcess)
vennlist <- get_vennlist(obj=lcw_lf_art_11_18_pro, factorNames="sample_type")
library(VennDiagram)
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("purple", "lightblue", "lightgreen"),
                      cat.col=c("purple", "lightblue", "lightgreen"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
grid::grid.draw(vennp)
# NOTE that ART has 430 ASVs, LF has 1007, and LCW has 1785

ART=126/430 # 29% of ART ASVs are unique to ART
LF=328/1007 # 33% of LF ASVs are unique to LF
LCW=1101/1785 # 62% of LCW ASVs are unique to LCW

# percentage of contribution from ART to LF
54/1007 # 5%
# percentage of contribution from LCW to LF
434/1007 # 43%
# number of ASVs shared by ART and LCW but not taken up by LF
59

# can also export venn list and use in Jvenn
capture.output(vennlist[["ART"]], file = "ART_taxa.txt")

# make data frames for jvenn
ART_taxa<-as.data.frame(vennlist[["ART"]])
ART_taxa$sample_type <- "ART"
colnames(ART_taxa)[1] <- "ASV"

LCW_taxa<-as.data.frame(vennlist[["LCW"]])
LCW_taxa$sample_type <- "LCW"
colnames(LCW_taxa)[1] <- "ASV"

LF_taxa<-as.data.frame(vennlist[["LF"]])
LF_taxa$sample_type <- "LF"
colnames(LF_taxa)[1] <- "ASV"

combined_taxa<-rbind(ART_taxa, LCW_taxa, LF_taxa)
#export as .csv file for Jvenn
write.csv(combined_taxa, "combined_taxa.csv", row.names = FALSE)


# Venn Diagram for 11 dph comparing LCW and LF_____
LCW_LF<-subset_samples(rarefied_physeq, sample_type != "ART")
LCW_LF_11<-subset_samples(LCW_LF, days_post_hatch == "11")

# create venn
library(MicrobiotaProcess)
vennlist <- get_vennlist(obj=LCW_LF_11, factorNames="sample_type")
library(VennDiagram)
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("orange", "blue"),
                      cat.col=c("orange", "blue"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
grid::grid.draw(vennp)
# 58.37% LF taxa shared with LCW

# Venn Diagram for 18 dph comparing LCW and LF_____
LCW_LF<-subset_samples(rarefied_physeq, sample_type != "ART")
LCW_LF_18<-subset_samples(LCW_LF, days_post_hatch == "18")

# create venn
library(MicrobiotaProcess)
vennlist <- get_vennlist(obj=LCW_LF_18, factorNames="sample_type")
library(VennDiagram)
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("orange", "blue"),
                      cat.col=c("orange", "blue"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
grid::grid.draw(vennp)
# 69.11% LF taxa shared with LCW



# SHARED TAXA between only ART and LF

# import list of 54 taxa from jvenn
common_elements_ART_LF <- read.csv("~/Probiotics/WSB 2024/data/common_elements_ART_LF.txt", 
                                   sep="")

# extract and manipulate taxonomy table
taxa_assign <- as.data.frame(rarefied_physeq@tax_table)

#rownames to column
taxa_assign$Common_elements_art_lf <- row.names(taxa_assign)

#assign taxonomy
taxa_assign_common_elements<-merge(common_elements_ART_LF, taxa_assign, by = "Common_elements_art_lf", 
      all.x = TRUE)




# remove names from metadata
# Define the IDs to remove from metadata 
ids_to_remove <- c("LCW.C.5d.3", "LCW.C.5d.6", "LCW.C.5d.9",
                   "LCW.P1.5d.12", "LCW.P1.5d.16", "LCW.P1.5d.2",
                   "LCW.P1.5d.4", "LCW.P2.5d.11", "LCW.P2.5d.15",
                   "LCW.P2.5d.5", "LCW.P2.5d.8", "LCW.R.5d.10",
                   "LCW.R.5d.7", "ART.c.11d.1", "ART.c.11d.2",
                   "ART.c.11d.3", "ART.c.18d.1", "ART.c.18d.2",
                   "ART.c.18d.3", "ART.c.5d.1", "ART.c.5d.2",
                   "ART.c.5d.3")

# Remove rows where 'ID' column matches any of the IDs in 'ids_to_remove'
subset_metadata <- metadata[!metadata$sample.id %in% ids_to_remove, ]


