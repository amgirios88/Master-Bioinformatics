# check R version for assessment purposes

R.version

# this assessment has been done under R version 4.4.1

# installing any package needed that has not been installed before in this computer

# install Rtools - we will use BiocManager to install DESeq2 because otherwise it returns an error because
# Rstudio is trying to install it and not finding it for this version (same for the rest)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

BiocManager::install("AnnotationDbi")

BiocManager::install("org.Hs.eg.db")

BiocManager::install("GO.db")

BiocManager::install("caret")

BiocManager::install("DESeq2")

BiocManager::install("car")

BiocManager::install("EnhancedVolcano")

# Loading all libraries needed beforehand (will add more as I keep adding them to 
# remove them from the middle of the coded)

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggfortify)
library(RColorBrewer)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(caret)
library(DESeq2)
library(car)
library(EnhancedVolcano)

# loading citations
citation('Biobase')
citation()
RStudio.Version()



# 1. Prepare the data --------------------------------------------------------

# all csv files are tab-separated (\t) with header as the first line

annot <- read.table("D:/Uni/MASTER/Classes/Stats/assessment/Annotations.csv", sep = "\t", header = TRUE)

# we need to know the structure of the data to be able to work with it

# check the data (first rows)
head(annot)
# check the structure of the data
str(annot)
# check the summary of the data
summary(annot)

# variable names are ID (chr), symbol(chr), biotype(chr), chromosome(chr), start(int), stop(int)

gout <- read.table("D:/Uni/MASTER/Classes/Stats/assessment/DE_GOUT_vs_HC.csv", sep = "\t", header = TRUE)

# check the data (first rows)
head(gout)
# check the structure of the data
str(gout)
# check the summary of the data
summary(gout)

# variable names are ID (chr), log2Fold (num), p (num), p.adj (num) - for gout vs HC

sa <- read.table("D:/Uni/MASTER/Classes/Stats/assessment/DE_SA_vs_HC.csv", sep = "\t", header = TRUE)

# check the data (first rows)
head(sa)
# check the structure of the data
str(sa)
# check the summary of the data
summary(sa)

# variable names are ID (chr), log2Fold (num), p (num), p.adj (num) - for SA vs HC

exp <- read.table("D:/Uni/MASTER/Classes/Stats/assessment/Expression_Table.csv", sep = "\t", header = TRUE)

# check the data (first rows)
head(exp)
# check the structure of the data
str(exp)
# check the summary of the data
summary(exp)

# we have here 28 variables, 1st is ID (chr) and the rest are numeric, which correspond to samples 1 to 9 for HC,
# 1 to 9 for gout, 1 to 9 for SA

sample <- read.table("D:/Uni/MASTER/Classes/Stats/assessment/Sample_Information.csv", sep = "\t", header = TRUE)

# check the data (first rows)
head(sample)
# check the structure of the data
str(sample)
# check the summary of the data
summary(sample)

# variable names are SAMPLE (chr), SAMPLE_GROUP (chr), SEX (chr), NEUTROPHILS (num)



# 2. Merge data with annotations --------------------------------------------------------

# we want to ensure we have gene biotype and other details in our dataframes for differential expression for 
# gout and SA

# we merge annotations dataframe with differential expression results for gout vs HC
gout_annot <- merge(gout, annot, by = "ID")

# we merge annotations dataframe with differential expression results for SA vs HC
sa_annot <- merge(sa, annot, by = "ID")

# we can now inspect our merged data; this is a good practice to check if we have merged the data correctly
head(gout_annot)
head(sa_annot)



# 3. Filter significant genes --------------------------------------------------------

# what are the most significant genes for each? to answer this question we will sort the dataframes by p.adj values

# we use the next code to subset from gout_annot to keep only values for which p.adj is less than 0.05
gout_sig <- subset(gout_annot, p.adj < 0.05)

# we repeat this for SA vs HC
sa_sig <- subset(sa_annot, p.adj < 0.05)

# in the results for gout vs HC, we have many NA values in p.adj column. After researching on this, I found that
# there are some reasons why this may happen, listed in the DEseq2 documentation in the bioconductor project website:
#       -If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change 
#       estimates, p value and adjusted p value will all be set to NA.
#       -If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set 
#       to NA. These outlier counts are detected by Cook’s distance.
#       -If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only 
#       the adjusted p value will be set to NA. (This is done to increase statistical power, as genes with extremely 
#       low expression across all samples generally provide little information and have a low chance of being 
#       statistically significant.)

# the method used for calculating the differential gene expression in the data provided was not mentioned, but
# DESeq2 is a popular package for this, so we can assume that the reasons provided may apply in general

# now that we have understood the reason behind the numerous NA values for gout data, we can continue



# 4. Exploring expression profile for each group ----------------------------------------

# to reply to this question, we will perform first a PCA to check if the groups are distinct based on their gene expression
# and then we will analyse the differentially expressed genes unique to gout and SA with DESeq2

# 4.1. PCA of gene expression for HC, Gout, and Sepsis ------------------------
#* @TASK4-Part1

# it is interesting to check if there are differences in gene expression profiles between the groups, because
# this will tell us if the groups are distinct based on their gene expression profiles, and then checking for
# significant differentially expressed genes unique to gout and sa that could be used for diagnostic purposes
# makes complete sense

# Are there any genes that are significantly different between Gout and SA? If so, what are they? 
# What do they look like when plotted?

#* *Exploring expression profile for each group*

# we need to transpose the expression dataframe and set gene names as column headers so we can merge this information
# with the clinical dataframe
exp_transposed <- as.data.frame(t(exp[, -1])) 
colnames(exp_transposed) <- exp$ID

# now we can merge transposed expression data with clinical data based on sample names
exp_clin <- merge(exp_transposed, sample, by.x = "row.names", by.y = "SAMPLE", all = TRUE)

# and we set appropriate column name for samples
colnames(exp_clin)[1] <- "SAMPLE"

# we can reorder columns to have clinical data first so the dataframe is easier to read
exp_clin <- exp_clin[, c("SAMPLE", "SAMPLE_GROUP", "SEX", "NEUTROPHILS", 
                         setdiff(colnames(exp_clin), c("SAMPLE", "SAMPLE_GROUP", "SEX", "NEUTROPHILS")))]

# we will prepare the data removing non-numeric columns (SAMPLE, SAMPLE_GROUP, SEX, NEUTROPHILS)
exp_pca <- exp_clin[, !(colnames(exp_clin) %in% c("SAMPLE", "SAMPLE_GROUP", "SEX", "NEUTROPHILS"))]

# standardize the data (important step before PCA to make the data comparable - PCA seeks to maximize the variance 
# of each component; a PCA without normalization would perform worse than one with normalization)
# note to myself about why is important to standardize the data before PCA: 
# https://stats.stackexchange.com/questions/69157/why-do-we-need-to-normalize-data-before-principal-component-analysis-pca
exp_pca_scaled <- scale(exp_pca)

# the PCA will be performed on the samples (rows) and we will see where these patients fit according to their gene
# expression profile
# we perform PCA on the standardized data
pca_result <- prcomp(exp_pca_scaled, center = TRUE, scale. = TRUE)

# this is the PCA summary to check variance explained by each component
summary(pca_result)

# to plot the data, we can extract PCA coordinates for the samples (rows)
pca_data <- as.data.frame(pca_result$x)

# we will add back the sample group information (GOUT, SEPSIS, HC)
pca_data$SAMPLE_GROUP <- exp_clin$SAMPLE_GROUP[match(rownames(pca_data), rownames(exp_clin))]
# [match(rownames(pca_data), rownames(exp_clin))] this is to make sure it matches the data in case the PCA has
# shuffled the samples during the process

# finally we plot the first two principal components
pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = SAMPLE_GROUP)) +
  geom_point(size = 3) + # this line of code says that we want to plot points with a size of 3
  labs(title = "PCA of Gene Expression for HC, Gout, and Sepsis", # we can add some labels to the plot
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()
pca

# we can save the plot
ggsave("D:/Uni/MASTER/Classes/Stats/assessment/pca_plot.png", pca, width = 8, height = 6, dpi = 300)

# each point corresponds to a sample (either Gout, HC, or Sepsis) and shows how that sample's 
# gene expression contributes to the principal components

# PCA is reducing the dimensionality of the gene expression data for each sample (based on 
# ~30,000 genes), but only samples are being plotted

# PCA plot helps answer whether the samples (patients) cluster into distinct groups based on their gene 
# expression profiles

# the plot shows that gout and sepsis samples are in different clusters, and gout is closer to HC than to 
# sepsis; we can now proceed to check the differentially expressed genes between gout and sepsis



# 4.2. Checking similarity of differentially expressed genes between gout and sa ------------------------

# we will use the dataframes we subset previously for significant genes (adjusted p-value < 0.05)
# we extract count of genes in each, which will be total significant differentially expressed genes in each
total_gout <- length(gout_sig$ID)
total_sa <- length(sa_sig$ID)

# now we check the overlap between the two dataframes
common_genes <- intersect(gout_sig$ID, sa_sig$ID)
common_length <- length(common_genes)

# now we get the length of unique genes (those in gout but not in sa, and vice versa)
unique_gout <- total_gout - common_length
unique_sa <- total_sa - common_length

# and we create a dataframe with the data for the significant differentially expressed genes (sdge)
summary_sdge <- data.frame(
  Category = c("Gout vs HC", "SA vs HC", "Overlap"),
  Total_Significant_Genes = c(total_gout, total_sa, common_length),
  Unique_Genes = c(unique_gout, unique_sa, NA)
)

print(summary_sdge)

# we can save this information to a csv file for further exploration for the report discussion
write.csv(summary_sdge, file = "D:/Uni/MASTER/Classes/Stats/assessment/summary_sdge.csv")

# the dataframe shows that there are 47 differentially expressed genes that
# overlap between the two, and 22 that are unique for Gout and 12999 that are unique for SA;
# this means that, while there are some overlapping, there are many differentially expressed 
# genes unique to each condition. 

# Knowing this, we can further explore which ones these are and if any of them may be used potentially 
# for diagnostic purposes



# 4.3. Identify non-overlapping significant genes that are differentially expressed ------------------------

# now we are only interested in the genes that are significant in gout but not in SA, so we use the subset again
# to subset gout_sig under the condition of ID that are NOT in sa_sig ID
gout_only <- subset(gout_sig, !(ID %in% sa_sig$ID))
# reorder based on p.adj
gout_only <- gout_only[order(gout_only$p.adj),]

# we repeat for sa
sa_only <- subset(sa_sig, !(ID %in% gout_sig$ID))
# reorder based on p.adj
sa_only <- sa_only[order(sa_only$p.adj),]


# we now check the data we obtained; the first rows are the most significant genes, since we have reordered the 
# dataframe according to p.adj values

# let's show only the first 5 rows, which are the 5 more significant genes (p.adj values are the lowest)
head(gout_only, 5)

# "ASPH" (ENSG00000198363), "IFT46" (ENSG00000118096), "SULT4A1" (ENSG00000130540), "CRIM1" (ENSG00000150938),  
# "TF" (ENSG00000091513) 

head(sa_only, 5)

# "AKR1B10" (ENSG00000198074), "KYNU" (ENSG00000115919), "PI3" (ENSG00000124102), "C10orf99" (ENSG00000188373), 
# "OASL" (ENSG00000135114)



# 4.4. Retrieve function information for the top genes *NOT INCLUDED IN REPORT* -----------------------------

#* *this section was finally not included in the report due to the word limit, but will leave it here for future reference*

# #* *Are these differentially expressed genes similar?*
# 
# # let's retrieve the function of the top significant unique genes so we can see if there are any genes that 
# # could be used for diagnostic purposes
# 
# #* *functions retrieved for all unique significant genes in gout and sa*
# # we will first get the function annot 
# top_genes <- unique(c(gout_only$ID, sa_only$ID))
# 
# # we had an error because 'keys' must be a character vector when attempting the select function from AnnotationDbi
# # so we need to ensure gene_symbols is a character vector
# top_genes <- as.character(top_genes)
# 
# # we now get the GO annotations for these genes; we want to make sure that we separate multiple GO terms that are 
# # stored in a single entry, and we flatten them into a single vector, as well as we remove duplicates to get a
# # list of unique GO terms; this is useful because we want to ensure that all GO terms associated with our genes
# # are considered but without duplicating any terms
# go_annots <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes, keytype = "ENSEMBL", 
#                                    columns = c("GO"))
# 
# # we need to ensure GO IDs are unique
# unique_go_ids <- unique(unlist(strsplit(go_annots$GO, split = ";")))
# 
# # now we retrieve GO term definitions for unique GO IDs
# go_definitions <- AnnotationDbi::select(GO.db, keys = unique_go_ids, keytype = "GOID", columns = "DEFINITION")
# 
# # now we combine GO annotations with definitions
# combined_annots <- left_join(go_annots, go_definitions, by = c("GO" = "GOID"))
# 
# # we need to summarize the GO terms for each gene symbol, because we have many rows for the same gene symbol with
# # all the functions associated with it; we want to collapse these into a single row
# collapsed_annots <- combined_annots %>%
#   group_by(ENSEMBL) %>%
#   summarise(DEFINITION = paste(unique(DEFINITION), collapse = "; ")) %>%
#   ungroup()
# 
# # finally, werge GO annotations with the gene expression data
# gout_func <- merge(gout_only, collapsed_annots, by.x = "ID", by.y = "ENSEMBL", all.x = TRUE)
# sa_func <- merge(sa_only, collapsed_annots, by.x = "ID", by.y = "ENSEMBL", all.x = TRUE)
# 
# # we add label gout or sepsis to each dataframe before merging
# gout_func$label <- "gout"
# sa_func$label <- "sepsis"
# 
# # we can now merge both dataframes to have all information in one, first 10 rows gout and last 10 sa
# allgenes_func <- rbind(gout_func, sa_func)
# 
# # we can now subset the dataframe to keep only the columns we are interested in to take a look at the functions
# gout_func <- gout_func[, c("ID", "symbol", "DEFINITION", "log2Fold", "p.adj")]
# sa_func <- sa_func[, c("ID", "symbol", "DEFINITION", "log2Fold", "p.adj")]
# 
# # we can save it to a csv file for further exploration for the report discussion
# write.csv(gout_func, file = "D:/Uni/MASTER/Classes/Stats/assessment/gout_func.csv")
# write.csv(sa_func, file = "D:/Uni/MASTER/Classes/Stats/assessment/sa_func.csv")
# 
# # we can take a look at these for the discussion section
# 
# # first we need to separate the dataframes based on string in DEFINITION separated with ;
# # we can use strsplit to split the strings and then unnest to separate them into different rows
# gout_func <- gout_func %>%
#   separate_rows(DEFINITION, sep = "; ")
# sa_func <- sa_func %>%
#   separate_rows(DEFINITION, sep = "; ")
# 
# # we will check the enriched biological processes for each gout and sepsis
# # first for gout sig
# gout_genes_func <- gout_sig$ID
# # now for sa sig. it has too many nodes for STRING so we will select the first 2000 genes
# # first we reorder the dataframe based on p.adj
# sa_sig <- sa_sig[order(sa_sig$p.adj),]
# sa_genes_func <- sa_sig[1:2000,]$ID
# # we print the first 5 lines in sa_genes_func
# head(sa_genes_func, 5)
# # we save these as txt
# write.table(gout_genes_func, file = "D:/Uni/MASTER/Classes/Stats/assessment/gout_genes_func.txt", row.names = FALSE)
# write.table(sa_genes_func, file = "D:/Uni/MASTER/Classes/Stats/assessment/sa_genes_func.txt", row.names = FALSE)
# 
# # we load the string enrichment analysis via GO in R
# enrichment_Process_goutsig <- read_excel("D:/Uni/MASTER/Classes/Stats/assessment/enrichment.Process_goutsig.xls")
# enrichment_Process_sasig <- read_excel("D:/Uni/MASTER/Classes/Stats/assessment/enrichment.Process_sasig.xls")
# 
# # Load enrichment data for your process
# # (assuming it’s already loaded in a similar structure as `enrichment_autism`)
# 
# # Change column "term description" to term_description if needed
# enrichment_Process_sasig <- enrichment_Process_sasig %>%
#   dplyr::rename(term_description = 'term description')
# 
# # Calculate the 90th percentile of the strength values to filter the dataframe
# percentile_90 <- quantile(enrichment_Process_sasig$strength, 0.90, na.rm = TRUE)
# threshold_sasig <- percentile_90
# 
# # Filter the data based on the threshold
# enrichment_Process_sasig_visu <- enrichment_Process_sasig %>% 
#   filter(strength > threshold_sasig)
# 
# # Calculate the total strength
# total_strength_sasig <- sum(enrichment_Process_sasig_visu$strength, na.rm = TRUE)
# 
# # Summarize data to identify duplicates and count occurrences
# summary_df_sasig <- enrichment_Process_sasig_visu %>%
#   group_by(term_description) %>%
#   summarize(
#     total_strength = sum(strength, na.rm = TRUE), # Sum strengths of duplicates
#     count = n(), # Count duplicates
#     mean_strength = mean(strength, na.rm = TRUE) # Mean strength of duplicates
#   ) %>%
#   ungroup()
# 
# # Merge this summary back with the original data to tag duplicates with their counts
# enrichment_Process_sasig_visu <- enrichment_Process_sasig_visu %>%
#   left_join(summary_df_sasig, by = "term_description")
# 
# # Create the label_without_percentage column, taking into account the PPI counts
# enrichment_Process_sasig_visu$label_without_percentage <- paste(
#   enrichment_Process_sasig_visu$term_description, 
#   "(", round(enrichment_Process_sasig_visu$mean_strength, 2), ", ", 
#   enrichment_Process_sasig_visu$count, " PPI)"
# )
# 
# # Prepare data for the bar plot with raw strength values
# agg_data_strength_sasig <- enrichment_Process_sasig_visu %>%
#   ungroup() %>%
#   dplyr::select(label_without_percentage, strength) %>%
#   distinct() %>%
#   mutate(label_without_percentage = factor(label_without_percentage, levels = unique(label_without_percentage[order(strength)])))
# 
# # Update the color palette to be color-blind friendly
# num_colors_needed_sasig <- length(unique(agg_data_strength_sasig$label_without_percentage))
# color_palette_sasig <- colorRampPalette(brewer.pal(min(num_colors_needed_sasig, 8), "Set2"))(num_colors_needed_sasig)
# colors_agg_strength_sasig <- setNames(color_palette_sasig, levels(agg_data_strength_sasig$label_without_percentage))
# 
# # Horizontal bar plot with raw strength values
# bar_plot_strength_sasig <- ggplot(agg_data_strength_sasig, aes(y = label_without_percentage, x = strength, fill = label_without_percentage)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = colors_agg_strength_sasig) +
#   labs(
#     y = "Biological Processes",
#     x = "Strength",
#     fill = "Biological Process") +
#   xlim(0, 1.0) + 
#   theme(axis.text.y = element_text(size = 10),
#         legend.position = "none")
# 
# # Display the plot
# print(bar_plot_strength_sasig)
# 
# # Save the plot
# ggsave("D:/Uni/MASTER/Classes/Stats/assessment/bar_plot_strength_sasig.png", 
#        plot = bar_plot_strength_sasig, height = 3, width = 18, dpi = 300)
# 
# # we retrieve the list of the biological processes that are enriched in the significant genes for sa
# enrich_sa <- enrichment_Process_sasig_visu$term_description
# 
# # we can save this information to a txt file for further exploration for the report discussion
# write.table(enrich_sa, file = "D:/Uni/MASTER/Classes/Stats/assessment/enrich_sa.txt", row.names = FALSE)
# 
# 
# #* *functions retrieved for top 5 unique significant genes in gout and sa*
# # for this, we will first combine both dataframes and get the unique genes
# top_genes <- unique(c(gout_only$ID[1:5], sa_only$ID[1:5]))
# 
# # we had an error because 'keys' must be a character vector when attempting the select function from AnnotationDbi
# # so we need to ensure gene_symbols is a character vector
# top_genes <- as.character(top_genes)
# 
# # we now get the GO annotations for these genes; we want to make sure that we separate multiple GO terms that are 
# # stored in a single entry, and we flatten them into a single vector, as well as we remove duplicates to get a
# # list of unique GO terms; this is useful because we want to ensure that all GO terms associated with our genes
# # are considered but without duplicating any terms
# go_annots <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes, keytype = "ENSEMBL", 
#                                    columns = c("GO"))
# 
# # we need to ensure GO IDs are unique
# unique_go_ids <- unique(unlist(strsplit(go_annots$GO, split = ";")))
# 
# # now we retrieve GO term definitions for unique GO IDs
# go_definitions <- AnnotationDbi::select(GO.db, keys = unique_go_ids, keytype = "GOID", columns = "DEFINITION")
# 
# # now we combine GO annotations with definitions
# combined_annots <- left_join(go_annots, go_definitions, by = c("GO" = "GOID"))
# 
# # we need to summarize the GO terms for each gene symbol, because we have many rows for the same gene symbol with
# # all the functions associated with it; we want to collapse these into a single row
# collapsed_annots <- combined_annots %>%
#   group_by(ENSEMBL) %>%
#   summarise(DEFINITION = paste(unique(DEFINITION), collapse = "; ")) %>%
#   ungroup()
# 
# # we want to add this information to the top genes dataframes, so we need to subset the top genes
# gout_only_top5 <- gout_only[1:5, ]
# sa_only_top5 <- sa_only[1:5, ]
# 
# # finally, werge GO annotations with the gene expression data
# gout_func <- merge(gout_only_top5, collapsed_annots, by.x = "ID", by.y = "ENSEMBL", all.x = TRUE)
# sa_func <- merge(sa_only_top5, collapsed_annots, by.x = "ID", by.y = "ENSEMBL", all.x = TRUE)
# 
# # we can now merge both dataframes to have all information in one, first 10 rows gout and last 10 sa
# top5_func <- rbind(gout_func, sa_func)
# 
# # we have a few NA values in the definition column, but we can take a look at the functions annotated
# # and get a better idea of what these genes do and if there are any genes that stand out for diagnostic purposes
# 
# # we can add a column label for gout and sa, and retrieve a smaller dataframe with only symbol, definition, and label
# # we create the label column
# top5_func$label <- "unknown"
# 
# # we know by how rbind works that the first 5 rows correspond to gout and the other 5 to sa, so with this information
# # we can add the label names
# top5_func$label[1:5] <- "gout"
# top5_func$label[6:10] <- "sepsis"
# 
# # we can now subset the dataframe to keep only the columns we are interested in to take a look at the functions
# top5_class <- top5_func[, c("ID", "label", "symbol", "DEFINITION", "log2Fold", "p.adj")]
# top5_class_nd <- top5_func[, c("ID", "label", "symbol", "log2Fold", "p.adj")]
# 
# # we can save it to a csv file for further exploration for the report discussion
# write.csv(top5_class, file = "D:/Uni/MASTER/Classes/Stats/assessment/top5_definitions.csv")
# write.csv(top5_class_nd, file = "D:/Uni/MASTER/Classes/Stats/assessment/top5_nd.csv")
# 
# # in general, there are some shared functions in cellular processes, extracellular functions, and immune response 
# # (which is to be expected)
# 
# # but there are also some functions found specifically in the top 5 differentially expressed genes in gout and sa;
# # gout has functions related to iron ion binding and transport:
# # sepsis has functions related to tryptophan and NAD+ Pathways, and specific immune pathways
# # note: this question was answered with the help of data exploration analysis tool (ChatGPT-4o) and double-checked
# # by human review for efficiency reasons
# 
# # we need to remember that the genes we are looking at are unique to each condition, but these are the top 5
# # significant ones and a more extensive research on every function of unique genes to each condition would
# # be needed to properly determine unique functions for each condition that could hint at diagnostic purposes



# 4.5. Are there any genes that are significantly different between Gout and SA? ----------------

# we will use the package DESeq2 to check for differentially expressed genes between gout and sa, using the expression
# data we have in the exp_clin dataframe (samples as rows and columns as genes + sex, neutrophils, sample group, and sample)

#* @TASK4-Part2

# documentation for DESeq2
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# first we need to subset the data to retain only Gout and Sepsis samples for DESeq2 analysis
countData <- exp_clin %>%
  filter(SAMPLE_GROUP %in% c("GOUT", "SEPSIS")) %>%
  dplyr::select(-SAMPLE_GROUP, -SEX, -NEUTROPHILS)  # we also need to remove non-gene columns
# we use dplyr::select to tell r that we want to use the select function from dplyr package, because another
# function is masking this one and r cannot find it otherwise

# now we want to transpose the countData so that genes are rows and samples are columns
countData_t <- t(countData)

# first row is the sample names so we remove it
countData_t <- countData_t[-1,]

countData <- as.matrix(countData_t)

# we need to define the condition factor for Gout and Sepsis, this is to specify the groups for the analysis
condition <- factor(c(rep("GOUT", 9), rep("SEPSIS", 9)))

# now we create colData dataframe for DESeq2
colData <- data.frame(row.names = colnames(countData), condition = condition)

# we were also getting an error because the data was not numeric so we will force a numeric conversion
countData <- apply(countData, 2, as.numeric)

# construct DESeqDataSet object since DESeq analysis requires this object; we also round numeric countData data and
# then a matrix that we can use as an object for DESeq2; colData is the table with sample information; design
# indicates how to model the samples, in out case condition is the only variable we are interested in
dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(countData)), colData = colData, design = ~ condition)

# now we run DESeq2 to perform differential expression analysis between the two conditions (Gout and Sepsis)
dds_analysis <- DESeq(dds)

# we extract the results (comparing conditions)
res <- results(dds_analysis)
summary(res)

# we convert the results to a dataframe so we can work with it
res <- as.data.frame(res)

# write results to a CSV file for future reference
write.csv(res, file = "D:/Uni/MASTER/Classes/Stats/assessment/DESeq2_results.csv")

# add row names to the 'res' dataframe.
res$ID <- rownames(countData_t)
head(res)

# set ID as the first column in the dataframe res
res <- res %>% relocate(ID, .before = 1)

# we will add the symbol information to the dataframe, using the annot dataframe
# we will merge the dataframes by ID but adding only the symbol column

res <- merge(x = res, y = annot[, c("ID", "symbol")], by = "ID", all.x = TRUE)
res <- res %>% relocate(symbol, .before = 2)

# now we extract the top 10 significant genes
significant_genes <- res[which(res$padj < 0.05), ]

# we reorder by padj values from lower to higher, so we make sure the first rows are the most significant genes
significant_genes <- significant_genes[order(significant_genes$padj),]
head(significant_genes, 10)

# we can save this information to a csv file for further exploration for the report discussion
write.csv(significant_genes, file = "D:/Uni/MASTER/Classes/Stats/assessment/significant_genes.csv")

# the top 10 significant genes are the ones with the lowest p.adj values, which are the most significant ones
# we can get a list of these genes
top10_genes <- significant_genes$ID[1:10]
top10_genes
top10_genes_s <- significant_genes$symbol[1:10]
top10_genes_s

# we add the symbol column from annot dataframe to the res dataframe
# these are in ID:
# "ENSG00000184330" "ENSG00000124102" "ENSG00000198074" "ENSG00000163220" "ENSG00000115919" "ENSG00000206073"
# "ENSG00000134827" "ENSG00000135373" "ENSG00000241794" "ENSG00000153802"
#in HGNC symbol:
# "S100A7A"   "PI3"       "AKR1B10"   "S100A9"    "KYNU"      "SERPINB4"  "TCN1"      "EHF"       "SPRR2A"   
# "TMPRSS11D"

# since the first 9 columns of our data were for gout samples and the next 9 were for sepsis
# the analysis shows how gene expression in SA samples compared to Gout samples. That is to say:
#     -Positive log2 fold changes will indicate genes that are upregulated in SA relative to Gout
#     -Negative log2 fold changes will indicate genes that are downregulated in SA relative to Gout

# this means that, i.e. for gene ENSG00000184330 with log2FoldChange of 10.057751, it is upregulated in SA relative to Gout

# we will retrieve the rows for the 10 top genes in this analysis for the results section
top10_genes_res <- res[res$ID %in% top10_genes, ]

# reorder by p.adj
top10_genes_res <- top10_genes_res[order(top10_genes_res$padj),]

# we can save this information to a csv file for further exploration for the report discussion
write.csv(top10_genes_res, file = "D:/Uni/MASTER/Classes/Stats/assessment/top10_genes_res.csv")

# volcano plot documentation: https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

volc_plot <- EnhancedVolcano::EnhancedVolcano(res,
                                 lab = res$symbol,
                                 x = 'log2FoldChange',
                                 y = 'pvalue')
volc_plot

# volcano plot: 
# log2 fold change (magnitude of gene expression change between the two conditions) on the x-axis 
# against the -log10 p-value (statistical significance) on the y-axis - Higher values on the 
# y-axis indicate stronger statistical evidence

# dot colors in the plot:
#    - Red points indicate genes that are both statistically significant (p-value) and show a 
#    significant fold change (either up or down). These genes are of most interest in our case
#    - Blue points represent genes with significant p-values but without a large fold change
#    - Green points show genes with significant fold changes but without a significant p-value
#    - Gray points are genes that are not statistically significant in terms of both fold 
#    change and p-value

# dashed lines:
# These represent thresholds. 
# The horizontal line is often placed at a p-value threshold,  
# and the vertical lines is the log2 fold change cutoffs (+1 or -1 for upregulated and 
# downregulated genes)

# we save the plot for further discussion
ggsave("D:/Uni/MASTER/Classes/Stats/assessment/volc_plot.png", width = 8, height = 6, dpi = 300)


#* *Comparison between gout VS SA and gout and SA vs healthy*
# we can compare the significant genes we obtained in the DESeq2 analysis with the unique significant genes we
# filtered from the differential gene expression data for gout vs healthy and sepsis vs healthy

# we will use the 10 significant genes from DESeq2 analysis in the dataframe top10_genes_s and check what results we 
# got from the unique significant genes in gout and sa dataframes

# in order to do this, we need to obtain the gene ID for top 10 significant differentially expressed genes in sepsis
# versus healthy (top10_genes - we have already done this in previous steps), and then retrieve the rows in the gout 
# and sa dataframes that correspond to these genes

# we check which of these genes are in unique gout dataframe and retrieve the rows (subsetting)
top10_in_gout <- gout[gout$ID %in% top10_genes, ]

top10_in_sa <- sa[sa$ID %in% top10_genes, ]

# we add symbol information to the dataframes from the annot dataframe
top10_in_gout <- merge(x = top10_in_gout, y = annot[, c("ID", "symbol")], by = "ID", all.x = TRUE)
top10_in_sa <- merge(x = top10_in_sa, y = annot[, c("ID", "symbol")], by = "ID", all.x = TRUE)

# recolacte to second column position
top10_in_gout <- top10_in_gout %>% relocate(symbol, .before = 2)
top10_in_sa <- top10_in_sa %>% relocate(symbol, .before = 2)

# we round the p and p.adj values to 2 decimal places for better readability
top10_in_gout$p <- round(top10_in_gout$p, 2)
top10_in_gout$p.adj <- round(top10_in_gout$p.adj, 2)

# we won't round values for sepsis because the values are too small

print(top10_in_gout)
print(top10_in_sa)

# we save this information to a csv file for further exploration for the report discussion
write.csv(top10_in_gout, file = "D:/Uni/MASTER/Classes/Stats/assessment/top10_in_gout.csv")
write.csv(top10_in_sa, file = "D:/Uni/MASTER/Classes/Stats/assessment/top10_in_sa.csv")

# this corroborates the findings from the DESeq2 analysis, as some of the top 10 significant genes identified are specifically 
# differentially expressed in sepsis versus healthy samples with log2Fold changes ranging from 2.79 to 9.06, and none was in 
# gout versus healthy samples

# this provides further confidence in the results obtained from the DESeq2 analysis, as the top 10 significant genes, and also 
# key genes in the SA expression profile

# we will now compare the 10 top genes in the DESeq2 analysis with the top 10 genes in sa dataframe, since we have performed a
# comparative analysis between gene expression in gout and sepsis samples with gout as reference
# also, the most significantly differentially expressed genes with higher absolute log2Fold changes are found in sa

# we need to retrieve the top 10 significant genes in sa  dataframes and compare them with the top 10 significant
# we first reorder sa dataframe based on p.adj from lower to higher
sa_only <- sa_only[order(sa_only$p.adj),]
# now we print the first 10 rows and save the dataframe only for these 10 first rows
head(sa_only, 10)
write.csv(sa_only[1:10,], file = "D:/Uni/MASTER/Classes/Stats/assessment/sa_top10.csv")

# by performing the DESeq2 analysis, we have identified a different set of genes that are found to be the most
# significantly differentially expressed between gout and sepsis samples, which are present in the significant genes
# for sepsis versus healthy, but not ranked in the same order; 



# 4.6. Retrieve function information for the top genes in DESeq2 analysis ----------------------------

# for this, we will get the top 10 genes
top10_genes <- unique(top10_genes_res$ID)

# we had an error because 'keys' must be a character vector when attempting the select function from AnnotationDbi
# so we need to ensure ID is a character vector
top10_genes <- as.character(top10_genes)

# we now get the GO annotations for these genes; we want to make sure that we separate multiple GO terms that are 
# stored in a single entry, and we flatten them into a single vector, as well as we remove duplicates to get a
# list of unique GO terms; this is useful because we want to ensure that all GO terms associated with our genes
# are considered but without duplicating any terms
go10_annots <- AnnotationDbi::select(org.Hs.eg.db, keys = top10_genes, keytype = "ENSEMBL", 
                                   columns = c("GO"))

# we need to ensure GO IDs are unique
unique_10go_ids <- unique(unlist(strsplit(go10_annots$GO, split = ";")))

# now we retrieve GO term definitions for unique GO IDs
go10_definitions <- AnnotationDbi::select(GO.db, keys = unique_10go_ids, keytype = "GOID", columns = "DEFINITION")

# now we combine GO annotations with definitions
combined_10annots <- left_join(go10_annots, go10_definitions, by = c("GO" = "GOID"))

# we need to summarize the GO terms for each gene symbol, because we have many rows for the same gene symbol with
# all the functions associated with it; we want to collapse these into a single row
collapsed_10annots <- combined_10annots %>%
  group_by(ENSEMBL) %>%
  summarise(DEFINITION = paste(unique(DEFINITION), collapse = "; ")) %>%
  ungroup()


# finally, werge GO annotations with the gene expression data
top10_func <- merge(top10_genes_res, collapsed_10annots, by.x = "ID", by.y = "ENSEMBL", all.x = TRUE)
top10_func <- top10_func %>% relocate(DEFINITION, .before = 3)

#we now delete the columns we do not need for the results section
top10_func <- top10_func[, c("ID", "symbol", "DEFINITION", "log2FoldChange", "pvalue", "padj")]

# we can save this information to a csv file for further exploration for the report discussion
write.csv(top10_func, file = "D:/Uni/MASTER/Classes/Stats/assessment/top10_func.csv")



# 5. Plot distributions for significant genes in gout and sepsis --------------------------------

#* @TASK2

# we can plot the distribution of log2Fold changes for the top 10 significant genes in gout and sa
# we will use boxplots and density plots to visualize the distribution of log2Fold changes for these genes

# use vRColorBrewer for colorblind-friendly colours

# defining a colorblind-friendly palette with `RColorBrewer`
colors <- brewer.pal(n = 3, name = "Set2")  # 'Set1' is colorblind-friendly

# boxplot for significant genes in Gout vs HC using colorblind-friendly colors
ggplot(gout_only, aes(x = "Gout vs HC", y = log2Fold)) +
  geom_boxplot(fill = colors[1], alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7, color = colors[2]) +
  labs(title = "Boxplot of Log2 Fold Change in Gout vs HC for significant genes", x = "", y = "Log2 Fold Change") +
  theme_minimal()

# boxplot for significant genes in SA vs HC
ggplot(sa_only, aes(x = "SA vs HC", y = log2Fold)) +
  geom_boxplot(fill = colors[2], alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7, color = colors[3]) +
  labs(title = "Boxplot of Log2 Fold Change in SA vs HC for significant genes", x = "", y = "Log2 Fold Change") +
  theme_minimal()

# create the boxplot for both distributions comparison, so it is easier to see
# got the idea from: https://stackoverflow.com/questions/45150841/boxplots-from-two-dataframes-in-r
# although I didn't use reshape, but I understood how to do that by binding the dataframes and using ggplot
# Combine the data for both Gout vs HC and SA vs HC

# label rows
gout_only$Condition <- "Gout vs HC"
sa_only$Condition <- "SA vs HC"

# use new dataframe for this
gout_selected <- gout_only %>%
  dplyr::select(ID, log2Fold, last_col())  # ID, log2Fold, and last column

sa_selected <- sa_only %>%
  dplyr::select(ID, log2Fold, last_col())

# merge by adding rows
boxplotdf <- bind_rows(gout_selected, sa_selected)

# create the boxplot using the colorblind-friendly palette we defined earlier
bp_goutsa <- ggplot(boxplotdf, aes(x = Condition, y = log2Fold, fill = Condition)) +
  geom_boxplot() +
  labs(title = "Significative Gene Expression Differences: Gout vs HC and SA vs HC",
       x = "Group", y = "log2 Fold Change") +
  scale_fill_brewer(palette = "Set3") +  # Using color-blind-friendly palette
  theme_minimal()
bp_goutsa

# violin plot may be more informative in this case so we will also plot a violin plot
# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

vp_goutsa <- ggplot(boxplotdf, aes(x = Condition, y = log2Fold, fill = Condition)) +
  geom_violin(alpha = 0.5, trim = FALSE) +  # Add violin plot with some transparency
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +  # Overlay boxplot without outliers
  labs(title = "              Significant Gene Expression Differences: 
       Gout vs Healthy and Septic arthritis (SA) vs Healthy",
       x = "Group", y = "log2 Fold Change") +
  scale_fill_brewer(palette = "Set3") +  # Using color-blind-friendly palette
  theme_minimal()
vp_goutsa

# we save the plot
ggsave("D:/Uni/MASTER/Classes/Stats/assessment/boxplot_goutsa.png", bp_goutsa, width = 8, height = 6, dpi = 300)
ggsave("D:/Uni/MASTER/Classes/Stats/assessment/violplot_goutsa.png", vp_goutsa, width = 8, height = 6, dpi = 300)

# density plot for Gout vs HC
dens_gout <- ggplot(gout_only, aes(x = log2Fold)) +
  geom_density(fill = colors[1], alpha = 0.5) +
  labs(title = "Density Plot of significative Log2 Fold Change in Gout vs HC", x = "Log2 Fold Change", y = "Density") +
  theme_minimal()
dens_gout

# density plot for SA vs HC
dens_sa <- ggplot(sa_only, aes(x = log2Fold)) +
  geom_density(fill = colors[2], alpha = 0.5) +
  labs(title = "Density Plot of significative Log2 Fold Change in SA vs HC", x = "Log2 Fold Change", y = "Density") +
  theme_minimal()
dens_sa

# we can save the plots
ggsave("D:/Uni/MASTER/Classes/Stats/assessment/density_gout.png", dens_gout, width = 8, height = 6, dpi = 300)
ggsave("D:/Uni/MASTER/Classes/Stats/assessment/density_sa.png", dens_sa, width = 8, height = 6, dpi = 300)

# the expression data for the unique significant genes in gout is much lower than in sepsis; in gout the expression 
# values fall between -1.5 and 3.5 while in sepsis they fall between -8.5 and 11.5, which is highly negatively/positively 
# expressed



# 6. Exploring the distribution of clinical data -------------------------------------------------

#* @TASK1

# we will subset sample dataframe according to each condition so we can plot the distribution of F and M for each one
sub_gout <- sample %>% filter(SAMPLE_GROUP == "GOUT")

# we want to make a dataframe with the summary data for gout: mean and standard deviation for neutrophils, and counts 
# for M and F
summary_data <- sample %>%
  group_by(SAMPLE_GROUP) %>%
  summarize(Mean_Neutrophils = mean(NEUTROPHILS, na.rm = TRUE),
            SD_Neutrophils = sd(NEUTROPHILS, na.rm = TRUE),
            Male_Count = sum(SEX == "M"),
            Female_Count = sum(SEX == "F"))

# we do fischer test for sex distribution differences because we have a small size of samples (counts)
sex_table <- table(sample$SEX, sample$SAMPLE_GROUP)
fisher_test_sex <- fisher.test(sex_table)
print(fisher_test_sex)

# ANOVA test for Neutrophils between the groups
anova_neutrophils <- aov(NEUTROPHILS ~ SAMPLE_GROUP, data = sample)
anova_summary <- summary(anova_neutrophils)

# we check if we can use ANOVA for neutrophils distribution between the groups

# residuals from the ANOVA model
residuals_anova <- residuals(anova_neutrophils)

# QQ plot
qqnorm(residuals_anova)
qqline(residuals_anova)

# Shapiro-Wilk test
shapiro_test <- shapiro.test(residuals_anova)
print(shapiro_test)
print(shapiro_test$p.value)

# the residuals are normally distributed

# Levene's test for homogeneity of variance
levene_test <- leveneTest(NEUTROPHILS ~ SAMPLE_GROUP, data = sample)
print(levene_test)

# the homogeneous variance assumption is met

# then we can accept ANOVA as a valid test for the data and see the results for sex and neutrophil counts
# significant difference (not matched data) or no significant difference (matched data)

print(summary_data)
print(fisher_test_sex)
print(anova_summary)

# we can save the summary data to a csv file for further exploration for the report discussion
write.csv(summary_data, file = "D:/Uni/MASTER/Classes/Stats/assessment/summary_data.csv")

#* *Conclusion for part 6*
# sex has shown no statistically significant differences between the groups (p-value=1)
# the mean neutrophils count in summary table is higher in sepsis samples than in gout and hc samples. 
# The difference in neutrophils count between hc, gout and sa in the anova test is shown to be statistically 
# significant
# therefore, we can conclude that the clinical data is not completely fairly distributed between the different 
# groups of samples; however this difference may be due to the underlying conditions and their impact on the 
# immune system and response. Since hc counts are lower than gout and sa, and sepsis has the higher counts, we 
# can assume that this may be the reason for the difference in neutrophil counts



# 7. Exploring the clinical data impact in gene expression --------------------------------------------

#* @TASK3

# now we want to explore if there is any relationship between the clinical data and the gene expression data, so we can take
# this into account when replying to our main research question 

# TRY USING GENE EXP FOR ALL SAMPLES FOR EACH GENE IN TOP 5

#* *Preparing the data*

# Subset gout samples and corresponding genes
gout_samples <- exp_clin[exp_clin$SAMPLE_GROUP == "GOUT", ]
gout_columns <- as.character(gout_only_top5$ID)
gout_genes <- gout_samples %>% dplyr::select(SAMPLE, SAMPLE_GROUP, SEX, NEUTROPHILS, any_of(gout_columns))

# Subset SA samples and corresponding genes
sa_samples <- exp_clin[exp_clin$SAMPLE_GROUP == "SEPSIS", ]
sa_columns <- as.character(sa_only_top5$ID)
sa_genes <- sa_samples %>% dplyr::select(SAMPLE, SAMPLE_GROUP, SEX, NEUTROPHILS, any_of(sa_columns))

# we have now two dataframes prepared for the linear model, one for gout and one for sa, with the clinical data and
# the expression data for each sample and gene of the top 5 significant genes that are uniquely differentially 
# expressed in each condition



#* *Fit a linear model to explore the relationship between clinical data and gene expression for specific genes*

# we ca apply the linear model for each gene separately for the 5 genes for gout and 5 for sa

# however, to apply the multiple linear model to explore the relationship between clinical data and gene 
# expression, we need to first check if the data meets the assumptions of the linear model

# we can check the distribution of the residuals, the linearity of the relationship between the dependent and
# independent variables, the homoscedasticity of the residuals, and the normality of the residuals with the 
# model diagnostics plots

#* *LM for gout genes*

# first retrieve a list of the ID for the top 10 genes uniquely differentially expressed in gout
gout_list <- gout_only_top5$ID
gout_list

# we need to apply an lm for each gene and, separately, to test the relationship with sex and neutrophil counts
# we also want to plot the diagnostics for each model to check if the assumptions of the linear model are met

# we will write a function to do all this looping for each gene in the list

separate_lms <- function(gene_id, data) {
  # linear model for SEX
  lm_sex <- lm(as.formula(paste(gene_id, "~ SEX")), data = data) # lm for each gene in the list and SEX

  # we apply anova to check if the model is significant
  anova_sex <- anova(lm_sex)
  
  # print the results
  print(summary(lm_sex))
  print(anova_sex)
  
  # diagnostics for SEX with the label for gene and SEX model
  print(autoplot(lm_sex) + ggtitle(paste(gene_id, "SEX")))
  
  # linear model for neutrophils
  lm_neutrophils <- lm(as.formula(paste(gene_id, "~ NEUTROPHILS")), data = data) # lm for each gene in 
                        
  anova_neutrophils <- anova(lm_neutrophils)
  
  print(summary(lm_neutrophils))
  print(anova_neutrophils)
  
  # diagnostics for Neutrophils with the label for gene and neutrophils model
  print(autoplot(lm_neutrophils) + ggtitle(paste(gene_id, "NEUTROPHILS")))
}

# Note to self: The ANOVA test compares the model including the predictor (e.g., SEX or NEUTROPHILS) against a 
# model without it (essentially, an intercept-only model). The F-statistic and the associated p-value from the 
# ANOVA test indicate whether the predictor variable has a statistically significant effect on the gene expression 
# values. A low p-value (typically < 0.05) suggests that the predictor contributes significantly to the model
# http://www.stat.yale.edu/Courses/1997-98/101/anovareg.htm
# While the summary of the linear model provides details on individual coefficients, the ANOVA test gives an 
# overall indication of the predictor’s contribution to the model's explanatory power. In this context, the ANOVA 
# output's F-statistic and p-value allow us to evaluate the model's fit and whether adding the predictor makes a 
# meaningful improvement over a simpler model


# now we can use this function to run the linear models for each gene in the gout list
for (gene in gout_list) {
  cat("Running models for gene:", gene, "\n") # this is to include a wee title for the gene and show the results 
                                              # in the next line
  separate_lms(gene, gout_genes)
}

# sex variable looks in general reasonably well-behaved, and it has shown to not be significant to explain the
# expression values of the genes in the list
# however, neutrophils variable does not seem to meet the assumptions, so we may need to apply a 
# non-parametric approach

#* *Non-parametric approach for neutrophils for gout top 5 genes*

# we will divide neutrophil counts into equally sized intervals and assign level scores

# we will use a spearman correlation test to check the relationship between neutrophils and gene expression
# since this is a non-parametric test that checks the relationship between two numerical variables

spearman1 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000198363, method = "spearman", exact=FALSE)
print(spearman1)

# we include "exact=FALSE" to avoid the error message "Cannot compute exact p-value with ties" that appears when
# there are duplicated values in the data, the calculated p-value is an approximation

# Spearman’s correlation coefficient is a statistical measure of the strength of a monotonic relationship between 
# paired data. Its interpretation is similar to that of Pearsons. Correlation is an effect size and so we can
# verbally describe the strength of the correlation using the following guide for the absolute value of:
# 
# - .00-.19 “very weak”
# - .20-.39 “weak”
# - .40-.59 “moderate”
# - .60-.79 “strong”
# - .80-1.0 “very strong”

# p-value assesses whether this observed correlation is statistically significant. Even is the rho value indicated a weak
# correlation, if p-value is not significant then there isn’t enough evidence to conclude a meaningful association between 
# neutrophil counts and the gene’s expression level in the gout samples

spearman2 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000118096, method = "spearman", exact=FALSE)
print(spearman2)

spearman3 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000130540, method = "spearman", exact=FALSE)
print(spearman3)

spearman4 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000150938, method = "spearman", exact=FALSE)
print(spearman4)

spearman5 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000091513, method = "spearman", exact=FALSE)
print(spearman5)

# neutrohpils do not appear to predict the gene expression values for the top 5 genes in gout, as the p-values are
# not significant for any of the genes

#* *LM for sa genes*

sa_list <- sa_only_top5$ID
sa_list

# we repeat the loop for sepsis list
for (gene in sa_list) {
  cat("Running models for gene:", gene, "\n")
  separate_lms(gene, sa_genes)
}

# again, sex variable looks in general reasonably well-behaved, and it has shown to not be significant to explain 
# the expression values of the genes in the list
# however, neutrophils variable does not seem to meet the assumptions, so we may need to apply a 
# non-parametric approach

#* *Non-parametric approach for neutrophils for sa top 5 genes*

# we will use a spearman correlation test again

spearman6 <- cor.test(sa_genes$NEUTROPHILS, sa_genes$ENSG00000198074, method = "spearman", exact=FALSE)
print(spearman6)

spearman7 <- cor.test(sa_genes$NEUTROPHILS, sa_genes$ENSG00000115919, method = "spearman", exact=FALSE)
print(spearman7)

spearman8 <- cor.test(sa_genes$NEUTROPHILS, sa_genes$ENSG00000124102, method = "spearman", exact=FALSE)
print(spearman8)

spearman9 <- cor.test(sa_genes$NEUTROPHILS, sa_genes$ENSG00000188373, method = "spearman", exact=FALSE)
print(spearman9)

spearman10 <- cor.test(sa_genes$NEUTROPHILS, sa_genes$ENSG00000135114, method = "spearman", exact=FALSE)
print(spearman10)

# again, neutrohpils do not appear to predict the gene expression values for the top 5 genes in gout, as the 
# p-values are not significant for any of the genes

#* *conclusion for part 7*
# the expression data is not significantly correlated with the clinical data for these genes
# at least for the data provided, age and other clinical data could be more relevant for the expression
# of these genes, but we do not have this data



# 8. Exploring the relationship between gene expression and clinical data for the top 10 genes in the DESeq2 analysis --------------

# we can additionally test the relationship between the gene expression and the clinical data for the top two 
# significant genes found to be differentially expressed in sepsis versus gout, which are ENSG00000184330 (S100A7A)
# and ENSG00000124102 (PI3)

# we need to add the data corresponding to these genes to the dataframes gout_genes and sa_genes for the
# function we wrote to work properly
# we need to find these genes in the columns of gout_samples and sa_samples and add them to the gout_genes and sa_genes dataframe
gout_merge <- gout_samples %>% dplyr::select(ENSG00000124102, ENSG00000184330)
gout_genes <- cbind(gout_genes, gout_merge)

# ENSG00000124102 was already in the dataframe

sa_merge <- sa_samples %>% dplyr::select(ENSG00000184330)
sa_genes <- cbind(sa_genes, sa_merge)

# we will extract the name of the two first genes in the dataframe significant_genes based on ID
top_genes <- significant_genes$ID[1:2]
top_genes # checking that they were selected correctly

# now we apply the function to these genes for gout and sa dataframes
for (gene in top_genes) {
  cat("Running models for gene:", gene, "\n") # this is to include a wee title for the gene and show the results 
  # in the next line
  separate_lms(gene, gout_genes)
}

# non-parametric test for neutrophil count
spearman11 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000124102, method = "spearman", exact=FALSE)
print(spearman11)
spearman12 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000184330, method = "spearman", exact=FALSE)
print(spearman12)

# Sex nor neutrophil count seem to predict the gene expression values for the top 2 genes found to be differentially 
# expressed in sepsis versus gout

for (gene in top_genes) {
  cat("Running models for gene:", gene, "\n") # this is to include a wee title for the gene and show the results 
  # in the next line
  separate_lms(gene, sa_genes)
}

# non-parametric test for neutrophil count
spearman13 <- cor.test(sa_genes$NEUTROPHILS, sa_genes$ENSG00000124102, method = "spearman", exact=FALSE)
print(spearman13)
spearman14 <- cor.test(sa_genes$NEUTROPHILS, sa_genes$ENSG00000184330, method = "spearman", exact=FALSE)
print(spearman14)

# Sex nor neutrophil count seem to predict the gene expression values for the top 2 genes in sepsis for the top 2
# genes found to be differentially expressed in sepsis versus gout

#* *REPEATING THIS FOR THE OTHER 8 GENES IN THE TOP 10 OF DESeq2 ANALYSIS*
#we will also check these for the other top 10 genes from the DESeq2 analysis

# we need to add the data corresponding to these genes to the dataframes gout_genes and sa_genes for the
# function we wrote to work properly
# we need to find these genes in the columns of gout_samples and sa_samples and add them to the gout_genes and sa_genes dataframe
gout_merge <- gout_samples %>% dplyr::select(ENSG00000198074, ENSG00000163220, ENSG00000115919, ENSG00000206073)
gout_genes <- cbind(gout_genes, gout_merge)

sa_merge <- sa_samples %>% dplyr::select(ENSG00000198074, ENSG00000163220, ENSG00000115919, ENSG00000206073)
sa_genes <- cbind(sa_genes, sa_merge)

# checking if there is any repeated gene column in either dataframe gout_genes or sa_genes
any(duplicated(colnames(gout_genes)))
any(duplicated(colnames(sa_genes)))

# we retrieve the duplicated genes in sa_genes
duplicated_genes <- colnames(sa_genes)[duplicated(colnames(sa_genes))]
duplicated_genes

# we remove the duplicated genes from sa_genes but only once, the first column found with that name
# first we get the index of each duplicated column's first occurrence
first_duplicates <- which(duplicated(colnames(sa_genes), fromLast = TRUE))
# what this does is this retrieves the index of the first occurrence of each duplicated column name
# so it only deletes those, and not all columns with the same name

# now we subset the dataframe by removing only the first occurrence of duplicated columns
sa_genes <- sa_genes[, -first_duplicates] # we use the index to remove the columns

# check again just in case
any(duplicated(colnames(sa_genes)))

# all good!

# we need to make a list with the next 4 genes in the top 10 genes from the DESeq2 analysis
top10_genes
# "ENSG00000198074" "ENSG00000163220" "ENSG00000115919" "ENSG00000206073"
top_genes2 <- significant_genes$ID[3:6]
top_genes2 # checking that they were selected correctly

# now we apply the function to these genes for gout and sa dataframes
for (gene in top_genes2) {
  cat("Running models for gene:", gene, "\n") # this is to include a wee title for the gene and show the results 
  # in the next line
  separate_lms(gene, gout_genes)
}

# non-parametric test for neutrophil count
spearman15 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000198074, method = "spearman", exact=FALSE)
print(spearman15)
spearman16 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000163220, method = "spearman", exact=FALSE)
print(spearman16)
spearman17 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000115919, method = "spearman", exact=FALSE)
print(spearman17)
spearman18 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000206073, method = "spearman", exact=FALSE)
print(spearman18)

# Sex nor neutrophil count seem to predict the gene expression values for the top 2 genes found to be differentially 
# expressed in sepsis versus gout

for (gene in top_genes2) {
  cat("Running models for gene:", gene, "\n") # this is to include a wee title for the gene and show the results 
  # in the next line
  separate_lms(gene, sa_genes)
}

# non-parametric test for neutrophil count
spearman19 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000198074, method = "spearman", exact=FALSE)
print(spearman19)
spearman20 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000163220, method = "spearman", exact=FALSE)
print(spearman20)
spearman21 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000115919, method = "spearman", exact=FALSE)
print(spearman21)
spearman22 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000206073, method = "spearman", exact=FALSE)
print(spearman22)

# sex nor neutrphils seem to predict the gene expression values for the top 4 genes found to be differentially
# expressed in sepsis versus gout

# and finally we repeat this for the last 4 genes in the top 10 genes from the DESeq2 analysis
# "ENSG00000134827" "ENSG00000135373" "ENSG00000241794" "ENSG00000153802"

gout_merge <- gout_samples %>% dplyr::select(ENSG00000134827, ENSG00000135373, ENSG00000241794, ENSG00000153802)
gout_genes <- cbind(gout_genes, gout_merge)

sa_merge <- sa_samples %>% dplyr::select(ENSG00000134827, ENSG00000135373, ENSG00000241794, ENSG00000153802)
sa_genes <- cbind(sa_genes, sa_merge)

# checking if there is any repeated gene column in either dataframe gout_genes or sa_genes
any(duplicated(colnames(gout_genes)))
any(duplicated(colnames(sa_genes)))

# all good!

# we need to make a list with the next 4 genes in the top 10 genes from the DESeq2 analysis
top10_genes
# "ENSG00000134827" "ENSG00000135373" "ENSG00000241794" "ENSG00000153802"
top_genes3 <- significant_genes$ID[7:10]
top_genes3 # checking that they were selected correctly

# now we apply the function to these genes for gout and sa dataframes
for (gene in top_genes3) {
  cat("Running models for gene:", gene, "\n") # this is to include a wee title for the gene and show the results 
  # in the next line
  separate_lms(gene, gout_genes)
}

# non-parametric test for neutrophil count
spearman23 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000134827, method = "spearman", exact=FALSE)
print(spearman23)
spearman24 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000135373, method = "spearman", exact=FALSE)
print(spearman24)
spearman25 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000241794, method = "spearman", exact=FALSE)
print(spearman25)
spearman26 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000153802, method = "spearman", exact=FALSE)
print(spearman26)

# sex nor neutrophil count seem to predict the gene expression values for the top 2 genes found to be differentially 
# expressed in sepsis versus gout

for (gene in top_genes3) {
  cat("Running models for gene:", gene, "\n") # this is to include a wee title for the gene and show the results 
  # in the next line
  separate_lms(gene, sa_genes)
}

# non-parametric test for neutrophil count
spearman27 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000134827, method = "spearman", exact=FALSE)
print(spearman27)
spearman28 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000135373, method = "spearman", exact=FALSE)
print(spearman28)
spearman29 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000241794, method = "spearman", exact=FALSE)
print(spearman29)
spearman30 <- cor.test(gout_genes$NEUTROPHILS, gout_genes$ENSG00000153802, method = "spearman", exact=FALSE)
print(spearman30)

# sex nor neutrphils seem to predict the gene expression values for the top 4 genes found to be differentially
# expressed in sepsis versus gout