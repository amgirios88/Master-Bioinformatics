if (!requireNamespace("amap", quietly = TRUE)) {
  install.packages("amap")
}
BiocManager::install("org.Mm.eg.db")


library(readxl)
library(dplyr)
library(ggplot2)
# in order to use melt() we need to load the reshape2 package
library(reshape2)
library(ggrepel)
library(RColorBrewer)
library(amap)
library(GO.db)
library(org.Hs.eg.db)
library(clusterProfiler)

# load data ----------------------------------------------------------------

Human_Background_GRCh38.p13 <- read.table("E:/Uni/MASTER/Classes/Data exploration R/Human_Background_GRCh38.p13.txt", 
                                          header=TRUE, quote="", comment.char="")

# 59421 obs and 6 variables: X.ID, SYMBOL, CHROMOSOME, START, STOP, BIOTYPE.
# cleaning tasks: X.ID to ID, BIOTYPE. to BIOTYPE, data in ID clean first ", data in BIOTYPE clean the last "

EM <- read.table("E:/Uni/MASTER/Classes/Data exploration R/EM.txt", header=TRUE, quote="", comment.char="")

# same obs, 10 variables: X.ID, Prolif_1, Prolif_2, Prolif_3, Senes_1, Senes_2, Senes_3, Senes_MtD_1, Senes_MtD_2, Senes_MtD_3.
# cleaning tasks: X.ID to ID, Senes_MtD_3.to Senes_MtD_3, data in ID clean first ", data in Senes_MtD_3 clean the last "

sample_sheet <- read.table("E:/Uni/MASTER/Classes/Data exploration R/sample_sheet.txt", header=TRUE, quote="", comment.char="")

# 9 obs, 2 variables: X.SAMPLE, SAMPLE_GROUP.
# cleaning tasks: X.SAMPLE to SAMPLE, SAMPLE_GROUP. to SAMPLE_GROUP, data in SAMPLE clean first ", data in SAMPLE_GROUP clean the last "

DE_Senes_vs_Prolif <- read.table("E:/Uni/MASTER/Classes/Data exploration R/DE_Senes_vs_Prolif.txt", header=TRUE, quote="", comment.char="")

DE_Senes_MtD_vs_Senes <- read.table("E:/Uni/MASTER/Classes/Data exploration R/DE_Senes_MtD_vs_Senes.txt", header=TRUE, quote="", comment.char="")

DE_Senes_MtD_vs_Prolif <- read.table("E:/Uni/MASTER/Classes/Data exploration R/DE_Senes_MtD_vs_Prolif.txt", header=TRUE, quote="", comment.char="")

# for all three exp tables, 4 variables, X.ID, log2fold, p, p.adj.
# cleaning tasks: X.ID to ID, p.adj. to p.adj, data in ID clean first ", data in p.adj clean the last "
# convert p.adj to num




# Data cleaning -----------------------------------------------------------

# Background table
colnames(Human_Background_GRCh38.p13)[1] <- "ID"
colnames(Human_Background_GRCh38.p13)[6] <- "BIOTYPE"
# now we remove the first " in the string data of ID column
Human_Background_GRCh38.p13$ID <- gsub("\"", "", Human_Background_GRCh38.p13$ID) 
# now we remove the last " in the string data of BIOTYPE column
Human_Background_GRCh38.p13$BIOTYPE <- gsub("\"", "", Human_Background_GRCh38.p13$BIOTYPE)
# we change the name to something shorter without creating a new dataframe
humanbg <- Human_Background_GRCh38.p13
rm(Human_Background_GRCh38.p13)

# EM table

colnames(EM)[1] <- "ID"
colnames(EM)[10] <- "Senes_MtD_3"
# now we remove the first " in the string data of ID column
EM$ID <- gsub("\"", "", EM$ID)
# now we remove the last " in the string data of Senes_MtD_3 column
EM$Senes_MtD_3 <- gsub("\"", "", EM$Senes_MtD_3)
# convert Senes_MtD_3 to numeric
EM$Senes_MtD_3 <- as.numeric(EM$Senes_MtD_3)

# Sample_sheet table

colnames(sample_sheet)[1] <- "SAMPLE"
colnames(sample_sheet)[2] <- "SAMPLE_GROUP"
# now we remove the first " in the string data of SAMPLE column
sample_sheet$SAMPLE <- gsub("\"", "", sample_sheet$SAMPLE)
# now we remove the last " in the string data of SAMPLE_GROUP column
sample_sheet$SAMPLE_GROUP <- gsub("\"", "", sample_sheet$SAMPLE_GROUP)
# we rename the dataframe as sample
sample <- sample_sheet
rm(sample_sheet)

# DE_Senes_vs_Prolif table

colnames(DE_Senes_vs_Prolif)[1] <- "ID"
colnames(DE_Senes_vs_Prolif)[4] <- "p.adj"
# now we remove the first " in the string data of ID column
DE_Senes_vs_Prolif$ID <- gsub("\"", "", DE_Senes_vs_Prolif$ID)
# now we remove the last " in the string data of p.adj column
DE_Senes_vs_Prolif$p.adj <- gsub("\"", "", DE_Senes_vs_Prolif$p.adj)
# convert p.adj to numeric
DE_Senes_vs_Prolif$p.adj <- as.numeric(DE_Senes_vs_Prolif$p.adj)
# we rename the dataframe
DE_SvP <- DE_Senes_vs_Prolif
rm(DE_Senes_vs_Prolif)

# DE_Senes_MtD_vs_Senes table

colnames(DE_Senes_MtD_vs_Senes)[1] <- "ID"
colnames(DE_Senes_MtD_vs_Senes)[4] <- "p.adj"
# now we remove the first " in the string data of ID column
DE_Senes_MtD_vs_Senes$ID <- gsub("\"", "", DE_Senes_MtD_vs_Senes$ID)
# now we remove the last " in the string data of p.adj column
DE_Senes_MtD_vs_Senes$p.adj <- gsub("\"", "", DE_Senes_MtD_vs_Senes$p.adj)
# convert p.adj to numeric
DE_Senes_MtD_vs_Senes$p.adj <- as.numeric(DE_Senes_MtD_vs_Senes$p.adj)
# we rename the dataframe
DE_SMvS <- DE_Senes_MtD_vs_Senes
rm(DE_Senes_MtD_vs_Senes)


# DE_Senes_MtD_vs_Prolif table

colnames(DE_Senes_MtD_vs_Prolif)[1] <- "ID"
colnames(DE_Senes_MtD_vs_Prolif)[4] <- "p.adj"
# now we remove the first " in the string data of ID column
DE_Senes_MtD_vs_Prolif$ID <- gsub("\"", "", DE_Senes_MtD_vs_Prolif$ID)
# now we remove the last " in the string data of p.adj column
DE_Senes_MtD_vs_Prolif$p.adj <- gsub("\"", "", DE_Senes_MtD_vs_Prolif$p.adj)
# convert p.adj to numeric
DE_Senes_MtD_vs_Prolif$p.adj <- as.numeric(DE_Senes_MtD_vs_Prolif$p.adj)
# we rename the dataframe
DE_SMvP <- DE_Senes_MtD_vs_Prolif
rm(DE_Senes_MtD_vs_Prolif)




# Data exploration --------------------------------------------------------

# Background table

summary(humanbg)
# data types are correct

# EM table

summary(EM)

# Prolif_1: median is 49.06, mean is 1338.62, this means that the data is skewed to the right
# Prolif_2: median is 48.80, mean is 1353.40, this means that the data is skewed to the right
# Prolif_3: median is 49.91, mean is 1324.03, this means that the data is skewed to the right

# Senes_1: median is 59.1, mean is 1372.7, this means that the data is skewed to the right
# Senes_2: median is 60.6, mean is 1441.6, this means that the data is skewed to the right
# Senes_3: median is 59.1, mean is 1416.6, this means that the data is skewed to the right

# Senes_MtD_1: median is 51.05, mean is 1323.44, this means that the data is skewed to the right
# Senes_MtD_2: median is 50.88, mean is 1311.10, this means that the data is skewed to the right
# Senes_MtD_3: median is 50.84, mean is 1305.17, this means that the data is skewed to the right

#* *DE_Senes_vs_Prolif table*

summary(DE_SvP)
# minimum of -7.8071, maximum of 8.2959, median 0.1275, mean 0.3291

# histogram
hist(DE_SvP$log2fold, breaks=50, col="blue", main="Histogram of log2fold in DE_Senes_vs_Prolif", 
     xlab="log2fold")

#* *DE_Senes_MtD_vs_Senes table*

summary(DE_SMvS)
# minimum of -10.81090, maximum of 10.00047, median -0.08625, mean -0.22888

# histogram
hist(DE_SMvS$log2fold, breaks=50, col="blue", main="Histogram of log2fold in DE_MtD_vs_Senes", 
     xlab="log2fold")

#* *DE_Senes_MtD_vs_Prolif table*

summary(DE_SMvP)
# minimum of -7.8561, maximum of 11.3684, median 0.0448, mean 0.1004

# histogram
hist(DE_SMvP$log2fold, breaks=50, col="blue", main="Histogram of log2fold in DE_MtD_vs_Prolif", 
     xlab="log2fold")


#* *Add new column to DE tables for gene names*

# add the new column of gene names from humanbg into DE_SvP by matching ID columns in both tables
# DE_SvP$SYMBOL <- humanbg$SYMBOL[match(DE_SvP$ID, humanbg$ID)]
# # and move it to after the ID column
# DE_SvP <- DE_SvP[, c(1, 5, 2, 3, 4)]
# 
# # do the others
# DE_SMvS$SYMBOL <- humanbg$SYMBOL[match(DE_SMvS$ID, humanbg$ID)]
# DE_SMvS <- DE_SMvS[, c(1, 5, 2, 3, 4)]
# 
# DE_SMvP$SYMBOL <- humanbg$SYMBOL[match(DE_SMvP$ID, humanbg$ID)]
# DE_SMvP <- DE_SMvP[, c(1, 5, 2, 3, 4)]


#* *Merge DE tables with background table to have gene names and other relevant data*

# we merge background dataframe with DEs
# I want to only merge the SYMBOL column from humanbg into the DE tables

DE_SvP_bg <- merge(DE_SvP, humanbg, by = "ID")[c("SYMBOL", "ID", "log2fold", "p.adj")]

DE_SMvS_bg <- merge(DE_SMvS, humanbg, by = "ID")[c("SYMBOL", "ID", "log2fold", "p.adj")]

DE_SMvP_bg <- merge(DE_SMvP, humanbg, by = "ID")[c("SYMBOL", "ID", "log2fold", "p.adj")]


#* *Filter the DE tables for significant genes*

DE_SvP_bg <- subset(DE_SvP_bg, p.adj < 0.05)

DE_SMvS_bg <- subset(DE_SMvS_bg, p.adj < 0.05)

DE_SMvP_bg <- subset(DE_SMvP_bg, p.adj < 0.05)




# Density plots -----------------------------------------------------------

# let's add the SYMBOL column from humanbg table into the EM table by ID

EM$SYMBOL <- humanbg$SYMBOL[match(EM$ID, humanbg$ID)]

# we need to melt the data into one column for values and another for sample

EM_melted <- melt(EM)

# now we can plot the density of gene expression for each group

dens_plot <- ggplot(EM_melted, aes(x = log10(value + 0.01), fill = variable)) +
  geom_density(alpha = 0.75) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Density plot of gene expression for each group",
       x = "Gene expression",
       y = "Density") +
  theme_minimal()

dens_plot

# the distribution of gene expression values for each group is shown in the density plot
# they look quite homogeneous and have bimodal distributions




# Principal components analysis ---------------------------------------------------------------------

#* *Exploring expression profile for each group*

# we need to transpose the expression dataframe and set gene names as column headers so we can merge this information
# with the clinical dataframe
# now we merge transposed EM with background data

EM_merged <- merge(EM, humanbg, by.x = "ID", by.y = "ID", all = TRUE)


# are there any duplicate genes SYMBOL in the EM_merged dataframe?

any(duplicated(EM_merged$SYMBOL))

# this won't allow us to set SYMBOL names to work with them, so we will work with the ID names and
# assign gene names associated at the end to show the results

# we can transpose now so it will be easier to work on the exp data
EM_transposed <- as.data.frame(t(EM_merged[, -1])) # we remove the ID column as we don't need it
colnames(EM_transposed) <- EM$ID # we set the gene names as column headers

# we clean genes with missing values

EM_transposed <- EM_transposed[, colSums(is.na(EM_transposed)) == 0]

# we will move the gene symbol row to the first row position (not row.names) using base r

new_order <- c(10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14)
EM_transposed <- EM_transposed[new_order, ]

# we will prepare the data removing non-numeric columns (SAMPLE, SAMPLE_GROUP, SEX, NEUTROPHILS)
EM_pca <- EM_transposed[-c(1, 11, 12, 13, 14), ]

# standardize the data (important step before PCA to make the data comparable - PCA seeks to maximize the variance 
# of each component; a PCA without normalization would perform worse than one with normalization)
# note to myself about why is important to standardize the data before PCA: 
# https://stats.stackexchange.com/questions/69157/why-do-we-need-to-normalize-data-before-principal-component-analysis-pca

# convert all rows to num

EM_pca <- as.data.frame(apply(EM_pca, 2, function(x) as.numeric(as.character(x))))

# we standardize the data

EM_pca_scaled <- scale(EM_pca)

# the PCA will be performed on the samples (rows) and we will see where these patients fit according to their gene
# expression profile
# we perform PCA on the standardized data
pca_result <- prcomp(EM_pca_scaled, center = TRUE, scale. = TRUE)

# this is the PCA summary to check variance explained by each component
summary(pca_result)

# to plot the data, we can extract PCA coordinates for the samples (rows)
pca_data <- as.data.frame(pca_result$x)

pca_data$SAMPLE_GROUP <- sample$SAMPLE_GROUP[match(rownames(pca_data), rownames(sample))]

# we have the sample group information in EM_pca_samp

# finally we plot the first two principal components
pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = SAMPLE_GROUP)) +
  geom_point(size = 3) + # this line of code says that we want to plot points with a size of 3
  labs(title = "PCA of Gene Expression for Senes with/out mitochondria and prolif", # we can add some labels to the plot
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()
pca

# each point corresponds to a sample (either Prolif, Senes, or senes_MtD) and shows how that sample's 
# gene expression contributes to the principal components

# group prolif does not play any role in explaining the variation on PC1, but it explains most of PC2
# groups senes and senes MtD explain most of PC1 and some of PC2, although not as much as prolif

# we need to get what percentage of each component PC1 and PC2 is explained by data points in the pca plot
# our dataframe for the pca analysis is pca_data

pca_variance <- pca_result$sdev^2
pca_var_percentage <- round(pca_variance / sum(pca_variance) * 100, 1)


pca <- pca + xlab(paste0("PC1 - ", pca_var_percentage[1], "%"))
pca <- pca + ylab(paste0("PC2 - ", pca_var_percentage[2], "%"))
pca

# different skin types explain 38.7% of variability in PC1, and 32.5% of variability in PC2

# we can save the plot
ggsave("E:/Uni/MASTER/Classes/Data exploration R/assessment results/pca_plot.png", pca, width = 8, height = 6, dpi = 300)




# Checking similarity of differentially expressed genes between the DE data --------

#* *Comparison between Proliferating vs MtD cells and Senescent vs MtD cells* 
# we will use the dataframes we subset previously for significant genes (adjusted p-value < 0.05)
# we extract count of genes in each, which will be total significant differentially expressed genes in each

total_SMvP <- length(DE_SMvP_bg$ID)
total_SMvS <- length(DE_SMvS_bg$ID)

# now we check the overlap and uniqueness between the three DE data (significant)

# between MtD vs Prolif and MtD vs Senes will tell us what genes are differentially expressed in senescent cells
# and proliferating cells from the mitochondria depleted cells; chekcing the uniqueness of differentially expressed
# genes will tell us which genes are differentially expressed specific for senescent cells and proliferating cells
common_genes_SP <- intersect(DE_SMvS_bg$ID, DE_SMvP_bg$ID)
common_length_SP <- length(common_genes_SP)

# now we get the length of unique genes (those in proliferating cells but not in senescent ones, and vice versa)
unique_SMvP <- total_SMvP - common_length_SP
unique_SMvS <- total_SMvS - common_length_SP

# and we create a dataframe with the data for the significant differentially expressed genes (sdge)
# between mitochondria depleted cells that are not found in senescent and proliferating cells
summary_sdgePS <- data.frame(
  Category = c("MtD vs Proliferating", "MtD vs Senescent", "Overlap"),
  Total_Significant_Genes = c(total_SMvP, total_SMvS, common_length_SP),
  Unique_Genes = c(unique_SMvP, unique_SMvS, NA)
)

print(summary_sdgePS)

# we can save the summary to a file
write.csv(summary_sdgePS, "E:/Uni/MASTER/Classes/Data exploration R/assessment results/summary_sdgePS.csv", row.names = FALSE)

# 8697 overlapping genes are found differentially expressed in senescent mitochondria depleted cells relative to both senescent and proliferating cells
# these genes may be related to functions that depend on the presence of the mitochondria, since its expression is differential in senescent cells
# without mitochondria when compared to either proliferating or senescent cells
# we can take a look at what these genes are

# but we don't want all of them, only the 10 most differentially expressed in absolute terms in both cases

common_genes_SMvS <- subset(DE_SMvS_bg, ID %in% common_genes_SP)
common_genes_SMvP <- subset(DE_SMvP_bg, ID %in% common_genes_SP)

# now we want to reorder the two dataframes by most DE in absolute terms

common_genes_SMvS <- common_genes_SMvS[order(-abs(common_genes_SMvS$log2fold)), ]
common_genes_SMvP <- common_genes_SMvP[order(-abs(common_genes_SMvP$log2fold)), ]

# and retrieve the first 10 genes that are most expressed in both dataframes

# add a rank column to each dataframe
common_genes_SMvS$Rank_SMvS <- seq_len(nrow(common_genes_SMvS))
common_genes_SMvP$Rank_SMvP <- seq_len(nrow(common_genes_SMvP))

# merge the two dataframes on the common gene IDs
merged_common_genes <- merge(common_genes_SMvS, common_genes_SMvP, by = "ID", suffixes = c("_SMvS", "_SMvP"))

# add a concatenated rank column
merged_common_genes$Rank_Concat <- paste(merged_common_genes$Rank_SMvS, merged_common_genes$Rank_SMvP, sep = "-")

# split Rank_Concat into two separate columns
rank_split <- do.call(rbind, strsplit(merged_common_genes$Rank_Concat, "-", fixed = TRUE))

# convert to a numeric matrix
rank_split_numeric <- apply(rank_split, 2, as.numeric)

# calculate the sum of the two ranks (meaning in the concatenated variable x-y we sum x and y)
merged_common_genes$Rank_Sum <- rowSums(rank_split_numeric)

# sort the dataframe by Rank_Sum in ascending order (so the lower sum appears first)
merged_common_genes <- merged_common_genes[order(merged_common_genes$Rank_Sum), ]

# what are the 10 genes that are most DE in both cases?
print(merged_common_genes[1:10, ])

# Save the sorted dataframe to a file if needed
write.csv(merged_common_genes, "E:/Uni/MASTER/Classes/Data exploration R/assessment results/sorted_merged_common_genes.csv", row.names = FALSE)

# we print only the symbols for the 10 most differentially expressed genes in both cases
# since these are overlapping genes, just by retrieving the symbols from one of the dataframes is enough

print(merged_common_genes[1:10, "SYMBOL_SMvS"])

# "OR7E19P", "EGLN3", "CXCR4", "CARMN", "VCAM1", "MT3", "CPXM1", "AK4", "MTATP6P2", "MTND4LP1"




#* *Comparison between Proliferating vs Senescent cells*
# we can now proceed to explore which genes are differentially expressed directly between senescent and proliferating cells

# DE data in this case is for differentially expressed senescent cells in comparison to proliferating cells 
# therefore the significant genes we obtain here will be the ones that are expressed differently for 
# senescent cells

# we need to calculate the difference between DE_SvP and DE_SvP_bg

unique_S <- setdiff(DE_SvP$ID, DE_SvP_bg$ID)
unique_S <- length(unique_S)
print(unique_S)

# there are 13071 genes that are not differentially expressed between senescent and proliferating cells

print(length(DE_SvP_bg$ID))

# and 11787 genes are differentially expressed in senescent cells in comparison to proliferating cells
# so these genes could potentially be used to differentiate senescent cells from proliferating cells 
# and so, potentially to tell if a cell is proliferating or senescent

# we can now reorder the DE_SvP_bg dataframe in ascending order by p.adj
DE_SvP_bg <- DE_SvP_bg[order(DE_SvP_bg$p.adj), ]

# now reorder DE_SvP_bg by p.adj first in ascending order and then by log2fold in descending order in absolute terms

DE_SvP_bg <- DE_SvP_bg[order(DE_SvP_bg$p.adj, -abs(DE_SvP_bg$log2fold)), ]

# and print the first 10 genes
print(DE_SvP_bg[1:10, ])
# these are the 10 genes more significantly and most differentially expressed in senescent cells in comparison
# with proliferating cells:
#     -Positive log2 fold changes will indicate genes that are upregulated in senescent cells relative to proliferating cells
#     -Negative log2 fold changes will indicate genes that are downregulated in senescent cells relative to proliferating cells

# we can get a list with the SYMBOL names of these genes

print(DE_SvP_bg[1:10, "SYMBOL"])
# "ARRDC4", "TXNIP", "OLR1", "HAS2", "KIF20A", "KIF4A", "BIRC5", "PRC1", "FIBIN", "PODXL" 

# these could potentially be used to differentiate proliferating cells from senescent cells

# we can save the list of genes to a file
write.csv(DE_SvP_bg[1:10,], "E:/Uni/MASTER/Classes/Data exploration R/assessment results/10DE_prolif.csv", row.names = FALSE)
# we could now manually explore more these genes to see if they can be used to differentiate proliferating cells from senescent cells




# MA plot -----------------------------------------
#* *MA PLOT*
# https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html (MA plot)

# first we get the main table for this 
# it needs to contain the mean expression and log2fold change for each gene

# we will add a mean expression for each gene and group into the EM dataframe
# this menans three means, one per group, we have three groups, so meanP for columns 2 to 4, meanS for 
# columns 5 to 7, and meanSMvD for columns 8 to 10

EM$meanSvP <- rowMeans(EM[, 2:7])
EM$meanSMvP <- rowMeans(EM[, c(2:4, 8:10)])
EM$meanSMvS <- rowMeans(EM[, 5:10])

# we will make three new dataframe with the information for ID + SYMBOL genes and mean expression for each
# DE table, and then we will add the log2fold data from each DE table

MA_SvP_df <- merge(DE_SvP, EM, by = "ID")
MA_SMvP_df <- merge(DE_SMvP, EM, by = "ID")
MA_SMvS_df <- merge(DE_SMvS, EM, by = "ID")

# now we need to subset for significant genes downregulated and significant genes upregulated
# we will use the p.adj column to subset for significant genes
MA_SvP_sig_up <- subset(MA_SvP_df, p.adj < 0.05 & log2fold > 1)
MA_SvP_sig_down <- subset(MA_SvP_df, p.adj < 0.05 & log2fold < -1)
MA_SMvP_sig_up <- subset(MA_SMvP_df, p.adj < 0.05 & log2fold > 1)
MA_SMvP_sig_down <- subset(MA_SMvP_df, p.adj < 0.05 & log2fold < -1)
MA_SMvS_sig_up <- subset(MA_SMvS_df, p.adj < 0.05 & log2fold > 1)
MA_SMvS_sig_down <- subset(MA_SMvS_df, p.adj < 0.05 & log2fold < -1)

# now the top 5 for each DE

# first reorder dataframes by log2fold in descending absolute order for each subset
MA_SvP_sig_up <- MA_SvP_sig_up[order(-abs(MA_SvP_sig_up$log2fold)), ]
MA_SvP_sig_down <- MA_SvP_sig_down[order(-abs(MA_SvP_sig_down$log2fold)), ]

# I printetd first the 5 most upregulated genes and then the 5 most downregulated genes
# but I could not plot them all without the labels overlapping, so I decided to do
# them separately according to the gene names that needed more adjustment
MA_SvP_sig_up_top5_1 = subset(MA_SvP_sig_up, SYMBOL %in% c("ADAMTS19", "PDPN"))
MA_SvP_sig_up_top5_2 = subset(MA_SvP_sig_up, SYMBOL %in% c("ITIH5", "C3", "FAM238C"))
MA_SvP_sig_down_top5_1 = subset(MA_SvP_sig_down, SYMBOL %in% c("DEPDC1-AS1"))
MA_SvP_sig_down_top5_2 = subset(MA_SvP_sig_down, SYMBOL %in% c("GABRA4", "WNT7B", "LCT-AS1", "NEK2"))


MA_SvP = ggplot(MA_SvP_df, aes(x=log10(meanSvP), y=log2fold)) + 
  # adds the dots
  geom_point(colour = "black") + 
  geom_point(data = MA_SvP_sig_down, colour = "blue") + 
  geom_point(data = MA_SvP_sig_up, colour = "red") + 
  # add the dashed lines
  geom_hline(yintercept=-1,linetype="dashed") + 
  geom_hline(yintercept=1,linetype="dashed") + 
  # add axis labels and the title
  labs(title = "MA plot for senescent vs proliferating cells", x= "Mean Expression (log10)", y= "Log2 fold change") +
 # add data labels
  geom_label_repel(data = MA_SvP_sig_down_top5_1, 
                   aes(label = SYMBOL), 
                   colour = "darkblue",
                   fill = NA,
                   max.overlaps = 5, 
                   size = 3, 
                   segment.color = "black", 
                   segment.size = 0.5,
                   nudge_y = -0.5,
                   nudge_x = -0.1,
                   box.padding = 0.5,
                   point.padding = 0.4 ) + 
  geom_label_repel(data = MA_SvP_sig_up_top5_1, 
                   aes(label = SYMBOL), 
                   colour = "darkred",
                   fill = NA,
                   max.overlaps = 3, 
                   size = 3, 
                   segment.color = "black", 
                   segment.size = 0.4,
                   nudge_y = -0.1,
                   nudge_x = -0.1,
                   box.padding = 0.9,
                   point.padding = 0.1) +
  geom_label_repel(data = MA_SvP_sig_down_top5_2, 
                 aes(label = SYMBOL), 
                 colour = "darkblue",
                 fill = NA,
                 max.overlaps = 5, 
                 size = 3, 
                 segment.color = "black", 
                 segment.size = 0.5,
                 nudge_y = -0.2,
                 nudge_x = 0.4,
                 box.padding = 0.5,
                 point.padding = 0.3 ) +
  geom_label_repel(data = MA_SvP_sig_up_top5_2, 
                   aes(label = SYMBOL), 
                   colour = "darkred",
                   fill = NA,
                   max.overlaps = 3, 
                   size = 3, 
                   segment.color = "black", 
                   segment.size = 0.5,
                   nudge_y = 0.7,
                   nudge_x = 0,
                   box.padding = 0.5,
                   point.padding = 0.4)
  
MA_SvP


MA_SMvP_sig_up <- MA_SMvP_sig_up[order(-abs(MA_SMvP_sig_up$log2fold)), ]
MA_SMvP_sig_down <- MA_SMvP_sig_down[order(-abs(MA_SMvP_sig_down$log2fold)), ]
MA_SMvP_sig_up_top5_1 = subset(MA_SMvP_sig_up, SYMBOL %in% c("CXCR4", "OR7E19P"))
MA_SMvP_sig_up_top5_2 = subset(MA_SMvP_sig_up, SYMBOL %in% c("LINC01164", "MRGPRX3", "EGLN3"))
MA_SMvP_sig_down_top5_1 = subset(MA_SMvP_sig_down, SYMBOL %in% c("CIDEA", "IL1RL2", "KRT19"))
MA_SMvP_sig_down_top5_2 = subset(MA_SMvP_sig_down, SYMBOL %in% c("CPXM1", "PPP1R14A"))

MA_SMvP_sig_up_top5 = MA_SMvP_sig_up [1:5,]
MA_SMvP_sig_down_top5 = MA_SMvP_sig_down [1:5,]


MA_SMvP = ggplot(MA_SMvP_df, aes(x=log10(meanSMvP), y=log2fold)) + 
  # adds the dots
  geom_point(colour = "black") + 
  geom_point(data = MA_SMvP_sig_down, colour = "blue") + 
  geom_point(data = MA_SMvP_sig_up, colour = "red") + 
  # add the dashed lines
  geom_hline(yintercept=-1,linetype="dashed") + 
  geom_hline(yintercept=1,linetype="dashed") + 
  # add axis labels and the title
  labs(title = "MA plot for senescent mitochondria depleted vs proliferating cells", x= "Mean Expression (log10)", y= "Log2 fold change") +
  # add data labels
  geom_label_repel(data = MA_SMvP_sig_down_top5_1,
                   aes(label = SYMBOL),
                   colour = "darkblue",
                   fill = NA,
                   max.overlaps = 5,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = -0.3,
                   nudge_x = -0.4,
                   box.padding = 0.5,
                   point.padding = 0.1 ) +
  geom_label_repel(data = MA_SMvP_sig_down_top5_2,
                   aes(label = SYMBOL),
                   colour = "darkblue",
                   fill = NA,
                   max.overlaps = 5,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = -0.2,
                   nudge_x = 0.3,
                   box.padding = 0.5,
                   point.padding = 0.1) +
  geom_label_repel(data = MA_SMvP_sig_up_top5_1,
                   aes(label = SYMBOL),
                   colour = "darkred",
                   fill = NA,
                   max.overlaps = 3,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = -0.2,
                   nudge_x = -0.4,
                   box.padding = 0.3,
                   point.padding = 0.1) +
  geom_label_repel(data = MA_SMvP_sig_up_top5_2,
                   aes(label = SYMBOL),
                   colour = "darkred",
                   fill = NA,
                   max.overlaps = 3,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = -0.7,
                   nudge_x = 0.5,
                   box.padding = 0.3,
                   point.padding = 0.1)

MA_SMvP


MA_SMvS_sig_up <- MA_SMvS_sig_up[order(-abs(MA_SMvS_sig_up$log2fold)), ]
MA_SMvS_sig_down <- MA_SMvS_sig_down[order(-abs(MA_SMvS_sig_down$log2fold)), ]
MA_SMvS_sig_up_top5_1 = subset(MA_SMvS_sig_up, SYMBOL %in% c("CXCR4", "OR7E19P", "EGLN3"))
MA_SMvS_sig_up_top5_2 = subset(MA_SMvS_sig_up, SYMBOL %in% c("ENSG00000233818", "INHBB"))
MA_SMvS_sig_down_top5_1 = subset(MA_SMvS_sig_down, SYMBOL %in% c("VCAM1", "CXCL10","MTATP6P2"))
MA_SMvS_sig_down_top5_2 = subset(MA_SMvS_sig_down, SYMBOL %in% c("CXCL6"))
MA_SMvS_sig_down_top5_3 = subset(MA_SMvS_sig_down, SYMBOL %in% c("MTCO3P12"))



MA_SMvS = ggplot(MA_SMvS_df, aes(x=log10(meanSMvS), y=log2fold)) + 
  # adds the dots
  geom_point(colour = "black") + 
  geom_point(data = MA_SMvS_sig_down, colour = "blue") + 
  geom_point(data = MA_SMvS_sig_up, colour = "red") + 
  # add the dashed lines
  geom_hline(yintercept=-1,linetype="dashed") + 
  geom_hline(yintercept=1,linetype="dashed") + 
  # add axis labels and the title
  labs(title = "MA plot for senescent mitochondria depleted vs senescent cells", x= "Mean Expression (log10)", y= "Log2 fold change") +
  # add data labels
  geom_label_repel(data = MA_SMvS_sig_down_top5_1,
                   aes(label = SYMBOL),
                   colour = "darkblue",
                   fill = NA,
                   max.overlaps = 5,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = -0.3,
                   nudge_x = 0.3,
                   box.padding = 0.5,
                   point.padding = 0.1 ) +
  geom_label_repel(data = MA_SMvS_sig_down_top5_2,
                   aes(label = SYMBOL),
                   colour = "darkblue",
                   fill = NA,
                   max.overlaps = 5,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = -0.9,
                   nudge_x = 0.6,
                   box.padding = 0.5,
                   point.padding = 0.1) +
  geom_label_repel(data = MA_SMvS_sig_down_top5_3,
                   aes(label = SYMBOL),
                   colour = "darkblue",
                   fill = NA,
                   max.overlaps = 5,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = 0.1,
                   nudge_x = -0.8,
                   box.padding = 0.5,
                   point.padding = 0.1) +
  geom_label_repel(data = MA_SMvS_sig_up_top5_1,
                   aes(label = SYMBOL),
                   colour = "darkred",
                   fill = NA,
                   max.overlaps = 3,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = 0.2,
                   nudge_x = -0.5,
                   box.padding = 0.3,
                   point.padding = 0.1) +
  geom_label_repel(data = MA_SMvS_sig_up_top5_2,
                   aes(label = SYMBOL),
                   colour = "darkred",
                   fill = NA,
                   max.overlaps = 3,
                   size = 3,
                   segment.color = "black",
                   segment.size = 0.5,
                   nudge_y = -0.2,
                   nudge_x = -0.9,
                   box.padding = 0.3,
                   point.padding = 0.1)

MA_SMvS




# Boxplots and violin plots for specific genes -----------------------------------------------

# we can do a violin plot for each sample and sample group in EM dataframe

# we will make a dataframe with the sample groups in EM, each gene, and the expression of each gene in each sample
# first transpose EM (we already have a transposed version EM_transposed, we clean it from other variables)

# we make SYMBOL.x the column names

colnames(EM_transposed) <- EM_transposed["SYMBOL.x", ]
# and we delent SYMBOL row and rows that we do not need
EM_transposed_SM <- EM_transposed[-c(1, 11, 12, 13, 14), ]

# we will add a column for the sample group in EM_transposed, so 1:3 are Prolif, 4:6 are Senes, and 7:9 are SenesMtD

EM_transposed_SM$SAMPLE_GROUP <- rep(c("Prolif", "Senes", "SenesMtD"), each = 3)

# I want to move this new column SAMPLE_GROUP to the first column position using relocate

EM_transposed_SM <- EM_transposed_SM %>% relocate(SAMPLE_GROUP, .before = 1)

EM_melted_SM <- melt(EM_transposed_SM, id.vars = "SAMPLE_GROUP")

# we will use this dataframe to plot facets boxplots and violin plots


#* *Senes MtD vs (Senes/Prolif)*
# genes commonly expressed differently in senescent cells with mitochondria depletion when compared to either
# prolif or senes with mitochondria

# "OR7E19P", "EGLN3", "CXCR4", "CARMN", "VCAM1", "MT3", "CPXM1", "AK4", "MTATP6P2", "MTND4LP1"
# which are: "ENSG00000087250" "ENSG00000088882" "ENSG00000121966" "ENSG00000129521" "ENSG00000162433" 
# "ENSG00000162692" "ENSG00000225980" "ENSG00000231501" "ENSG00000249669" "ENSG00000270307"

EM_transposed_SM <- EM_transposed[-c(1, 11, 12, 13, 14), ]

# we felter EM_transposed dataframe for only these genes and convert values to num

EM_box_SM <- subset(EM_melted_SM, variable %in% c("OR7E19P", "EGLN3", "CXCR4", "CARMN", "VCAM1", "MT3", "CPXM1", "AK4", "MTATP6P2", "MTND4LP1"))
EM_box_SM$value <- as.numeric(EM_box_SM$value)

# facets boxplots for these genes and both groups

EM_boxplot_SM <- ggplot(EM_box_SM, aes(x = SAMPLE_GROUP, y = value, fill = SAMPLE_GROUP)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free") +  # Separate facet for each gene
  labs(title = "Gene expression for top 10 differentially expressed genes in senescent 
       mitochondria depleted vs proliferating/senescent cells",
       x = "Sample group",
       y = "Gene expression") +
  theme(axis.text.x = element_text(angle = 25, size = 8),
        plot.title = element_text(size=10, hjust = 0.5, face = "bold"))


EM_boxplot_SM

# we can save the plot

ggsave("E:/Uni/MASTER/Classes/Data exploration R/assessment results/boxplot_SM.png", EM_boxplot_SM, width = 8, height = 6, dpi = 300)


#* *Senes vs Prolif*
# genes that are differentially expressed in senescent cells when compared to proliferating cells
# "ARRDC4", "TXNIP", "OLR1", "HAS2", "KIF20A", "KIF4A", "BIRC5", "PRC1", "FIBIN", "PODXL"

EM_transposed_PS <- EM_transposed[-c(1, 11, 12, 13, 14), ]

EM_transposed_PS$SAMPLE_GROUP <- rep(c("Prolif", "Senes", "SenesMtD"), each = 3)

EM_transposed_PS <- EM_transposed_PS %>% relocate(SAMPLE_GROUP, .before = 1)

EM_melted_PS <- melt(EM_transposed_PS, id.vars = "SAMPLE_GROUP")


EM_box_PS <- subset(EM_melted_PS, variable %in% c("ARRDC4", "TXNIP", "OLR1", "HAS2", "KIF20A", "KIF4A", "BIRC5", "PRC1", "FIBIN", "PODXL"))
EM_box_PS$value <- as.numeric(EM_box_PS$value)

# facets boxplots for these genes and both groups

EM_boxplot_PS <- ggplot(EM_box_PS, aes(x = SAMPLE_GROUP, y = value, fill = SAMPLE_GROUP)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free") +  # Separate facet for each gene
  labs(title = "Gene expression for top 10 differentially expressed genes in senescent 
       cells vs proliferating cells",
       x = "Sample group",
       y = "Gene expression") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
        plot.title = element_text(size=10, hjust = 0.5, face = "bold"))


EM_boxplot_PS




# Heatmaps ----------------------------------------------------------------


EM_heatmap <- EM_transposed[-c(1, 11, 12, 13, 14), ]

# EM_heatmap$SAMPLE_GROUP <- rep(c("Prolif", "Senes", "SenesMtD"), each = 3)
# 
# EM_heatmap <- EM_heatmap %>% relocate(SAMPLE_GROUP, .before = 1)
# 
# #rownames as sample names and delete sample group column
# rownames(EM_heatmap) <- EM_heatmap$SAMPLE_GROUP

heatmap_matrix <- as.matrix(EM_heatmap)
y.dist = Dist(heatmap_matrix, method="spearman")
y.clust = hclust(y.dist, method="average")
y.tree = as.dendrogram(y.clust)
y.tree_reorder = reorder(y.tree,0,FUN="average")
y_order = order.dendrogram(y.tree_reorder)
heatmap_matrix_clustered = heatmap_matrix[y_order,]

colours = c("darkgreen", "lightgreen", "skyblue","pink","darkred", "yellow")
palette = colorRampPalette(colours)(500)

heatmap_matrix_clustered = melt(heatmap_matrix_clustered)
heatmap_matrix_clustered$value <- as.numeric(heatmap_matrix_clustered$value)

range(heatmap_matrix_clustered$value, na.rm = TRUE)
summary(heatmap_matrix_clustered$value)
head(heatmap_matrix_clustered)

heatmap_matrix_clustered$value <- log1p(heatmap_matrix_clustered$value) 

heatmap = ggplot(heatmap_matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = palette) +
  theme_minimal() +
  labs(x = "Genes", y = "Samples", fill = "Log Expression")
heatmap




# Over-Representation Analysis (ORA) --------------------------------------

sig_genes_SMvS <- subset(DE_SMvS_bg, abs(log2fold) > 1)
sig_genes_SMvP <- subset(DE_SMvP_bg, abs(log2fold) > 1)
sig_genes_SvP <- subset(DE_SvP_bg, abs(log2fold) > 1) 

# convert IDs to Entrez for cluster profiles
sig_SMvS_entrez = bitr(sig_genes_SMvS, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
sig_SMvP_entrez = bitr(sig_genes_SMvP, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
sig_SvP_entrez = bitr(sig_genes_SvP, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)


