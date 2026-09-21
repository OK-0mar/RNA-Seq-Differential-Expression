############################
# Title: DESeq2 Differential Expression Analysis
# Script by : Omar Khaled
############################


# 1. Calling DESeq2 
library(DESeq2)


# 2. Loading Data

data <- as.matrix(read.csv("bladder_counts.csv", row.names = 1))
pheno <- read.csv("bladder_sample_sheet.csv", row.names = 1)
table(pheno$Sample.Type)


# 3. Explore the Data

dim(data)


# 4. Explore the Data Distribution Using Histogram Plot

hist(x = data, col = "orange", main = "Histogram")


# 5. Scaling the Data Using log2 Transformation to Better Visulization
# we use (+1) to avoid the infinity character when we log zero values

hist(x = log2(data+1), col = "orange", main = "Histogram")

# 6. Making Sure that col of Gene Expression Matrix are the Same of pheno rows

pheno <- pheno[colnames(data),]


# 7. Preparing Data for DESeq2
# saving the names of genes in variable because apply function will turn them into null

genes <- row.names(data)

# converting data values to integers for DESeq2
data <- apply(data,2,as.integer)

# view data
head(data)

# bring back the row names
row.names(data) <- genes


# 8. Differential Expression Analysis Using DESeq2

# specify how many conditions according to phenotype table
cond1 <- "Solid Tissue Normal"
cond2 <- "Primary Tumor"

# create DESeq dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = pheno, design = ~ Sample.Type)

# run DESeq2 workflow
dds.run <- DESeq(dds)

# to make the res based on specified conditions
res <- results(dds.run, contrast = c("Sample.Type", cond1, cond2))


# 9. Remove Nulls

res <- as.data.frame(res[complete.cases(res), ])


# 10. Chose Statistically Significant Differentially Expressed Genes (DEGS)

desq.deg <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 2,]


# 11. Saving the Output

write.csv(as.matrix(desq.deg), file = "desq.deg.csv", quote = F, row.names = T)

