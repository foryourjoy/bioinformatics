install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db",force = TRUE)
BiocManager::install("apeglm")
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)

################################################################################
### 1. Create DESeq results ####################################################
################################################################################

countData <- read.csv("gene_count_matrix.csv", row.names = 1)
countData %>% head(10)
countData %>% is.na() %>% sum()
colData <- data.frame(condition = factor(c("Normal","Normal","Normal","Cancer", "Cancer", "Cancer")),  row.names = colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition )
dds %>% head()
dds <- DESeq(dds)

resultsNames(dds)

deseq2_result <- results(dds, name="condition_Normal_vs_Cancer")
summary(deseq2_result)
nrow(deseq2_result)

################################################################################
### 2. variation shrinkage plot ################################################
################################################################################

resLFC <- lfcShrink(dds, coef="condition_Normal_vs_Cancer", type="apeglm")
resLFC
plotMA(resLFC, ylim=c(-2,2))


################################################################################
### 3. t-test ##################################################################
################################################################################

group <- factor(c("Normal", "Normal","Normal", "Cancer", "Cancer","Cancer"))

p_values <- apply(countData, 1, function(x) { 
  normal <- as.numeric(x[1:3]) # 정상인 샘플 발현량
  cancer <- as.numeric(x[4:6]) # 암환자 샘플 발현량
  
  # Check if the values in either group are constant
  if (var(normal) == 0 || var(cancer) == 0) {
    return(NA) # Return NA if values are constant
  } else {
    return(t.test(normal, cancer)$p.value) # Perform t-test and return p-value
  }
})

t_test_results <- data.frame( Gene = rownames(countData), P_Value = p_values ) %>% arrange(P_Value)
write.csv(t_test_results, "t_test_results.csv", row.names = FALSE)
summary(t_test_results)

significant_genes <- subset(t_test_results, P_Value < 0.05)
write.csv(significant_genes, "t_test_results_deg.csv", row.names = FALSE)
nrow(significant_genes)
significant_genes %>% head(10)

## plot
library(ggplot2)
top_genes <- significant_genes %>%
  arrange(P_Value) %>%
  head(6)

top_gene_data <- countData[rownames(countData) %in% top_genes$Gene, ]
top_gene_data_long <- data.frame(Gene = rownames(top_gene_data), top_gene_data) %>%
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(Group = ifelse(as.numeric(gsub("[^0-9]", "", Sample)) <= 3, "Normal", "Cancer"))

top_gene_data_long <- merge(top_gene_data_long, top_genes, by.x = "Gene", by.y = "Gene")

# Plot (2x3) matrix of boxplots with p-values
ggplot(top_gene_data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Gene, nrow = 2, ncol = 3, scales = "free") +  # 2x3 matrix of plots
  geom_text(aes(x = 1.5, y = max(Expression), label = paste0("p = ", round(P_Value, 4))),
            vjust = -0.5, size = 3) +  # Annotate with p-value
  labs(title = "Differential Gene Expression for Top 6 Significant Genes",
       x = "Group", y = "Expression Level") +
  theme_minimal() +
  scale_fill_manual(values = c("Normal" = "lightblue", "Cancer" = "pink"))


################################################################################
### 4. DESeq2 ##################################################################
################################################################################

write.csv(deseq2_result, "deseq2_results.csv", row.names = FALSE)
deseq2_genes <- deseq2_result[!is.na(deseq2_result$padj) & deseq2_result$padj < 0.05, ]
deseq2_genes <- data.frame(Gene = rownames(deseq2_genes), deseq2_genes, row.names = NULL)
nrow(deseq2_genes)
deseq2_genes %>% head(10)

write.csv(deseq2_genes, "deseq2_results_deg.csv", row.names = FALSE)

## plot
# Sort deseq2_genes by padj and select top 6 genes
top_genes <- deseq2_genes %>%
  arrange(padj)

top_genes %>% head(10)

# Extract the corresponding expression data for the top genes
top_gene_data <- countData[rownames(countData) %in% top_genes$Gene, ]

# Reshape data for ggplot
top_gene_data_long <- data.frame(Gene = rownames(top_gene_data), top_gene_data) %>%
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(Group = ifelse(as.numeric(gsub("[^0-9]", "", Sample)) <= 3, "Normal", "Cancer"))

# Merge p-values for annotation
top_gene_data_long <- merge(top_gene_data_long, top_genes, by.x = "Gene", by.y = "Gene")

# Plot (2x3) matrix of boxplots with p-values
ggplot(top_gene_data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Gene, nrow = 2, ncol = 3, scales = "free") +  # 2x3 matrix of plots
  geom_text(aes(x = 1.5, y = max(Expression), label = paste0("p = ", round(padj, 4))),
            vjust = -0.5, size = 3) +  # Annotate with p-value
  labs(title = "Differential Gene Expression for Top 6 Significant Genes",
       x = "Group", y = "Expression Level") +
  theme_minimal() +
  scale_fill_manual(values = c("Normal" = "lightblue", "Cancer" = "pink"))


################################################################################
### 5. DESeq2 VS t-test ########################################################
################################################################################

deseq2_gene_names <- deseq2_genes%>% filter(padj < 0.05) %>% select(Gene)
ttest_genes <- t_test_results %>% filter(P_Value < 0.05) %>% select(Gene)
common_genes <- intersect(deseq2_gene_names, ttest_genes)

deseq2_only <- setdiff(deseq2_gene_names, ttest_genes)
ttest_only <- setdiff(ttest_genes, deseq2_gene_names)


################################################################################
### 6. GO analysis #############################################################
################################################################################

### for deseq2 only

library(stringr)
gene_names <- str_extract(deseq2_only$Gene, "(?<=\\|).+")
gene_names <- unique(gene_names)

go_deseq2 <- enrichGO( 
  gene = gene_names, 
  OrgDb = org.Hs.eg.db,  
  keyType = "SYMBOL",  
  ont = "BP",  
  pvalueCutoff = 0.05)

if (!is.null(go_deseq2)) {
  # Convert to data frame
  go_results <- as.data.frame(go_deseq2)
  
  # Save to CSV
  write.csv(go_results, file = "GO_enrichment_results_deseq2.csv", row.names = FALSE)
  print("GO enrichment results saved to 'GO_enrichment_results.csv'")
} else {
  message("No significant GO terms found.")
}


# Load necessary libraries for plotting
library(ggplot2)

# Check if the go_deseq2 result is not NULL
if (!is.null(go_deseq2)) {
  # Convert to data frame
  go_results <- as.data.frame(go_deseq2)
  
  # Select top 10 enriched terms based on p-value
  top_terms <- head(go_results[order(go_results$pvalue), ], 10)
  
  # Create a bar plot
  ggplot(top_terms, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +  # Flip coordinates for better readability
    labs(title = "Top 10 Enriched GO Terms",
         x = "GO Terms",
         y = "-log10(p-value)") +
    theme_minimal()
}

####################
### for t-test only

gene_names_t <- str_extract(ttest_only$Gene, "(?<=\\|).+")
gene_names_t <- unique(gene_names_t)
length(gene_names_t)

go_ttest <- enrichGO( 
  gene = gene_names_t, 
  OrgDb = org.Hs.eg.db,  
  keyType = "SYMBOL",  
  ont = "BP",  
  pvalueCutoff = 0.05)

if (!is.null(go_ttest)) {
  # Convert to data frame
  go_results_t <- as.data.frame(go_ttest)
  
  # Save to CSV
  write.csv(go_results_t, file = "GO_enrichment_results_ttest.csv", row.names = FALSE)
} else {
  message("No significant GO terms found.")
}


# Load necessary libraries for plotting
library(ggplot2)

# Check if the go_deseq2 result is not NULL
if (!is.null(go_ttest)) {
  # Convert to data frame
  go_results_t <- as.data.frame(go_ttest)
  
  # Select top 10 enriched terms based on p-value
  top_terms_t <- head(go_results_t[order(go_results_t$pvalue), ], 10)
  
  # Create a bar plot
  ggplot(top_terms_t, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +  # Flip coordinates for better readability
    labs(title = "Top 10 Enriched GO Terms",
         x = "GO Terms",
         y = "-log10(p-value)") +
    theme_minimal()
}


################################################################################
### 7. FPKM & TPM from abundance.txt file ######################################
################################################################################
library(dplyr)

abd_s1 = read.table('sample1_abundance.txt',sep='\t', header = TRUE)
abd_s1 %>% head()
abd_s1 %>% arrange(desc(TPM)) %>% head(10)

abd_s4 = read.table('sample4_abundance.txt',sep='\t', header = TRUE)
abd_s4 %>% head()
abd_s4 %>% arrange(desc(TPM)) %>% head(10)
