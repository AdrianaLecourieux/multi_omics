library(biomaRt)
library(ensembldb)

data <- readRDS("data_tcga.Rds")
mRNa_data <- data$mRNA
prot_data <- data$prot
miRNA_data <- data$miRNA

sample_info_data <- data$sample_info

coef.var <- function(x){
  c.var = sd(x)/mean(x)
}

# mRNA ----
data_mRNA <- mRNa_data[,colSums(mRNa_data) > 10]
coef.mRNA <- as.numeric(lapply(data_mRNA, coef.var))
data.filtered_mRNA <- data_mRNA[,abs(coef.mRNA) > 0.2]
data.scaled_mRNA <- scale(data.filtered_mRNA, center = TRUE, scale = TRUE)

# Protein ----
data_prot <- prot_data[,colSums(prot_data) > 10]
coef.prot <- as.numeric(lapply(data_prot, coef.var))
data.filtered_prot <- data_prot[,abs(coef.prot) > 0.2]
data.scaled_prot <- scale(data.filtered_prot, center = TRUE, scale = TRUE)

# most variable genes ----
max_mRNA <- quantile(coef.mRNA, 0.95)
most_variable_genes <- colnames(mRNa_data[which(coef.mRNA >= max_mRNA)])

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters="ensembl_gene_id",
                attributes=c("ensembl_gene_id", "hgnc_symbol"),
                values=most_variable_genes, mart=mart)

most_variable_genes_id <- G_list$hgnc_symbol

prots <- colnames(prot_data)

