library(biomaRt)
library(ensembldb)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(rentrez)

data <- readRDS("data_tcga.Rds")
data_mRNA <- data$mRNA
data_prot <- data$prot
data_miRNA <- as.data.frame(data$miRNA)

sample_info_data <- data$sample_info


gene_symbol <- mapIds(org.Hs.eg.db, keys=colnames(data_mRNA),
                      column="SYMBOL", keytype="ENSEMBL",
                      multiVals="first")

colnames(data_mRNA) <- gene_symbol

prot_id <- mapIds(org.Hs.eg.db, keys=colnames(data_mRNA),
                  column="UNIPROT", keytype="SYMBOL",
                  multiVals="first")

symbol <- mapIds(org.Hs.eg.db, keys=colnames(data_prot),
                 column="SYMBOL", keytype="UNIPROT",
                 multiVals="first")

na_id <- complete.cases(symbol)

gene_from_prot <- symbol[na_id]

intersect(colnames(data_prot), prot_id)
intersect(colnames(data_mRNA), gene_from_prot)

not_uniprot <- names(symbol[!na_id])


r_search <- sapply(not_uniprot, function(x) entrez_search(db="protein", term=paste0(x, "[Protein Name]"))$ids[1])

r_search <- r_search[lengths(r_search) != 0]

entrez_id <- sapply(r_search, function(x) as.character(entrez_summary(db="protein", id=x)$taxid))

entrez_id <- na.omit(entrez_id)

gene_from_non_uniprot <- mapIds(org.Hs.eg.db, keys=entrez_id, column="SYMBOL", keytype="ENTREZID", multiVals="first")

gene_from_non_uniprot <- na.omit(gene_from_non_uniprot)

intersect(colnames(data_mRNA), gene_from_non_uniprot)

