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

not_uniprot_df <- data.frame(prot = not_uniprot, code = 0)
for (i in seq(length(not_uniprot))){
  r_search <- entrez_search(db="protein", term=paste0(not_uniprot[i], " AND Homo sapiens[Organism] "))$ids
  r_search <- r_search[lengths(r_search) != 0]
  
  multi_summs <- entrez_summary(db="protein", id=r_search[1:10])
  
  entrez_id <- as.character(extract_from_esummary(multi_summs, "extra"))
}
r_search <- sapply(not_uniprot, function(x) entrez_search(db="protein", term=paste0(x, " AND Homo sapiens[Organism] "))$ids)
r_search <- r_search[lengths(r_search) != 0]

multi_summs <- lapply(r_search, function(x) entrez_summary(db="protein", id=unlist(x, use.names=FALSE))[1:10])
multi_summs <- multi_summs[lengths(multi_summs) != 0]

entrez_id <- sapply(multi_summs, function(x) as.character(extract_from_esummary(x, "extra")))

gene_from_non_uniprot <- mapIds(org.Hs.eg.db, keys=entrez_id, column="SYMBOL", keytype="ENTREZID", multiVals="first")

gene_from_non_uniprot <- na.omit(gene_from_non_uniprot)

intersect(colnames(data_mRNA), gene_from_non_uniprot)

uniprot_of_non_uniprot <- mapIds(org.Hs.eg.db, keys=entrez_id, column="UNIPROT", keytype="ENTREZID", multiVals="first")

uniprot_of_non_uniprot <- na.omit(uniprot_of_non_uniprot)

uniprot_df <- data.frame(entrezid = names(uniprot_of_non_uniprot), uniprot = uniprot_of_non_uniprot)
entrez_df <- data.frame(entrezid = as.vector(entrez_id), name = names(entrez_id))

replace_df <- merge(uniprot_df, entrez_df, by="entrezid")
