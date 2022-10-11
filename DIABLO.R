library(mixOmics)
library(biomaRt)
library(igraph)

# settings for parallelization
BPPARAM <- BiocParallel::SnowParam(workers = max(parallel::detectCores()-1, 2))

data <- readRDS("data_tcga.Rds")
data_mRNA <- as.data.frame(data$mRNA)
data_prot <- as.data.frame(data$prot)
data_miRNA <- as.data.frame(data$miRNA)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters="ensembl_gene_id",
                attributes=c("ensembl_gene_id", "hgnc_symbol"),
                values=colnames(data_mRNA), mart=mart)

colnames(data_mRNA) <- G_list$hgnc_symbol


sample_info <- data$sample_info

coef.var <- function(x){
  c.var = sd(x)/mean(x)
}

# mRNA ----
# remove genes without enought expression
data_mRNA <- data_mRNA[, colSums(data_mRNA) > 10]

coef_mRNA <- as.numeric(lapply(data_mRNA, coef.var))
# remove genes without enough variability
data_filtered_mRNA <- data_mRNA[, abs(coef_mRNA) > 0.2]
# scale
data_scaled_mRNA <- scale(data_filtered_mRNA, center = TRUE, scale = TRUE)

# Protein ----
# remove genes without enought expression
data_prot <- data_prot[, colSums(data_prot) > 10]

coef_prot <- as.numeric(lapply(data_prot, coef.var))
# remove genes without enough variability
data_filtered_prot <- data_prot[, abs(coef_prot) > 0.2]
# scale
data_scaled_prot <- scale(data_filtered_prot, center = TRUE, scale = TRUE)

# miRNA ----
# remove genes without enought expression
data_miRNA <- data_miRNA[, colSums(data_miRNA) > 10]

coef_miRNA <- as.numeric(lapply(data_miRNA, coef.var))
# remove genes without enough variability
data_filtered_miRNA <- data_miRNA[, abs(coef_miRNA) > 0.2]
# scale
data_scaled_miRNA <- scale(data_filtered_miRNA, center = TRUE, scale = TRUE)

# data ----
X <- list(mRNA = data_scaled_mRNA,
          protein = data_scaled_prot,
          miRNA = data_scaled_miRNA)

Y <- as.factor(sample_info$Y)

# PLSDA and SPLSDA ----
result.diablo.tcga <- block.plsda(X, Y)
plotIndiv(result.diablo.tcga)
plotVar(result.diablo.tcga)

# TODO: update the values with sPCA ?
list.keepX <- list(miRNA = c(16, 17), mRNA = c(18, 5), protein = c(5, 5)) 

result.sparse.diablo.tcga <-  block.splsda(X, Y, keepX = list.keepX) 

plotLoadings(result.sparse.diablo.tcga, ncomp = 1) 
plotIndiv(result.sparse.diablo.tcga)
plotVar(result.sparse.diablo.tcga)

# DIABLO ----
# square matrix filled with 0.1s
design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
# set diagonal to 0s
diag(design) <- 0

# basic diablo model
basic.diablo.model <- block.splsda(X = X, Y = Y, ncomp = 5, design = design) 

perf.diablo <- perf(basic.diablo.model, validation = 'Mfold', 
                    folds = 10, nrepeat = 10) 

plot(perf.diablo)

# optimal ncomp value
ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 

# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote

test.keepX <- list(mRNA = c(5:9, seq(10, 18, 2), seq(20, 30, 5)), 
                   miRNA = c(5:9, seq(10, 18, 2), seq(20, 30, 5)),
                   protein = c(5:9, seq(10, 18, 2), seq(20, 30, 5)))

# run the feature selection tuning
tune.TCGA <- tune.block.splsda(X = X, Y = Y, ncomp = ncomp,
                               test.keepX = test.keepX, design = design,
                               validation = 'Mfold', folds = 10, nrepeat = 1,
                               dist = "centroids.dist", BPPARAM = BPPARAM)

# optimal values of features to retain
list.keepX <- tune.TCGA$choice.keepX

# final diablo model
final.diablo.model <- block.splsda(X = X, Y = Y, ncomp = ncomp, 
                                   keepX = list.keepX, design = design)

# final design matrix
final.design <- final.diablo.model$design

# the features selected to form the first component
selected.vars <- selectVar(final.diablo.model, block = 'mRNA', comp = 1)$mRNA$name 

# plot ----
save_and_plot <- function(dest_file, plot_function, ...){
  svg(dest_file)
  plot_function(...)
  dev.off()
}

# plotDiablo(final.diablo.model, ncomp = 1)
save_and_plot("figures_diablo/final_diablo.svg", plotDiablo, final.diablo.model, ncomp = 1)

# plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
#           title = 'DIABLO Sample Plots')
save_and_plot("figures_diablo/diablo_individuals.svg", plotIndiv,
              final.diablo.model, ind.names = FALSE,
              legend = TRUE, title = 'DIABLO Sample Plots')

# plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
save_and_plot("figures_diablo/diablo_arrow.svg", plotArrow, final.diablo.model,
              ind.names = FALSE, legend = TRUE, title = 'DIABLO')

# plotVar(final.diablo.model, var.names = FALSE, 
#         style = 'graphics', legend = TRUE,
#         pch = c(16, 17, 15), cex = c(2,2,2), 
#         col = c('darkorchid', 'brown1', 'lightgreen'))
save_and_plot("figures_diablo/diablo_variables.svg", plotVar, final.diablo.model,
              var.names = FALSE, style = 'graphics', legend = TRUE,
              col = c('darkorchid', 'brown1', 'lightgreen'),
              pch = c(16, 17, 15), cex = c(2,2,2))

# circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
#            color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
#            color.cor = c("chocolate3","grey20"), size.labels = 1.5)
save_and_plot("figures_diablo/diablo_circos.svg", circosPlot,
              final.diablo.model,cutoff = 0.7, line = TRUE,
              color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
              color.cor = c("chocolate3","grey20"), size.labels = 1.5) 

# network(final.diablo.model, blocks = c(1,2,3),
#         color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
save_and_plot("figures_diablo/diablo_relevance_network.svg", network,
              final.diablo.model, blocks = c(1,2,3), cutoff = 0.65,
              color.node = c('darkorchid', 'brown1', 'lightgreen'))

# save graph to gml
# my.network = network(final.diablo.model, blocks = c(1,2,3),
#                      color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
# write.graph(my.network$gR, file = "myNetwork.gml", format = "gml")

# plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median')
save_and_plot("figures_diablo/diablo_loading_vectors.svg", plotLoadings,
              final.diablo.model, comp = 2, contrib = 'max', method = 'median')

# cimDiablo(final.diablo.model)
save_and_plot("figures_diablo/diablo_cim.svg", cimDiablo, final.diablo.model)

# perf ----
# run repeated CV performance evaluation
perf.diablo = perf(final.diablo.model, validation = 'Mfold', 
                   M = 10, nrepeat = 10, 
                   dist = 'centroids.dist') 

perf.diablo$MajorityVote.error.rate

# auc.splsda <- auroc(final.diablo.model, roc.block = "miRNA", roc.comp = 2, print = FALSE)
save_and_plot("figures_diablo/diablo_perf.svg", auroc, final.diablo.model,
              roc.block = "miRNA", roc.comp = 2, print = FALSE)
