install.packages("BiocManager")
install.packages("devtools")
BiocManager::install("netOmics")
BiocManager::install("GO.db")
library(devtools)
library(mixOmics)
library(ggplot2)
library(netOmics)
library(GO.db)


#donnees_150_échatillons_patients

data_tcga <- readRDS("C:/Users/m_graa/Desktop/projet_omique/data_tcga.Rds")
View(data_tcga)

prot_data <- data_tcga[["prot"]]
m_mRNA <- data_tcga[["mRNA"]]
miRNA <- data_tcga[["miRNA"]]

View(data_tcga)
View(miRNA)
View(sample_data)
View(prot_data)


#les_samples

sample_data <- data_tcga$sample_info
dim(sample_data)

#presentation_des_distributions
summary(as.factor(sample_data$Y))
boxplot(m_mRNA, col = as.factor(sample_data$Y), main = "boxplot mRNA")
boxplot(miRNA, col = as.factor(sample_data$Y), main = "boxplot miRNA")
boxplot(prot_data, col = as.factor(sample_data$Y), main = "boxplot protdata")




#filtration_des_genes_peu_variants

#mRNA

data_rna <- m_mRNA[,colSums(m_mRNA) > 400]
coef.var <- function(x){
  c.var = sd(x)/mean(x)
}
coef.mRNA <- as.numeric(lapply(m_mRNA, coef.var))
hist(coef.mRNA)


rna_filtered <- m_mRNA[,abs(coef.mRNA) > 0.15]

dim(data_rna)
dim(rna.filtered)

#PROT

data_prot <- prot_data[,colSums(prot_data) > 95]
coef.var <- function(x){
  c.var = sd(x)/mean(x)
}
coef.prot <- as.numeric(lapply(prot_data, coef.var))
hist(coef.prot)


prot_filtered <- prot_data[,abs(coef.prot) > -50]


dim(data_prot)
dim(prot_filtered)


#miRNA_sans_filtre


colSums(miRNA)

hist(colSums(miRNA), breaks = 20)


coef_var <- function(x){
  c_var = sd(x)/mean(x)
}
coef_miRNA <- as.numeric(lapply(as.data.frame(miRNA), coef_var))
hist(coef_miRNA)

miRNA_filtered <- miRNA[,coef_miRNA > 0.05]



dim(miRNA_filtered)
dim(miRNA)



#scale

data.scaled_mRNA <- scale(rna.filtered, center = TRUE, scale = TRUE)

data.scaled_prot <- scale(prot_filtered, center = TRUE, scale = TRUE)

data.scaled_miRNA <- scale(miRNA_filtered, center = TRUE, scale = TRUE)  

boxplot(data.scaled_mRNA, col = as.factor(sample_data$Y), main = "boxplot mRNA after scale")
boxplot(data.scaled_prot, col = as.factor(sample_data$Y), main = "boxplot prot after scale")
boxplot(data.scaled_miRNA, col = as.factor(sample_data$Y), main = "boxplot prot after scale")


#pls

#pls1_miRNA_mRNA

pls.result <- pls(data.scaled_miRNA,data.scaled_mRNA ) # run the method
plotIndiv(pls.result)   # plot the samples
plotVar(pls.result)  # plot the variables

#pls2_miRNA_prot

pls.result <- pls(data.scaled_miRNA,data.scaled_prot ) 
plotIndiv(pls.result) 
plotVar(pls.result)

#pls3_mRNA_PRot

pls.result <- pls(data.scaled_mRNA,data.scaled_prot ) 
plotIndiv(pls.result)   
plotVar(pls.result)

#sPLS# select arbitrary values of features to keep

list.keepX = c(27, 27) 
list.keepY = c(27, 27)

# generate three pairwise PLS models

pls1 <- spls(data.scaled_miRNA,data.scaled_mRNA , 
             keepX = list.keepX, keepY = list.keepY) 
pls1
names(pls1)
pls1$variates$X
plotIndiv(pls1)
pls1.perf.res <- perf(pls1, validation = "loo")
pls1.perf.res





pls2 <- spls(data.scaled_miRNA,data.scaled_prot , 
             keepX = list.keepX, keepY = list.keepY)

pls2
names(pls2)
pls2$variates$X
plotIndiv(pls2)
pls2.perf.res <- perf(pls2, validation = "loo")
pls2.perf.res
#plot(pls2.perf.res$Q2.total)


pls3 <- spls(data.scaled_mRNA,data.scaled_prot , 
             keepX = list.keepX, keepY = list.keepY)

pls3
names(pls3)
pls3$variates$X
plotIndiv(pls3)
pls3.perf.res <- perf(pls3, validation = "loo")
pls3.perf.res



# plot features of first PLS
plotVar(pls1, cutoff = 0.5, title = "(a) miRNA vs mRNA", 
        legend = c("miRNA", "mRNA"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of second PLS
plotVar(pls2, cutoff = 0.5, title = "(b) miRNA vs proteomics", 
        legend = c("miRNA", "proteomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of third PLS
plotVar(pls3, cutoff = 0.5, title = "(c) mRNA vs proteomics", 
        legend = c("mRNA", "proteomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))




#TP#Faire une première analyse sPLS avec les gènes et les protéines (10 features/bloc composante 1, 5 pour composante 2)


prot_gene <- spls(X = data.scaled_prot , Y = data.scaled_mRNA, ncomp = 2, keepX = c(1,5), keepY = c(1, 5))

prot_gene
names(prot_gene)

prot_gene$variates$X

plotIndiv(prot_gene)

plotVar(prot_gene)

prot_gene_res <- perf(prot_gene, validation = "loo")
prot_gene_res

















