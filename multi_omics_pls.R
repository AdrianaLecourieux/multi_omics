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
View(m_mRNA)

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


#analyse_integration_pls_comparaison

#pls1_comparaison_miRNA_mRNA

pls1.result <- pls(data.scaled_miRNA,data.scaled_mRNA , ncomp = 2 ) # run the method
plotIndiv(pls1.result)   # plot the samples
plotVar(pls1.result , var.names = FALSE)  # plot the variables

#pls2_comparaison_miRNA_prot

pls2.result <- pls(data.scaled_miRNA,data.scaled_prot , ncomp = 2 ) 
plotIndiv(pls2.result) 
plotVar(pls2.result , var.names = FALSE)

#pls3_comparaison_mRNA_prot

pls3.result <- pls(data.scaled_mRNA,data.scaled_prot , ncomp = 2 ) 
plotIndiv(pls3.result)   
plotVar(pls3.result , legend = c("mRNA", "prot") ,  var.names = FALSE)




#Initial sPLS model_basic_sPLS

#sPLS_model_miRNA_mRNA

pls1 <- spls(data.scaled_miRNA,data.scaled_mRNA , 
             keepX = list.keepX, keepY = list.keepY) 
pls1
names(pls1)
pls1$variates$X
plotIndiv(pls1)

plotVar(pls1, cutoff = 0.5, title = "(a) miRNA vs mRNA", 
        legend = c("miRNA", "mRNA"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))


#sPLS_model_miRNA_prot

pls2 <- spls(data.scaled_miRNA,data.scaled_prot , 
             keepX = list.keepX, keepY = list.keepY) 
pls2
names(pls2)
pls1$variates$X
plotIndiv(pls2)

plotVar(pls2, cutoff = 0.5, title = "(a) miRNA vs PROT", 
        legend = c("miRNA", "prot"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))


#sPLS_model_mRNA_prot

pls3 <- spls(data.scaled_mRNA,data.scaled_prot , 
             keepX = list.keepX, keepY = list.keepY) 
pls3
names(pls3)
pls1$variates$X
plotIndiv(pls3)

plotVar(pls2, cutoff = 0.5, title = "(a) mRNA vs PROT", 
        legend = c("mRNA", "prot"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

#Tuning sPLS
#Selecting the number of components

spls_prot_gene <- spls(data.scaled_mRNA,data.scaled_prot , ncomp = 10 , mode = 'regression')


per_spls_prot_gene <- perf(spls_prot_gene, validation = 'Mfold',
                        folds = 10, nrepeat = 5) 

plot(per_spls_prot_gene , criterion = 'Q2.total')

#vu le resultat du graphe on va réduire le nombre de components à 2



list.keepX = c(seq(1, 10, 1)) 
list.keepY = c(seq(1, 10, 1))

spls_final<- tune.spls(data.scaled_mRNA,data.scaled_prot, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 10, # use 10 folds
                             mode = 'regression', measure = 'cor') 


plot(spls_final) 

# determiner la valeur brute du keepX et du keepY

optimal.keepX <- spls_final$choice.keepX
print(optimal.keepX)


optimal.keepY <- spls_final$choice.keepY
print(optimal.keepY)


# MODEL FINAL, la sPLS
#etant_donné qu'au niveau du tp on doit analyser sPLS avec les gènes et les protéines (10 features/bloc composante 1, 5 pour composante 2)


final<- spls(data.scaled_mRNA,data.scaled_prot, ncomp = 2, 
                         keepX = c(10, 5),
                         keepY = c(10, 5),
                         mode = "regression")



#final_plots

#corrélation_sample_plot


plotVar(final, cex = c(3,4), var.names = c(FALSE, TRUE) )


#Cluster Image Map (CIM)

cim(final, comp = 1:2, xlab = "gene", ylab = "prot" , )











