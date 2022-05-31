if(!require("BiocManager",quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(version="3.15")
BiocManager::install("GEOquery")
BiocManager::install('WGCNA', force=TRUE)
BiocManager::install('graph',force=TRUE)
install.packages("psych")
install.packages("caret")
install.packages('randomForest')
install.packages("ggm")
library(GEOquery)
library(psych)
library(caret)
library(randomForest)
library(ggm)
library(WGCNA)
library(GGMselect)
library(glasso)
options(strinsAsFactors=FALSE)
allowWGCNAThreads()


setwd("Users/suminlee/Desktop/Probability_Graph_Model/term_project/anal")

# download .gz file and load the data
gse49710 <- getGEO("GSE49710",GSEMatrix=TRUE)[[1]]

# gse49710 <- getGEO("./GSE49710_series_matrix.txt.gz", GSEMatrix=TRUE)[[1]]
d_phenomenon <- pData(phenoData(gse49710))
d_gene <- pData(featureData(gse49710))
d_express <- as.data.frame(exprs(gse49710))

#========================
# print basic information
#========================
sprintf("Number of patients: %d",
        dim(d_phenomenon)[1])
sprintf("Number of genes: %d",
        length(d_gene$GeneSymbol))
assigned <- (d_gene$GeneSymbol != "")
sprintf("Number of assigned genes: %d/%d (%.1f%%)",
        sum(assigned),
        length(d_gene$GeneSymbol),
        sum(assigned)/length(d_gene$GeneSymbol)*100.0)
sprintf("Number of genes with missing data: %d",
        dim(d_express)[1]-dim(na.omit(d_express))[1])
sprintf("Number of assigned genes without missing data: %d (omitted: %d)",
        dim(na.omit(d_express[assigned,]))[1],
        sum(assigned)-dim(na.omit(d_express[assigned,]))[1])
d_express <- na.omit(d_express[assigned,])
d_gene <- d_gene[assigned,]
splrintf("Unique genes: %d",
         lenght(unique(d_gene$GeneSymbol)))
boxplot(d_express[,300:315])

#========================
# high risk prediction
# with mycn amplification
#========================
# MYCN-high_risk, MYCN-progression, high_risk-death_from_disease correlation
mycn <- d_phenomenon$'mycn status:ch1'
mycn[mycn == "N/A"] <- "0"
mycn <- as.integer(mycn)
risk <- as.integer(d_phenomenon$'high risk:ch1')
progression <- as.integer(d_phenomenon$'progression:ch1')
death <- as.integer(d_phenomenon$'death from disease:ch1')

sprintf("Cohen Kappa (mycn, risk) = %.4f",
        cohen.kappa(cbind(mycn, risk))$kappa)
sprintf("Cohen Kappa (mycn, progression) = %.4f",
        cohen.kappa(cbind(mycn, progression))$kappa)
sprintf("Cohen Kappa (risk, death) = %.4f",
        cohen.kappa(cbind(risk, death))$kappa)

# mycn amplification as high risk predictor
confusionMatrix(as.factor(mycn),as.factor(risk))

# duplication example "CUL4A"
cul4a <- d_express[d_gene$GeneSymbol == "CUL4A",]
cor(t(cul4a))

#========================
# high risk prediction 
# with Random Forest
#========================
# leave the first data if gene is duplicated
d_express_unique <- d_express[!duplicated(d_gene$GeneSymbol),]
# split train-test set
train_size <- 0.8*dim(d_express_unique)[2]
train_set <- d_express_unique[,1:train_size]
test_set <- d_express_unique[,train_size:dim(d_express_unique)[2]]

# random forest (for gene selection)
rf <- randomForest(x=t(train_set), y=as.factor(risk[1:train_size]))
train_pred <- predict(rf, t(train_set))
confusionMatrix(train_pred, as.factor(risk[1:train_size]))
test_pred <- predict(rf, t(test_set))
confusionMatrix(test_pred, as.factor(risk[train_size:dim(d_express_unique)[2]]))

impt <- importance(rf)
top_n <- 1000
sum(head(sort(impt/sum(impt),decreasing=TRUE), n = top_n))
sorted_gene <- sort(impt/sum(impt),decreasing=TRUE,index.return=TRUE)
top_n_id <- head(sorted_gene$ix, n= top_n)
top_n_importance <- head(sorted_gene$x, n= top_n)

#========================
# Graph construction
#========================
d_express_top <- d_express_unique[c(top_n_id),]

# covariance matrix
cov_mat = cov(d_express_unique)

# ggm package
# need adjacency matrix
# glasso_ggm <- fitCovGraph(cov_mat)

# glasso package
inv_cov <- glasso(d_express, rho=0.1)

# GGMselect package
data_mat <- data.matrix(d_express_unique)
ggm <- selectFast(data_mat, family = "LA", verbose = TRUE)




#========================
# WCGNA
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#========================
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(d_express_top, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
net = blockwiseModules(d_express_top, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "NeuralBalstomaTOM", 
                       verbose = 3)

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "network.RData")


#========================
# high risk prediction
# with PGM
#========================
# TODO


