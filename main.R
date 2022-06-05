pack_s<-list("BiocManager", "psych", "caret", "randomForest", "ggm",
               "GGMselect", "glasso","plot.matrix")
bpack_s<-list("GEOquery","WGCNA","graph")
for (p in pack_s){
  if(!p%in%rownames(installed.packages())){
    install.packages(p)
  }
}
for (b in bpack_s){
  if(!b%in%rownames(installed.packages())){
    BiocManager::install(b)
  }
}
library(GEOquery)
library(psych)
library(caret)
library(randomForest)
library(ggm)
library(WGCNA)
library(GGMselect)
library(glasso)
library(plot.matrix)
library(RColorBrewer)
options(strinsAsFactors=FALSE)
allowWGCNAThreads()


# download .gz file and load the data
gse49710 <- getGEO("GSE49710",GSEMatrix=TRUE)[[1]]

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
sprintf("Unique genes: %d",
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
sprintf("gene duplication example, CUL4A")
cor(t(cul4a))

#========================
# high risk prediction 
# with Random Forest
#========================
# leave the first data if gene is duplicated
d_express_unique <- d_express[!duplicated(d_gene$GeneSymbol),]
# split train-test set
# balance high-risk and low-risk ratio of train-test sets
high_risk_express <- d_express_unique[risk==1]
low_risk_express <- d_express_unique[risk==0]
high_train_size <- as.integer(0.8*sum(risk==1))
low_train_size <- as.integer(0.8*sum(risk==0))
train_set <- cbind(high_risk_express[,1:high_train_size],
                     low_risk_express[,1:low_train_size])
test_set <- cbind(high_risk_express[,(high_train_size+1):sum(risk==1)], 
                    low_risk_express[,(low_train_size+1):sum(risk==0)])
train_risk <- c(rep(c(1),high_train_size), rep(c(0),low_train_size))
test_risk <- c(rep(c(1),sum(risk==1)-high_train_size), 
                 rep(c(0),sum(risk==0)-low_train_size))

# random forest (for gene selection)
if (file.exists("rf.RData")){
  load("rf.RData")
} else {
  rf <- randomForest(x=t(train_set), y=as.factor(train_risk));
  save(rf, file="rf.RData")
}
train_pred <- predict(rf, t(train_set))
confusionMatrix(train_pred, as.factor(train_risk))
test_pred <- predict(rf, t(test_set))
confusionMatrix(test_pred, as.factor(test_risk))

impt <- importance(rf)
top_n <- 1000
sum(head(sort(impt/sum(impt),decreasing=TRUE), n = top_n))
sorted_gene <- sort(impt/sum(impt),decreasing=TRUE,index.return=TRUE)
top_n_id <- head(sorted_gene$ix, n= top_n)
top_n_importance <- head(sorted_gene$x, n= top_n)

#========================
# Graph construction
#========================
# top 1000 genes
#========================
d_express_topn <- d_express_unique[c(top_n_id),]

# standardization
topn_scale <- t(scale(d_express_topn))

# covariance matrix
topn_cov_mat <- cov(topn_scale)

# glasso package
if (file.exists("topn_inv_cov.RData")){
  load("topn_inv_cov.RData")
};
else {
  topn_inv_cov <- glasso(topn_scale, rho=0.1, nobs=dim(topn_scale)[1])
  save(topn_inv_cov,file="topn_inv_cov.RData")
}
# visualize inverse covariance matrix's zero pattern up to top100 genes
# more than 400 dimension too mash
jpeg(file='topn_inv_cov.jpeg')
plot(topn_inv_cov[1:200,1:200], border=NA) 
dev.off()

# GGMselect package
m_topn_cov <- data.matrix(topn_cov_mat)
topn_ggm <- selectFast(m_topn_cov, family = "LA", verbose = TRUE)

#========================
# entire genes
#========================
#unique_scale <- scale(t(d_express))
#unique_cov_mat <- cov(uqnie_scale)
#unique_inv_cov <- glasso(unique_scale, rho=0.1, nobs=dim(unique_scale)[1])
#m_unique_cov <- data.matrix(unique_cov_mat)
#unique_ggm <- selectFast(m_unqiue_cov, family = "LA", verbose = TRUE)

#========================
# high risk prediction
# with PGM
#========================
# posterior P(G|x)
# GGM assumes multivariate Gaussian
#========================

# low-risk group
topn_lowrisk <- low_risk_express[c(top_n_id),1: low_train_size]
topn_lowrisk_scale <- scale(t(topn_lowrisk))
topn_lowrisk_cov <- cov(topn_lowrisk_scale)
if (file.exists("topn_lowrisk_inv_cov.RData")){
  load("topn_lowrisk_inv_cov.RData")
} else {
  topn_lowrisk_inv_cov <- glasso(topn_lowrisk, rho=0.1, nobs=dim(topn_lowrisk_scale)[1]);
  save(topn_lowrisk_inv_cov, "topn_lowrisk_inv_cov.RData");
}


# high-risk group
topn_highrisk <- high_risk_express[c(top_n_id),1:high_train_size]
topn_highrisk_scale <- scale(t(topn_highrisk))
topn_highrisk_cov <- cov(topn_highrisk_scale)
if (file.exists("topn_highrisk_inv_cov.RData")){
  load("topn_highrisk_inv_cov.RData")
} else {
  topn_highrisk_inv_cov <- glasso(topn_highrisk, rho=0.1, nobs=dim(topn_highrisk_scale)[1]);
  save(topn_highrisk_inv_cov, "topn_highrisk_inv_cov.RData");
}


# Posterior
multiGaussian <- function(x, cov, inv_cov){
  dimension <- length(x)
  det <- determinant(cov)
  det <- sqrt(exp(det$modulus[1]))
  result <- exp(-0.5*(x%*%inv_cov%*%x))/(sqrt(2*pi)**dimension)/det
  return(result)
}
low_transform <- function(x){
  low_mu <- attr(topn_lowrisk_scale,"scaled:center")
  low_sd <- attr(topn_lowrisk_scale,"scaled:scale")
  result <- (x-low_mu)/low_sd
  return(result)
}
high_transform <- function(x){
  low_mu <- attr(topn_highrisk_scale,"scaled:center")
  low_sd <- attr(topn_highrisk_scale,"scaled:scale")
  result <- (x-high_mu)/high_sd
  return(result)
}

GGMpred <- function(x){
  # probability of high risk given x
  # P(high), P(low)
  p_high <- high_train_size/(high_train_size+low_train_size)
  high_x <- high_transform(x)
  post_high <- multiGaussian(high_x, topn_highrisk_cov, topn_highrisk_inv_cov)*p_high
  
  p_low <- low_train_size/(high_train_size+low_train_size)
  low_x <- low_transform(x)
  post_low <- multiGaussian(low_x, topn_lowrisk_cov, topn_lowrisk_inv_cov)*p_low
  
  return(p_high/(p_high+p_low))
}


#========================
# factor analysis
#========================
fa <- factanal(train_set[1:1000,], factors=2)


#========================
# expected calibration error
#========================







quit(save="yes")



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



