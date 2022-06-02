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


setwd("~/Source/work/PGM_project")

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

train_set   <- d_express_unique[,1:train_size]
train_label <- as.factor(risk[1:train_size])
test_set   <- d_express_unique[,train_size:dim(d_express_unique)[2]]
test_label <- as.factor(risk[ train_size:dim(d_express_unique)[2] ])

# random forest (for gene selection)
rf <- randomForest(x=t(train_set), y=as.factor(risk[1:train_size]))
train_pred <- predict(rf, t(train_set))
confusionMatrix(train_pred, train_label)

test_pred <- predict(rf, t(test_set))
confusionMatrix(test_pred, test_label)

impt <- importance(rf)
top_n <- 1000
sum(head(sort(impt/sum(impt),decreasing=TRUE), n = top_n))
sorted_gene <- sort(impt/sum(impt),decreasing=TRUE,index.return=TRUE)
top_n_id <- head(sorted_gene$ix, n= top_n)
top_n_importance <- head(sorted_gene$x, n= top_n)

# New tensors sorted by importance
d_express_top <- d_express_unique[c(top_n_id),]

## Feature classifier using important features

## Create random index for fair comparison 
seed_key <- strtoi("00000001")  # Seed for random sampling
set.seed(seed_key)
random_sample <- sample(x=1:length(risk), size=length(risk))

# Random indices;
rnd_tr_idx <- c(random_sample[         1:train_size  ])
rnd_ts_idx <- c(random_sample[train_size:length(risk)])

trunc_train_set   <- d_express_top[, rnd_tr_idx]
trunc_train_label <- as.factor(risk[rnd_tr_idx])
trunc_test_set   <- d_express_top[, rnd_ts_idx]
trunc_test_label <- as.factor(risk[rnd_ts_idx])

small_rf <- randomForest(x=t(trunc_train_set), y=trunc_train_label)

trunc_train_pred <- predict(small_rf, t(trunc_train_set))
confusionMatrix(trunc_train_pred, trunc_train_label)

trunc_test_pred  <- predict(small_rf, t(trunc_test_set))
confusionMatrix(trunc_test_pred,  trunc_test_label )
## As expected, using important features gives at least sufficiently enough result

#========================
# Graph construction
#========================
# covariance matrix

# Pearson correlation matrix; normalized by autocorrelation
pearson_cor <- function(arr) {
  cov_mat = cov(arr)
  sig_arr = data.matrix(sqrt(diag(cov_mat)))
  norm_cov_mat <- cov_mat / (sig_arr %*% t(sig_arr))
}


# min value of normalized covariance is 0.77756, retrieved by min(*)
# Could we make the estimator based on the correlation?

# Feature pruning is required to remove the neglectable features.
# but how?
# Covariance based methods are not possible; 10000 features could be used.

# Then, based on the std values, we could remove some features.

mu  <- apply(d_express_unique, 1, mean)
std <- apply(d_express_unique, 1, sd)

# TOP-N function
get_top_n <- function(target, top_n) {
  sorted_list = sort(target, decreasing=TRUE, index.return=TRUE)
  top_n_idx = head(sorted_list$ix, n = top_n)
  target[c(top_n_idx)]
}

top_std <- get_top_n(std, top_n = (2 * top_n))

d_express_std <- d_express_unique[c(names(top_std)), ]

## Features (for description)
# feature-1: top-n-id / features-2: features selected by std
features_1 <- as.character(top_n_id)
features_2 <- names(d_express_std)

diff_feature_2 <- features_2[!(features_1 %in% features_2)] # only in std based

## we could use the features above (std-based) in RF;
diff_rf <- randomForest(t(d_express_unique[, rnd_tr_idx]), trunc_train_label)
confusionMatrix(predict(diff_rf, t(d_express_unique[, rnd_ts_idx])), trunc_test_label)
# and results are similar.

## But comparing the histogram of Pearson correlation matrix,
# important-based selection gives more flat graph, while std-based filtering gives more pointy graph.
# hist(pearson_cor(t(d_express_top)))
# hist(pearson_cor(t(d_express_std)))
# Meaning, the std-based filtering might result in non-related important features.

# From this, it seems the Random Forest method are powerful


# ggm package
# need adjacency matrix
# glasso_ggm <- fitCovGraph(cov_mat)

# glasso package
inv_cov <- glasso(d_express, rho=0.1)  
# Matrix size of 19860 x 498 - cannot run the command with 16 Gb RAM

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


