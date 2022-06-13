pack_s<-list("BiocManager", "psych", "caret", "randomForest", "ggm",
               "glasso","plot.matrix")
bpack_s<-list("GEOquery","graph","csdR")
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
library(glasso)
library(plot.matrix)
library(csdR)
library(magrittr)
library(igraph)
library(glue)
library(dplyr)
set.seed(45394534)
setwd('/Users/suminlee/Desktop/Probability_Graph_Model/term_project/PGM_project')

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
         length(unique(d_gene$GeneSymbol)))
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
#confusionMatrix(as.factor(mycn),as.factor(risk))

# duplication example "CUL4A"
cul4a <- d_express[d_gene$GeneSymbol == "CUL4A",]
sprintf("gene duplication example, CUL4A")
cor(t(cul4a))

#========================
# train-test set split for
# Random Forest and GGM
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

#========================
# MYCN performance 
# on test set
#========================
sprintf("MYCN amplification: test set")
test_mycn <- c(mycn[risk==1][(high_train_size+1):sum(risk==1)], 
               mycn[risk==0][(low_train_size+1):sum(risk==0)])
confusionMatrix(as.factor(test_mycn),as.factor(test_risk))


#========================
# high risk prediction
# with Random forest
#========================
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
confusionMatrix(as.factor(test_pred), as.factor(test_risk))

impt <- importance(rf)
sorted_gene <- sort(impt/sum(impt),decreasing=TRUE,index.return=TRUE)
# for Diff.Net.Anal. top 1000
top_n <- 1000
sum(head(sort(impt/sum(impt),decreasing=TRUE), n = top_n))
top_n_id <- head(sorted_gene$ix, n= top_n)
top_n_importance <- head(sorted_gene$x, n= top_n)
# for GGMpred. top 140
top_n_min <- min(high_train_size, low_train_size)
sum(head(sort(impt/sum(impt),decreasing=TRUE), n = top_n_min))
top_n_min_id <- head(sorted_gene$ix, n= top_n_min)
top_n_min_importance <- head(sorted_gene$x, n= top_n_min)

#========================
# Graph construction
#========================
# top 1000 genes
#========================
#d_express_topn <- d_express_unique[c(top_n_id),]

# covariance matrix
#topn_cov_mat <- cov(t(d_express_topn))

# glasso package
#if (file.exists("topn_inv_cov.RData")){
#  load("topn_inv_cov.RData")
#} else {
#  topn_inv_cov <- glasso(topn_cov_mat, rho=0.1, nobs=length(risk))
#  save(topn_inv_cov,file="topn_inv_cov.RData")
#}
# visualize inverse covariance matrix's zero pattern up to top100 genes
# more than 400 dimension too mash
#jpeg(file='topn_inv_cov.jpeg', width=1000, height=1000, res=200)
#plot(topn_inv_cov$wi[1:200,1:200], border=NA) 
#dev.off()


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
# top 140 genes
#========================
# construct inverse covariance matrix
#========================
# low-risk group
# top 1000
#top1000_lowrisk <- low_risk_express[c(top_n_id),1: low_train_size]
#top1000_lowrisk_cov <- cov(t(top1000_lowrisk))
#if (file.exists("top1000_lowrisk_inv_cov.RData")){
#  load("top1000_lowrisk_inv_cov.RData")
#} else {
#  top1000_lowrisk_inv_cov <- glasso(top1000_lowrisk_cov, rho=0.1, nobs=low_train_size);
#  save(top1000_lowrisk_inv_cov, file="top1000_lowrisk_inv_cov.RData");
#}
#print(top1000_lowrisk_inv_cov$loglik)

# top 140
top140_lowrisk <- low_risk_express[c(top_n_min_id),1: low_train_size]
top140_lowrisk_cov <- cov(t(top140_lowrisk))
low_mu <- attr(scale(t(top140_lowrisk)), "scaled:center")
if (file.exists("top140_lowrisk_inv_cov.RData")){
  load("top140_lowrisk_inv_cov.RData")
} else {
  top140_lowrisk_inv_cov <- glasso(top140_lowrisk_cov, rho=0.1, nobs=low_train_size);
  save(top140_lowrisk_inv_cov, file="top140_lowrisk_inv_cov.RData");
}
print(top140_lowrisk_inv_cov$loglik)

# high-risk group
# top 1000
#top1000_highrisk <- high_risk_express[c(top_n_id),1:high_train_size]
#top1000_highrisk_cov <- cov(t(top1000_highrisk))
#if (file.exists("top1000_highrisk_inv_cov.RData")){
#  load("top1000_highrisk_inv_cov.RData")
#} else {
#  top1000_highrisk_inv_cov <- glasso(top1000_highrisk_cov, rho=0.1, nobs=high_train_size);
#  save(top1000_highrisk_inv_cov, file="top1000_highrisk_inv_cov.RData");
#}
#print(top1000_highrisk_inv_cov$loglik)

# top 140
top140_highrisk <- high_risk_express[c(top_n_min_id),1:high_train_size]
top140_highrisk_cov <- cov(t(top140_highrisk))
high_mu <-attr(scale(t(top140_highrisk)), "scaled:center")
if (file.exists("top140_highrisk_inv_cov.RData")){
  load("top140_highrisk_inv_cov.RData")
} else {
  top140_highrisk_inv_cov <- glasso(top140_highrisk_cov, rho=0.1, nobs=high_train_size);
  save(top140_highrisk_inv_cov, file="top140_highrisk_inv_cov.RData");
}
print(top140_highrisk_inv_cov$loglik)

#========================
# GGM prediction
# posterior P(G|x)
# GGM assumes multivariate Gaussian
#========================
multiGaussian <- function(x, mu, inv_cov){
  dimension <- dim(inv_cov)[1]
  sqrt_det <- sqrt(det(inv_cov))
  result <- sqrt_det*exp(-0.5*((x-mu)%*%inv_cov%*%(x-mu)))/(sqrt(2*pi)**dimension)
  # high dimensional problem...
  #result <- abs((-0.5*((x-mu)%*%inv_cov%*%(x-mu))) - log(sqrt_det) - log(sqrt(2*pi)**dimension))
  return(result)
}

GGMpred <- function(x){
  # probability of high risk given x
  # P(high), P(low)
  p_high <- high_train_size/(high_train_size+low_train_size)
  post_high <- multiGaussian(x, high_mu, top140_highrisk_inv_cov$wi)*p_high
  
  p_low <- low_train_size/(high_train_size+low_train_size)
  post_low <- multiGaussian(x, low_mu, top140_lowrisk_inv_cov$wi)*p_low
  #print(paste('post_high:',format(post_high),', post_low:',format(post_low)))
  return(post_high/(post_high+post_low))
}

#for test sample 1, 
#post_high = 1.449216e-84, post_low = 1.89427e-142
# return 1

top140_test_set <- test_set[c(top_n_min_id),]
ggm_pred_prob <- c(1:ncol(top140_test_set))
for (i in 1:ncol(top140_test_set)){
  ggm_pred_prob[i] <-GGMpred(top140_test_set[,i])
}
jpeg(file='ggm_pred_hist.jpeg', width=1000, height=1000, res=200)
hist(ggm_pred_prob)
dev.off()

ggm_pred <- c(1:length(ggm_pred_prob))
for (i in 1:length(ggm_pred_prob)){
  if (ggm_pred_prob[i] > 0.5){
    ggm_pred[i] <- 1
  } else {
    ggm_pred[i] <- 0
  }
}

sprintf("GGM test set prediction performance")
confusionMatrix(as.factor(ggm_pred), as.factor(test_risk))

# GGM with train set
top140_train_set <- train_set[c(top_n_min_id),]
ggm_train_pred_prob <- c(1:ncol(top140_train_set))
for (i in 1:ncol(top140_train_set)){
  ggm_train_pred_prob[i] <-GGMpred(top140_train_set[,i])
}
ggm_train_pred <- c(1:length(ggm_train_pred_prob))
for (i in 1:length(ggm_train_pred_prob)){
  if (ggm_train_pred_prob[i] > 0.5){
    ggm_train_pred[i] <- 1
  } else {
    ggm_train_pred[i] <- 0
  }
}

sprintf("GGM train set prediction performance")
confusionMatrix(as.factor(ggm_train_pred), as.factor(train_risk))

#========================
# expected calibration error
#========================
ECE <- function(pred, label, file_name){
  n_bin = 10
  h <- hist(pred)
  N <- length(pred)
  count <- h$count
  breaks <- h$breaks
  breaks[length(breaks)] <- breaks[length(breaks)] + 0.1 
  pr <- c(0:(n_bin-1))/n_bin + 0.05
  bins <- c(1:n_bin)*0
  ece <- 0
  for (i in 1:n_bin) {
    cond <- (breaks[i] <= pred) & (pred < breaks[i+1])
    if (sum(cond) != 0){
      bins[i] <- mean(label[cond]);
      ece <- ece + abs(mean(pred[cond])-bins[i])*count[i]/N;
    }
  }
  
  # plot
  jpeg(file=paste(file_name,".jpeg",sep=""),width=1000, height=1000, res=200)
  # positive ratio
  plot(c(0,1),c(0,1),type='l',lty=3,col="blue",xlab="",ylab="")
  lines(pr,bins, type='b')
  title(main=paste(file_name,format(ece,digits=3)),
        xlab='predicted value',
        ylab='positive ratio')
  dev.off()
  return(ece)
}

## RF
rf_pred <- predict(rf, t(test_set))
rf_pred_prob <- predict(rf, t(test_set),"prob")[,2]

## histogram
jpeg(file='rf_pred_hist.jpeg', width=1000, height=1000, res=200)
ax = pretty(0:10,n=10)
rf_pred_hist <- hist(rf_pred_prob, breaks=ax, plot=FALSE)
rf_pred_hist_positive <- hist(rf_pred_prob[test_risk==1], breaks=ax, plot=FALSE)
plot(rf_pred_hist)
plot(rf_pred_hist_positive,add=TRUE)
legend("topright",c("non-HR","HR"),col=c("grey","darkgrey"),lwd=6)
dev.off()

rf_ece <- ECE(rf_pred_prob, test_risk, 'Random_Foreest_ece')
print(rf_ece)

## GGM
## histogram
jpeg(file='ggm_pred_hist.jpeg', width=1000, height=1000, res=200)
ax = pretty(0:10,n=10)
ggm_pred_hist <- hist(ggm_pred_prob, breaks=ax, plot=FALSE)
ggm_pred_hist_positive <- hist(ggm_pred_prob[test_risk==1], breaks=ax, plot=FALSE)
plot(ggm_pred_hist)
plot(ggm_pred_hist_positive,add=TRUE)
legend("topright",c("non-HR","HR"),col=c("grey","darkgrey"),lwd=6)
dev.off()

ggm_ece <- ECE(ggm_pred_prob, test_risk, 'GGM_ece')
print(ggm_ece)

#========================
# MYCN amplification with RF or GGM
# for more accurate and 
# reliable prediction
#========================
# first filtering MYCN - high specificity
non_mycn_risk <- test_risk[test_mycn == 0]

# RF
non_mycn_test_set <- test_set[,test_mycn == 0]
non_mycn_rf_pred <- predict(rf, t(non_mycn_test_set))
non_mycn_rf_pred_prob <- predict(rf, t(non_mycn_test_set),"prob")[,2]
# accuracy
mean(non_mycn_rf_pred == non_mycn_risk)
rf_total_acc <- (sum(test_mycn==1) + sum(non_mycn_rf_pred==non_mycn_risk))/length(test_risk)
print(paste('non_mycn rf_total_acc:',format(rf_total_acc)))
print(paste('non_mycn RF correct sum:',format(sum(non_mycn_rf_pred==non_mycn_risk)),
            '/', format(sum(test_mycn == 0))))
non_mycn_rf_ece <- ECE(non_mycn_rf_pred_prob, non_mycn_risk, 'RF_nonMYCN')
print(paste('non_mycn_rf_ece:',format(non_mycn_rf_ece)))
confusionMatrix(non_mycn_rf_pred, as.factor(non_mycn_risk))

# GGM
top140_non_mycn_set <- top140_test_set[,test_mycn==0]
non_mycn_ggm_pred_prob <- c(1:ncol(top140_non_mycn_set))
for (i in 1:ncol(top140_non_mycn_set)){
  non_mycn_ggm_pred_prob[i] <-GGMpred(top140_non_mycn_set[,i])
}
non_mycn_ggm_pred <- c(1:length(non_mycn_ggm_pred_prob))
for (i in 1:length(non_mycn_ggm_pred_prob)){
  if (non_mycn_ggm_pred_prob[i] > 0.5){
    non_mycn_ggm_pred[i] <- 1
  } else {
    non_mycn_ggm_pred[i] <- 0
  }
}
confusionMatrix(as.factor(non_mycn_ggm_pred), as.factor(non_mycn_risk))

# accuracy
mean(non_mycn_ggm_pred==non_mycn_risk)
ggm_total_acc <- (sum(test_mycn==1) + sum(non_mycn_ggm_pred==non_mycn_risk))/length(test_risk)
print(paste('non mycn ggm_total_acc:',format(ggm_total_acc)))
print(paste('non_mycn GGM correct sum:',format(sum(non_mycn_ggm_pred==non_mycn_risk)),
            '/', format(sum(test_mycn == 0))))
non_mycn_ggm_ece <- ECE(non_mycn_ggm_pred_prob, non_mycn_risk, 'GGM_nonMYCN')
print(paste('non_mycn_ggm_ece:',format(non_mycn_ggm_ece)))


#========================
# Diff. Net. Anal
#========================
# top 140
csd_result_140 <- run_csd(x_1=t(top140_highrisk), x_2=t(top140_lowrisk),
                      n_it=10, verbose=FALSE)

pairs_to_pick <- 100
c_filter <- partial_argsort(csd_result_140$cVal, pairs_to_pick)
c_frame <- csd_result_140[c_filter, ]
s_filter <- partial_argsort(csd_result_140$sVal, pairs_to_pick)
s_frame <- csd_result_140[s_filter, ]
d_filter <- partial_argsort(csd_result_140$dVal, pairs_to_pick)
d_frame <- csd_result_140[d_filter, ]

csd_filter <- c_filter %>%
  union(s_filter) %>%
  union(d_filter)
csd_frame <- csd_results_140[csd_filter, ]

c_network <- graph_from_data_frame(c_frame, directed = FALSE)
s_network <- graph_from_data_frame(s_frame, directed = FALSE)
d_network <- graph_from_data_frame(d_frame, directed = FALSE)
E(c_network)$edge_type <- "C"
E(s_network)$edge_type <- "S"
E(d_network)$edge_type <- "D"
combined_network <- igraph::union(c_network, s_network, d_network)
# Auxillary function for combining
# the attributes of the three networks in a proper way
join_attributes <- function(graph, attribute) {
  ifelse(
    test = is.na(edge_attr(graph, glue("{attribute}_1"))),
    yes = ifelse(
      test = is.na(edge_attr(graph, glue("{attribute}_2"))),
      yes = edge_attr(graph, glue("{attribute}_3")),
      no = edge_attr(graph, glue("{attribute}_2"))
    ),
    no = edge_attr(graph, glue("{attribute}_1"))
  )
}
E(combined_network)$edge_type <- join_attributes(combined_network, "edge_type")
layout <- layout_nicely(combined_network)
E(combined_network)$color <- recode(E(combined_network)$edge_type,
                                    C = "darkblue", S = "green", D = "darkred"
)
jpeg(file='Diff_Net_Anal.jpeg', width=1000, height=1000, res=200)
plot(combined_network, layout = layout,
     vertex.size = 3, edge.width = 2, vertex.label.cex = 0.001)
dev.off()

# find network's hub genes
d_network_t <- table(c(d_frame$Gene1, d_frame$Gene2))
hub_gene <- as.numeric(names(d_network_t)[which.max(d_network_t)])
print(d_gene$GeneSymbol[hub_gene])


d_genes <- as.numeric(union(d_frame$Gene1, d_frame$Gene2))
s_genes <- as.numeric(union(s_frame$Gene1, s_frame$Gene2))

# =======================
# TODO
#========================
# GGM using d_genes
#========================


# =======================
# TODO
#========================
# GGM using s_genes
#========================





