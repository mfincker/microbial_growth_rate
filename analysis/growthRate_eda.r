library(ggplot2)
library(caret)
setwd('/Users/maeva/Documents/classes/fall 2016/CS 229/project')

#############
# GR labels #
#############
gr <- read.delim('./formattedData/gr_1206.levelled.data', header=F)
colnames(gr) <- c('growthRate')
summary(gr)

gr$growthRate <- factor(gr$growthRate)
summary(gr)
ggplot(data=gr, aes(x=growthRate)) + 
  geom_bar(stat="count", fill="#56B4E9", colour="black") + 
  stat_count(aes(label=..count..), vjust=2, geom="text", color="white")+
  xlab("incubation time (days)")

########
# Data #
########
data <- read.delim('./formattedData/pfam_16S_length_1206.filtered.data', header=F)

################
# Feature name #
################
features <- read.delim('./formattedData/features_1206.filtered.name', header= F)
colnames(data) <- t(features)

####################
# Adding log2(16S) #
####################
data$log16S <- log2(data$`16S`)

#############################
#############################
      # PCA analysis #
#############################

# data transformation:
trans = preProcess(data, c("BoxCox", "center", "scale"))
#trans = preProcess(data, c( "center", "scale"))

data.Trans = data.frame(
  trans = predict(trans, data))
PC = predict(trans, data)
data.pca = prcomp(data.Trans) 

plot(data.pca, type = "l")


data.pca.coord = as.data.frame(data.pca$x)
ggplot(data = data.pca.coord,
       aes(x = PC1,y = PC2)) + 
  geom_point(aes(colour = as.factor(gr$growthRate)))

###################################
# Feature selection by projection #
###################################

# To explain 95% of the variance we need 392 variables
# To explain 90% of the variance we need 299 variables
# To explain 80% of the variance we need 181 variables
trans2 = preProcess(data, 
                   method=c("BoxCox", "center", 
                            "scale", "pca"), thresh = 0.9)
PC2 = predict(trans2, data)

# Saving PCA projected data
write.table(PC2, 
            file='./formattedData/pfam_16S_length_1207.PCAprojected90.tsv', 
            quote=FALSE, sep='\t', col.names = F, row.names = F)

ggplot(data = PC2,
       aes(x = PC1,y = PC2)) + 
  geom_point(aes(colour = as.factor(gr$growthRate)))

biplot(data.pca)

# Important components of PC1 and PC2
sortedPC1.weights <- sort(abs(trans2$rotation[,1]), decreasing = T, index.return = T)
sortedPC1.weights <- sortedPC1.weights$ix[1:30]
features[sortedPC1.weights,]

sortedPC2.weights <- sort(abs(trans2$rotation[,2]), decreasing = T, index.return = T)
sortedPC2.weights <- sortedPC2.weights$ix[1:30]
features[sortedPC2.weights,]

##################################
##################################
# Feature selection with BioComb #
##################################

library("Biocomb")

data_labelled <- data
data_labelled[,ncol(data_labelled) + 1] <- as.factor(gr$growthRate)
colnames(data_labelled) <- c(colnames(data_labelled[,1:ncol(data_labelled)-1]), "label")

method="FastFilter"
disc<-"MDL"
thr=0
thr.cons=0.05
attrs.nominal=numeric()
max.f=300
features.fcbf = select.process(data_labelled,method=method,disc.method=disc,
                                   threshold=thr,attrs.nominal=attrs.nominal,
                                   max.no.features = max.f)

features[features.fcbf,]

method="auc"
features.auc = select.process(data_labelled,method=method,disc.method=disc,
                               threshold=thr,attrs.nominal=attrs.nominal,
                               max.no.features = max.f)

####################
## HUM froze
method="HUM"
features.HUM = select.process(data_labelled,method=method,disc.method=disc,
                                    threshold=thr,attrs.nominal=attrs.nominal,
                                    max.no.features = max.f)
####################

method="Chi-square"
features.chi2 = select.process(data_labelled,method=method,disc.method=disc,
                                    threshold=thr,attrs.nominal=attrs.nominal,
                                    max.no.features = max.f)

method="InformationGain"
features.infoGain = select.process(data_labelled,method=method,disc.method=disc,
                               threshold=thr,attrs.nominal=attrs.nominal,
                               max.no.features = max.f)

method="symmetrical.uncertainty"
features.symUncertainty = select.process(data_labelled,method=method,disc.method=disc,
                                   threshold=thr,attrs.nominal=attrs.nominal,
                                   max.no.features = max.f)

method="Relief"
features.relief = select.process(data_labelled,method=method,disc.method=disc,
                                         threshold=thr,attrs.nominal=attrs.nominal,
                                         max.no.features = max.f)

method="CFS"
features.CFS = select.process(data_labelled,method=method,disc.method=disc,
                                 threshold=thr,attrs.nominal=attrs.nominal,
                                 max.no.features = max.f)

method="CorrSF"
features.corrSF = select.process(data_labelled,method=method,disc.method=disc,
                              threshold=thr,attrs.nominal=attrs.nominal,
                              max.no.features = max.f)

method="Chi2-algorithm"
features.chi2Algo = select.process(data_labelled,method=method,disc.method=disc,
                                 threshold=thr,attrs.nominal=attrs.nominal,
                                 max.no.features = max.f)


######################
######################
# Softmax regression #
######################

library("softmaxreg2")
## split data into a train and test set
index     <- 1:nrow(data)
testindex <- sample(index, trunc(length(index)/10))
testset   <- data.Trans[testindex,]
testlabel <- gr$growthRate[testindex]
trainset  <- data.Trans[-testindex,]
trainlabel <- gr$growthRate[-testindex] 

# Softmax
softmax_model = softmaxReg(trainset[,features.fcbf], trainlabel, 
                          maxit = 100, type = "class",
                           algorithm = "sgd", rate = 0.05, batch = 20)

yFitMat = softmax_model$fitted.values
yFit = c()
for (i in 1:length(trainlabel)) {
  yFit = c(yFit, which(yFitMat[i,]==max(yFitMat[i,]))-1)
}
table(trainlabel, yFit)

yPred = predict.softmax(softmax_model, testset[,features.fcbf]) - 1
table(testlabel, yPred)

# Softmax with L2
softmax_model.l2 = softmaxReg(trainset, trainlabel, 
                           maxit = 5, type = "class",
                           algorithm = "sgd", rate = 0.05, batch = 50, L2 = TRUE)

yFitMat = softmax_model.l2$fitted.values
yFit = c()
for (i in 1:length(trainlabel)) {
  yFit = c(yFit, which(yFitMat[i,]==max(yFitMat[i,]))-1)
}
table(trainlabel, yFit)
yPred = predict(softmax_model.l2, testset) - 1
table(testlabel, yPred)

# Softmax with L1
softmax_model.l1 = softmaxReg(trainset, trainlabel, 
                              maxit = 50, type = "class",
                              algorithm = "sgd", rate = 0.05, batch = 50, 
                              L1 = TRUE)

yFitMat = softmax_model.l1$fitted.values
yFit = c()
for (i in 1:length(trainlabel)) {
  yFit = c(yFit, which(yFitMat[i,]==max(yFitMat[i,]))-1)
}
table(trainlabel, yFit)
yPred = predict(softmax_model.l1, testset) - 1
table(testlabel, yPred)

# Softmax with L1 and L2
softmax_model.l1l2 = softmaxReg(trainset, trainlabel, 
                              maxit = 50, type = "class",
                              algorithm = "sgd", rate = 0.05, batch = 50, 
                              L1 = TRUE, L2 = TRUE)

yFitMat = softmax_model.l1l2$fitted.values
yFit = c()
for (i in 1:length(trainlabel)) {
  yFit = c(yFit, which(yFitMat[i,]==max(yFitMat[i,]))-1)
}
table(trainlabel, yFit)
yPred = predict.softmax(softmax_model.l1l2, testset) - 1
table(testlabel, yPred)
# Same L1 and L2 with Adagrad
softmax_model.l1l2.ada = softmaxReg(trainset, trainlabel, 
                                maxit = 50, type = "class",
                                algorithm = "adagrad", rate = 0.05, batch = 50, 
                                L1 = TRUE, L2 = TRUE)

yFitMat = softmax_model.l1l2.ada@fitted.values
yFit.l1l2.ada = c()
for (i in 1:length(trainlabel)) {
  yFit.l1l2.ada = c(yFit.l1l2.ada, which(yFitMat[i,]==max(yFitMat[i,]))-1)
}
table(trainlabel, yFit.l1l2.ada)
yPred = predict.softmax(softmax_model.l1l2.ada, testset) - 1
table(testlabel, yPred)

# Different L2 and L1
softmax_model.l1l2.0.01 = softmaxReg(trainset, trainlabel, 
                                maxit = 50, type = "class",
                                algorithm = "sgd", rate = 0.05, batch = 50, 
                                L1 = TRUE, L2 = TRUE, penalty = 0.01,
                                penaltyL1 = 0.01)

yPred.0.01 = predict.softmax(softmax_model.l1l2.0.01, testset) - 1
table(testlabel, yPred.0.01)

