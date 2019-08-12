library(caret)
library(LiblineaR)
library(softmaxreg2)

# Workflow

setwd('~/Documents/classes/fall 2016/CS 229/project')

# Classifiers:
# - softmax with L1/L2 reg - optimize L1, L2
# - svm one vs one with RBF Kernel - optimize C, gamma
# - svm one vs all without Kernel, with L1 / L2 reg - optimize L1, L2, not necessary
# - random forest (reg via. tree depth) - optimize tree depth
# - simple nnet?

# Feature selection:
# - auc
# - fcbf
# - infoGain
# - symUncertainty
# - chi2
# - PCA projected
# - mrmr?
# Ven diagram of selected features?

# Comparison:
# - accuracy (vs. epoch)
# - confusion matrix (training / test)
# - precision / recall graph

# Probabilistic / weights interpretation

# PCA analysis

# What if I consider only slow (4, 5) / fast (0, 1, 2, 3)


# Utility functions
precision_recall <- function(yPred, yTrue, classes.num) {
  confusion = as.data.frame(table(pred = yPred, true = yTrue))
  p = 0
  r = 0
  for (i in 0:(classes.num -1)) {
    precisionData = confusion[confusion$pred == i,]
    if (dim(precisionData)[1] > 0) {
      precision.i = precisionData[precisionData$true == i,]
      if (dim(precision.i)[1]>0 && sum(precisionData$Freq) > 0){
        p = p + precisionData[precisionData$true == i,]$Freq / 
          sum(precisionData$Freq)
      }
    }
    
    recallData = confusion[confusion$true == i,]
    if (dim(recallData)[1] > 0) {
      recall.i = recallData[recallData$pred == i,]
      if (dim(recall.i)[1] > 0 && sum(recallData$Freq) > 0) {
        r = r + recallData[recallData$pred == i,]$Freq / sum(recallData$Freq)
      }
    }
  }
  p = p/ classes.num
  r = r/ classes.num
  return(c("p" = p, "r" = r))
}

# Loading data

## Labels
gr <- read.delim('./formattedData/gr_1206.levelled.data', header=F)
colnames(gr) <- c('label')
summary(gr)

## Data
data <- read.delim('./formattedData/pfam_16S_length_1206.filtered.data', header=F)

## Features
features <- read.delim('./formattedData/features_1206.filtered.name', header= F)
colnames(data) <- t(features)

## Data transformation
trans = preProcess(data, c("BoxCox", "center", "scale"))
data.trans = data.frame(trans = predict(trans, data))

## Labelled data
data.labelled = cbind(data, gr)
data.trans.labelled = cbind(data.trans, gr)

## Split test / train
index     <- 1:nrow(data)
testindex <- sample(index, trunc(length(index)/10))
testset.labelled   <- data.trans.labelled[testindex,]
trainset.labelled  <- data.trans.labelled[-testindex,]

testset <- data.trans[testindex,]
testlabel <- factors(gr$label[testindex])
trainset  <- data.trans[-testindex,]
trainlabel <- factors(gr$label[-testindex])

## Softmax optimization on full data
# Non regularized
softmax.noReg.fcbf = softmaxReg(trainset[,features.fcbf], trainlabel, 
                           maxit = 100, type = "class",
                           algorithm = "sgd", rate = 0.05, batch = 100)
pr.softmax.noReg.fcbf = precision_recall(yPred = predict.softmax(softmax.noReg.fcbf, testset[,features.fcbf]), 
                 yTrue = testlabel, 6)
# L1 regularized
tryTypes=c(6)
tryCosts=c(1000,10, 8, 4, 3, 2, 1, 0.6, 0.3, 0.1, 0.001)
bestCost=NA
bestBcr=0
bestType=NA
bestPr = c()

class.num = length(unique(trainlabel))
k = 10
bcr.list = c()
softmax.pr.df = data.frame("p" = double(), "r" = double(), 
                           "costL1" = double(), "costL2" = double(), "fs" = character())

softmax.pr.df$fs = "none"  
for(ty in tryTypes){
  for(co in tryCosts){
    # k-fold crossvalidation
    for (i in seq(0, k)) {
      i.min = max(1, i*nrow(trainset)%/%k)
      i.max = min((i+1)*nrow(trainset)%/%k -1, nrow(trainset))
      cat("i min: ",i.min,", i max: ",i.max,"\n")
      train.i = trainset[-seq(i.min, i.max),]
      train.label.i = trainlabel[-seq(i.min,i.max)]
      train.cv.i = trainset[seq(i.min, i.max),]
      train.cv.label.i = trainlabel[seq(i.min, i.max)]
      
      m=LiblineaR(data=train.i[,features.fcbf] ,
                  target=train.label.i,
                  type=ty,cost=co,verbose=FALSE)
      
      yPred = predict(m, train.cv.i[,features.fcbf])$predictions
      pr = precision_recall(yPred, train.cv.label.i, classes.num = class.num )
      print(pr)
      softmax.pr.df = rbind(softmax.pr.df, data.frame("p" = pr[1], "r" = pr[2], 
                                                      "costL1" = co, "costL2" = 0, "fs" = "fcbf"))
    }
    #bcr.mean = mean(bcr.list)
    cat("Results for C=",co," : ",
        "precision: ",mean(softmax.pr.df[softmax.pr.df$costL1 == co,]$p),
        "recall: ",mean(softmax.pr.df[softmax.pr.df$costL1 == co,]$r),sep="")
  }
}

softmax.pr.plot <- data.frame("p" = double(), "r" = double(), 
                      "costL1" = double(), "costL2" = double(), "fs" = character())
softmax.pr.plot = rbind(softmax.pr.plot, c(pr.softmax.noReg.fcbf, 0, 0, 0, 0, "fcbf"))
for (i in tryCosts) {
  softmax.pr.plot = rbind(softmax.pr.plot, 
                          c(colMeans(softmax.pr.df[softmax.pr.df$costL1 == i & softmax.pr.df$fs == "fcbf",c(1,2,3,4)]),
                            apply(softmax.pr.df[softmax.pr.df$costL1 == i & softmax.pr.df$fs == "fcbf",], 2, sd)[c(1,2)],
                          softmax.pr.df[softmax.pr.df$costL1 == i & softmax.pr.df$fs == "fcbf",c(5)][1]))
}
colnames(softmax.pr.plot) <- c("p", "r", "costL1", "costL2", "sd.p", "sd.r", "fs")
ggplot(data = softmax.pr.df, aes(x = p, y = r )) +
  geom_point(aes(colour = factor(costL1)))+
  xlab("precision") + ylab("recall") + 
  theme_minimal(base_size = 10) + theme(legend.position="bottom") +
  scale_colour_discrete(name ="L1 cost:")

# The best L1 cost term is 0.3 - 0.6

# Retraining on the full trainingset
softmax.L1.full = LiblineaR(data=trainset,
          target=trainlabel,
          type=6,cost=0.6,verbose=FALSE)

table(pred = predict(softmax.L1.full, testset)$predictions, true = testlabel)

yPred = predict(softmax.L1.full, testset)$predictions
precision_recall(predict(softmax.L1.full.fcbf, testset[,features.fcbf])$predictions, testlabel, classes.num = class.num)
print(pr.softmax.L1)

ggplot(data = softmax.pr.plot, aes(x = p, y = r)) +
  geom_point(aes(colour = factor(costL1), shape=factor(fs)), size = 4) + 
 ylab("recall") + xlab("precision")+
  theme_minimal(base_size = 15) + theme(legend.position="right") +
  scale_colour_discrete(name ="L1 cost:") +
  scale_shape(name = "feature selection") +
  labs(title = "Precision - recall of softmax classification\nwith L1 regularization")+
  theme(plot.title = element_text(hjust = 0.5))



  geom_point(aes(x = pr.softmax.noReg[1], y = pr.softmax.noReg[2]), colour = "black") +
    geom_point(aes(x = pr.softmax.L1[1], y = pr.softmax.L1[2]), colour = "red") 

length(softmax.L1.full$W[abs(softmax.L1.full$W) > 0.0000000000001])
# Weights of 850 different features are non zero over the 6 labels
# Label:          2   0   1   4   3   5 
# non 0 weights: 310 663 323 159  85  54
l0 <- softmax.L1.full$W["0",]
l0 <- names(l0[abs(l0) > 0.00001])
l0 <- l0[1:length(l0)-1]
l1 <- softmax.L1.full$W["1",]
l1 <- names(l1[abs(l1) > 0.00001])
l1 <- l1[1:length(l1)-1]
l2 <- softmax.L1.full$W["2",]
l2 <- names(l2[abs(l2) > 0.00001])
l2 <- l2[1:length(l2)-1]
l3 <- softmax.L1.full$W["3",]
l3 <- names(l3[abs(l3) > 0.00001])
l3 <- l3[1:length(l3)-1]
l4 <- softmax.L1.full$W["4",]
l4 <- names(l4[abs(l4) > 0.00001])
l4 <- l4[1:length(l4)-1]
l5 <- softmax.L1.full$W["5",]
l5 <- names(l5[abs(l5) > 0.00001])
l5 <- l5[1:length(l5)-1]


features.L1 <- union(union(union(union(l1, l2), l3), l4), l5)



#######################################################
#######################################################
#######################################################
#######################################################
## Softmax optimization on selected features
#######################################################

softmax.pr.plot.features = data.frame("p" = double(), "r" = double(), 
                                      "costL1" = double(), "fs" = character())

#######################################
## auc-selected features
# no regularization
softmax.noReg.auc = softmaxReg(trainset[,features.auc], trainlabel, 
                           maxit = 100, type = "class",
                           algorithm = "sgd", rate = 0.05, batch = 100)
pr = precision_recall(yPred = predict.softmax(softmax.noReg.auc, 
                                              testset[,features.auc]), 
                                    yTrue = testlabel, 6)
softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "auc", "costL1" = 0))

# l1 regularization
softmax.L1.auc=LiblineaR(data=trainset[,features.auc] ,
            target=trainlabel,
            type=6,cost=0.5,verbose=FALSE)

pr = precision_recall(yPred = predict(softmax.L1.auc, 
                                              testset[,features.auc])$predictions, 
                                              yTrue = testlabel, 6)

softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "auc", "costL1" = 0.5))


#######################################
## chi2-selected features

# no regularization
softmax.noReg.chi2 = softmaxReg(trainset[,features.chi2], trainlabel, 
                               maxit = 100, type = "class",
                               algorithm = "sgd", rate = 0.05, batch = 100)
pr = precision_recall(yPred = predict.softmax(softmax.noReg.chi2, 
                                              testset[,features.chi2]), 
                      yTrue = testlabel, 6)
softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "chi2", "costL1" = 0))

# l1 regularization
softmax.L1.chi2=LiblineaR(data=trainset[,features.chi2] ,
                         target=trainlabel,
                         type=6,cost=0.5,verbose=FALSE)

pr = precision_recall(yPred = predict(softmax.L1.chi2, 
                                      testset[,features.chi2])$predictions, 
                      yTrue = testlabel, 6)

softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "chi2", "costL1" = 0.5))


#######################################
## fcbf-selected features
# no regularization
softmax.noReg.fcbf = softmaxReg(trainset[,features.fcbf], trainlabel, 
                                maxit = 100, type = "class",
                                algorithm = "sgd", rate = 0.05, batch = 100)
pr = precision_recall(yPred = predict.softmax(softmax.noReg.fcbf, 
                                              testset[,features.]), 
                      yTrue = testlabel, 6)
softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "fcbf", "costL1" = 0))

# l1 regularization
softmax.L1.fcbf=LiblineaR(data=trainset[,features.fcbf] ,
                          target=trainlabel,
                          type=6,cost=0.5,verbose=FALSE)

pr = precision_recall(yPred = predict(softmax.L1.fcbf, 
                                      testset[,features.fcbf])$predictions, 
                      yTrue = testlabel, 6)

softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "fcbf", "costL1" = 0.5))

#######################################
## infoGain-selected features
# no regularization
softmax.noReg.infoGain = softmaxReg(trainset[,features.infoGain], trainlabel, 
                                maxit = 100, type = "class",
                                algorithm = "sgd", rate = 0.05, batch = 100)
pr = precision_recall(yPred = predict.softmax(softmax.noReg.infoGain, 
                                              testset[,features.]), 
                      yTrue = testlabel, 6)
softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "infoGain", "costL1" = 0))

# l1 regularization
softmax.L1.infoGain=LiblineaR(data=trainset[,features.infoGain] ,
                          target=trainlabel,
                          type=6,cost=0.5,verbose=FALSE)

pr = precision_recall(yPred = predict(softmax.L1.infoGain, 
                                      testset[,features.infoGain])$predictions, 
                      yTrue = testlabel, 6)

softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "infoGain", "costL1" = 0.5))


#######################################
## symUncertainty-selected features
# no regularization
softmax.noReg.symUncertainty = softmaxReg(trainset[,features.symUncertainty], trainlabel, 
                                    maxit = 100, type = "class",
                                    algorithm = "sgd", rate = 0.05, batch = 100)
pr = precision_recall(yPred = predict.softmax(softmax.noReg.symUncertainty, 
                                              testset[,features.]), 
                      yTrue = testlabel, 6)
softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "symUncertainty", "costL1" = 0))

# l1 regularization
softmax.L1.symUncertainty=LiblineaR(data=trainset[,features.symUncertainty] ,
                              target=trainlabel,
                              type=6,cost=0.5,verbose=FALSE)

pr = precision_recall(yPred = predict(softmax.L1.symUncertainty, 
                                      testset[,features.symUncertainty])$predictions, 
                      yTrue = testlabel, 6)

softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr[1], "r" = pr[2],
                                             "fs" = "symUncertainty", "costL1" = 0.5))



softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr.softmax.L1[1], 
                                             "r" = pr.softmax.L1[2],
                                             "fs" = "none", "costL1" = 0.5))

softmax.pr.plot.features <- rbind(softmax.pr.plot.features,
                                  data.frame("p" = pr.softmax.noReg[1], 
                                             "r" = pr.softmax.noReg[2],
                                             "fs" = "none", "costL1" = 0))


ggplot(data = softmax.pr.plot.features, aes(x = p, y = r)) +
  geom_jitter(aes(colour = factor(fs), shape=factor(costL1)), size = 3) +
  xlab("precision") + ylab("recall") +
  scale_color_discrete(name="Algorithm for feature selection:",
                       breaks=c("none", "fcbf", "auc", "chi2", "symUncertainty", "infoGain"),
                       labels=c("none", "FastFilter","AUC", "Chi2 Gain",
                                "Symmetric uncertainty", "Information Gain")) +
  scale_shape_discrete(name = "L1 regularization:",
                       breaks = c(0, 0.5),
                       labels= c("without", "with")) +
  scale_colour_brewer(palette="Set1")+
  ggtitle("Precision - recall plot of Softmax regression\ncomparing different feature selection algorithms\n and L1 regularization")

# Ven diagram of selected features
#function() {library(VennDiagram)
draw.quintuple.venn(area1 = 30, area2 = 30, area3 = 28, area4 = 30, area5 = 30,
                    n12 = length(intersect(features.auc, features.chi2)),
                    n13 = length(intersect(features.auc, features.fcbf)),
                    n14 = length(intersect(features.auc, features.infoGain)),
                    n15 = length(intersect(features.auc, features.symUncertainty)),
                    n23 = length(intersect(features.chi2, features.fcbf)),
                    n24 = length(intersect(features.chi2, features.infoGain)),
                    n25 = length(intersect(features.chi2, features.symUncertainty)),
                    n34 = length(intersect(features.fcbf, features.infoGain)),
                    n35 = length(intersect(features.fcbf, features.symUncertainty)),
                    n45 = length(intersect(features.infoGain, features.symUncertainty)),
                    n123 = length(intersect(features.auc, intersect(features.chi2, features.fcbf))),
                    n124 = length(intersect(features.auc, intersect(features.chi2, features.infoGain))),
                    n125 = length(intersect(features.auc, intersect(features.chi2, features.symUncertainty))),
                    n134 = length(intersect(features.auc, intersect(features.fcbf, features.infoGain))),
                    n135 = length(intersect(features.auc, intersect(features.fcbf, features.symUncertainty))),
                    n145 = length(intersect(features.auc, intersect(features.infoGain, features.symUncertainty))),
                    n234 = length(intersect(features.chi2, intersect(features.fcbf, features.infoGain))),
                    n235 = length(intersect(features.chi2, intersect(features.fcbf, features.symUncertainty))),
                    n245 = length(intersect(features.chi2, intersect(features.infoGain, features.symUncertainty))),
                    n345 = length(intersect(features.fcbf, intersect(features.infoGain, features.symUncertainty))),
                    n1234 = length(intersect(features.auc, intersect(features.chi2, intersect(features.fcbf, features.infoGain)))),
                    n1235 = length(intersect(features.auc, intersect(features.chi2, intersect(features.fcbf, features.symUncertainty)))),
                    n1245 = length(intersect(features.auc, intersect(features.chi2, intersect(features.infoGain, features.symUncertainty)))),
                    n1345 = length(intersect(features.auc, intersect(features.fcbf, intersect(features.infoGain, features.symUncertainty)))),
                    n2345 = length(intersect(features.chi2, intersect(features.fcbf, intersect(features.infoGain, features.symUncertainty)))),
                    n12345 = length(intersect(features.auc, intersect(features.chi2, intersect(features.fcbf, intersect(features.infoGain, features.symUncertainty))))),
                    category = c("AUC", "Chi2", "FastFilter", "Information Gain" ,
                                 "Symmetric Uncertainty"),
                    euler.d = TRUE
)}


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#################################
# SVM - one vs. one, RBF kernel #
#################################
library(e1071)
svm.pr = data.frame("p" = double(), "r" = double(), 
                    "fs" = character(), "kernel" = character(),
                    "gamma" = double(), "cost" = double(), "weights" = integer())

addToSvmPr <- function(model, fs = "none", class.num = 6, weights = TRUE){
  if (fs == "none") {
    t = testset
  } else if (fs == "fcbf") {
    t = testset[,features.fcbf]
  } else if (fs == "infoGain") {
    t = testset[,features.infoGain]
  } else if (fs == "L1") {
    t = testset[,features.L1]
  }
  
  pred <- predict(model, t)
  
  pr = precision_recall(yPred = pred, yTrue = testlabel, class.num)
  svm.pr <<- rbind(svm.pr, data.frame(
    "p" = pr[1], "r" = pr[2], "fs" = fs, "kernel" = model$kernel, 
    "gamma" = model$gamma, "cost" = model$cost, "weigths" = weights
  ))
}

########## 
# Linear kernel ######################################
# Full data
svm.linear.noWeights <- svm(trainset, as.factor(trainlabel), kernel = "linear")
addToSvmPr(model = svm.linear.noWeights, weights = 0)

# with weigths
wi = c("0" = 0.5, "1" = 1, "2" = 1.2, "3" = 3.4, "4" = 0.7, "5" = 5.7)
svm.linear <- svm(trainset, as.factor(trainlabel), kernel = "linear",class.weights = wi)
addToSvmPr(model = svm.linear)

####### 
# RBF kernel ######################################
# Full data
svm.RBF.noWeights <- svm(trainset, as.factor(trainlabel), kernel = "radial")
addToSvmPr(model = svm.RBF.noWeights, weights = 0)

  # with weigths
wi = c("0" = 0.5, "1" = 1, "2" = 1.2, "3" = 3.4, "4" = 0.7, "5" = 5.7)
svm.RBF <- svm(trainset, as.factor(trainlabel), kernel = "radial",class.weights = wi)
addToSvmPr(model = svm.RBF)

# Plot
ggplot(data = svm.pr[svm.pr$gamma < 2e-4,], aes(x = p, y = r)) +
  geom_jitter(aes(colour = factor(gamma), shape=factor(kernel)), size = 3)

# Optimizing for cost
for (c in c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) {
  cat("cost: ", c, "\n")
  m <- svm(trainset, as.factor(trainlabel), kernel = "radial",class.weights = wi,
           cost = c)
  print(table(pred = predict(m, testset), true = testlabel))
  #addToSvmPr(model = m)
}

# Optimizing for gamma
for (g in c(5000, 7000, 8000, 9000)) {
  cat("gamma: ", g, "\n")
  m <- svm(trainset, as.factor(trainlabel), kernel = "radial",class.weights = wi,
           cost = 100, gamma = 1/g)
  addToSvmPr(model = m)
}

# Best model
m <- svm(trainset, as.factor(trainlabel), kernel = "radial",class.weights = wi,
         cost = 7000, gamma = 0.03571428)
table(pred = predict(m, testset), true = testlabel)
precision_recall(predict(m, testset), testlabel, 6)


ggplot(svm.pr[c(1,2, 27, 28),]) + 
  geom_point(aes(x = p,y = r, colour = as.factor(svm.pr[c(1,2, 27, 28),]$fs),
                 shape = as.factor(svm.pr[c(1,2, 27, 28),]$weigths)), size = 4) +
  theme_minimal(base_size = 25) +
  xlab("precision") + ylab("recall") +
  scale_shape_discrete(name = "weighted",
                       breaks = c(0,1),
                       labels = c("no", "yes")) +
  scale_colour_discrete(name = "feature selection",
                       breaks = c("none","fcbf"),
                       labels = c("none", "fast filter"))


ggplot(svm.pr[svm.pr$fs == "none" & svm.pr$weigths == 1 & svm.pr$kernel==2 & abs(svm.pr$cost - 100) < 0.000001,]) + 
  geom_jitter(aes(x = p,y = r, 
                 colour = as.factor(svm.pr[svm.pr$fs == "none" & svm.pr$weigths == 1 & svm.pr$kernel==2 & abs(svm.pr$cost - 100) < 0.000001,]$gamma)), size = 4) +
  theme_minimal(base_size = 25) +
  xlab("precision") + ylab("recall") +
  scale_colour_discrete(name = "RBF gamma")


  

#######################################################
# FCBF features
svm.RBF.noWeights <- svm(trainset[,features.fcbf], as.factor(trainlabel), kernel = "radial")
addToSvmPr(model = svm.RBF.noWeights, weights = 0, fs = "fcbf")

# with weigths
wi = c("0" = 0.5, "1" = 1, "2" = 1.2, "3" = 3.4, "4" = 0.7, "5" = 5.7)
svm.RBF <- svm(trainset[,features.fcbf], as.factor(trainlabel), kernel = "radial",class.weights = wi)
addToSvmPr(model = svm.RBF, fs = "fcbf")

# Optimizing for cost
for (c in c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) {
  cat("cost: ", c, "\n")
  m <- svm(trainset[,features.fcbf], as.factor(trainlabel), kernel = "radial",class.weights = wi,
           cost = c)
  print(table(pred = predict(m, testset[,features.fcbf]),
              true = testlabel))
 # addToSvmPr(model = m, fs="fcbf")
}

# Optimizing for gamma
for (g in c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) {
  cat("gamma: ", g, "\n")
  m <- svm(trainset[,features.fcbf], as.factor(trainlabel), kernel = "radial",class.weights = wi,
           cost = 1000, gamma = 1/g)
  print(table(pred = predict(m, testset[,features.fcbf]),
              true = testlabel))
  #addToSvmPr(model = m, fs = "fcbf")
}

#####################
# Model comparison
svm_full <-  svm(trainset, as.factor(trainlabel), kernel = "radial",class.weights = wi, cost = 100)
svm_fcbf <-  svm(trainset[,features.fcbf], as.factor(trainlabel), kernel = "radial",class.weights = wi, cost = 1000)
softmax_full = LiblineaR(data=trainset,target=trainlabel,type=6,cost=0.6,verbose=FALSE)
softmax_fcbf = LiblineaR(data=trainset[,features.fcbf],target=trainlabel,type=6,cost=0.6,verbose=FALSE)
rf_full = rf
rf_fcbf = rf.fcbf.20





###################

# Plot
ggplot(data = svm.pr[svm.pr$fs == "fcbf",], aes(x = r, y = p)) +
  geom_jitter(aes(colour = factor(cost)), size = 3)

table(pred = predict(m, testset[,features.fcbf]), true = testlabel)

###################################
# InfoGain features
svm.RBF.noWeights <- svm(trainset[,features.infoGain], as.factor(trainlabel), kernel = "radial")
addToSvmPr(model = svm.RBF.noWeights, weights = 0, fs = "infoGain")

# with weigths
wi = c("0" = 0.5, "1" = 1, "2" = 1.2, "3" = 3.4, "4" = 0.7, "5" = 5.7)
svm.RBF <- svm(trainset[,features.infoGain], as.factor(trainlabel), kernel = "radial",class.weights = wi)
addToSvmPr(model = svm.RBF, fs = "infoGain")

# Optimizing for cost
for (c in c(600, 800, 1200, 1400)) {
  cat("cost: ", c, "\n")
  m <- svm(trainset[,features.infoGain], as.factor(trainlabel), kernel = "radial",class.weights = wi,
           cost = c)
  print(table(pred = predict(m, testset[,features.infoGain]),
              true = testlabel))
  #addToSvmPr(model = m, fs="infoGain")
}

# Optimizing for gamma
for (g in c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) {
  cat("gamma: ", g, "\n")
  m <- svm(trainset[,features.infoGain], as.factor(trainlabel), kernel = "radial",class.weights = wi,
           cost = 1000, gamma = 1/g)
  print(table(pred = predict(m, testset[,features.infoGain]),
        true = testlabel))
  #addToSvmPr(model = m, fs = "infoGain")
}

# Plot
ggplot(data = svm.pr[svm.pr$fs == "infoGain" ,], aes(x = r, y = p)) +
  geom_jitter(aes(colour = factor(gamma)), size = 3)

table(pred = predict(m, testset[,features.infoGain]), true = testlabel)

########################################################
# Features selected via L1-softmax
svm.RBF.noWeights <- svm(trainset[,features.L1], as.factor(trainlabel), kernel = "radial")
cat("No weights")
print(table(pred = predict(svm.RBF.noWeights, testset[,features.L1]),
            true = testlabel))
addToSvmPr(model = svm.RBF.noWeights, weights = 0, fs = "L1")

# with weigths
wi = c("0" = 0.5, "1" = 1, "2" = 1.2, "3" = 3.4, "4" = 0.7, "5" = 5.7)
svm.RBF <- svm(trainset[,features.L1], as.factor(trainlabel), kernel = "radial",class.weights = wi)
addToSvmPr(model = svm.RBF, fs = "L1")
cat("default, with weights")
print(table(pred = predict(svm.RBF, testset[,features.L1]),
            true = testlabel))

ggplot(data= svm.pr[svm.pr$fs == "L1",], aes(x = p, y = r)) + geom_point(aes(colour = factor(gamma)))

# Optimizing for cost
for (c in c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) {
  cat("cost: ", c, "\n")
  m <- svm(trainset[,features.L1], as.factor(trainlabel), kernel = "radial",class.weights = wi,
           cost = c)
  print(table(pred = predict(m, testset[,features.L1]),
              true = testlabel))
  addToSvmPr(model = m, fs="L1")
}

# Optimizing for gamma
for (g in c(800, 850, 900)) {
  cat("gamma: ", g, "\n")
  m <- svm(trainset[,features.L1], as.factor(trainlabel), kernel = "radial",class.weights = wi,
           cost = 10000, gamma = 1/g)
  print(table(pred = predict(m, testset[,features.L1]),
              true = testlabel))
  addToSvmPr(model = m, fs = "L1")
}


###########################################
###########################################
###########################################
###########################################
# Random forest
###############
rf.pr <- data.frame("p" = double(), "r" = double(), 
                    "mtry" = double(), "fs" = character())

pr.rf <- precision_recall(yPred = predict(rf, testset),
                          yTrue = testlabel, 6)
pr.rf.fcbf <- precision_recall(yPred = predict(rf.fcbf, testset),
                          yTrue = testlabel, 6)

rf <- randomForest(trainset, as.factor(trainlabel), importance=TRUE, ntree = 2500)
table(pred = predict(rf, testset), true = testlabel)

rf.fcbf<- randomForest(trainset[,features.fcbf], as.factor(trainlabel), importance=TRUE, ntree = 2500)
table(pred = predict(rf.fcbf, testset[,features.fcbf]), true = testlabel)

rf.10 <- randomForest(trainset, as.factor(trainlabel), importance=TRUE, ntree = 2500, mtry = 10)
table(pred = predict(rf.10, testset), true = testlabel)

rf.fcbf.10<- randomForest(trainset[,features.fcbf], as.factor(trainlabel), importance=TRUE,mtry = 10, ntree = 2500)
table(pred = predict(rf.fcbf.10, testset[,features.fcbf]), true = testlabel)

pr.rf.10 <- precision_recall(yPred = predict(rf.10, testset),
                          yTrue = testlabel, 6)
pr.rf.fcbf.19 <- precision_recall(yPred = predict(rf.fcbf.10, testset),
                               yTrue = testlabel, 6)

rf.1 <- randomForest(trainset, as.factor(trainlabel), importance=TRUE, ntree = 2500, mtry = 1)
table(pred = predict(rf.1, testset), true = testlabel)

rf.fcbf.1<- randomForest(trainset[,features.fcbf], as.factor(trainlabel), importance=TRUE,mtry = 1, ntree = 2500)
table(pred = predict(rf.fcbf.1, testset[,features.fcbf]), true = testlabel)

pr.rf.1 <- precision_recall(yPred = predict(rf.1, testset),
                             yTrue = testlabel, 6)
pr.rf.fcbf.1 <- precision_recall(yPred = predict(rf.fcbf.1, testset),
                                  yTrue = testlabel, 6)

rf.20 <- randomForest(trainset, as.factor(trainlabel), importance=TRUE, ntree = 2500, mtry = 20)
table(pred = predict(rf.20, testset), true = testlabel)

rf.fcbf.20<- randomForest(trainset[,features.fcbf], as.factor(trainlabel), importance=TRUE,mtry = 20, ntree = 2500)
table(pred = predict(rf.fcbf.20, testset[,features.fcbf]), true = testlabel)

pr.rf.20 <- precision_recall(yPred = predict(rf.20, testset),
                             yTrue = testlabel, 6)
pr.rf.fcbf.20 <- precision_recall(yPred = predict(rf.fcbf.20, testset),
                                  yTrue = testlabel, 6)


rf.pr <- rbind(rf.pr, data.frame("p" = pr.rf.1[1], "r" = pr.rf.1[2], "mtry" = 1, "fs" = "none"))
rf.pr <- rbind(rf.pr, data.frame("p" = pr.rf.fcbf.1[1], "r" = pr.rf.fcbf.1[2], "mtry" = 1, "fs" = "fcbf"))
rf.pr <- rbind(rf.pr, data.frame("p" = pr.rf[1], "r" = pr.rf[2], "mtry" = 5, "fs" = "none"))
rf.pr <- rbind(rf.pr, data.frame("p" = pr.rf.fcbf[1], "r" = pr.rf.fcbf[2], "mtry" = 5, "fs" = "fcbf"))
rf.pr <- rbind(rf.pr, data.frame("p" = pr.rf.10[1], "r" = pr.rf.10[2], "mtry" = 10, "fs" = "none"))
rf.pr <- rbind(rf.pr, data.frame("p" = pr.rf.fcbf.19[1], "r" = pr.rf.fcbf.19[2], "mtry" = 10, "fs" = "fcbf"))
rf.pr <- rbind(rf.pr, data.frame("p" = pr.rf.20[1], "r" = pr.rf.20[2], "mtry" = 20, "fs" = "none"))
rf.pr <- rbind(rf.pr, data.frame("p" = pr.rf.fcbf.20[1], "r" = pr.rf.fcbf.20[2], "mtry" = 20, "fs" = "fcbf"))

ggplot(rf.pr) + 
  geom_jitter(aes(x = p,y = r, 
                  colour = as.factor(rf.pr$mtry), shape = as.factor(rf.pr$fs)), size = 4) +
  theme_minimal(base_size = 25) +
  xlab("precision") + ylab("recall") +
  scale_colour_discrete(name = "Number of variable\nsampled") +
  scale_shape(name = "Feature selection",
              breaks = c("none", "fcbf"),
              labels = c("none", "fast filter"))



