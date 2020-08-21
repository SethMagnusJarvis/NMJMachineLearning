library(caret)
library(tidyverse)
library(skimr)
library(ggplot2)
library(RANN) 
library(plotROC)
library(MLeval)
set.seed(1337)

#########################################################################################################################
#########################################Initialise Data#################################################################
#########################################################################################################################

#Load in initial data
MLDataset <- read_csv("AlanMLDatabaseTrueFinal.csv")
#Rename All columns and remove some unwanted values
colnames(MLDataset) <- c("Sample", "Muscle", "Fenotype", "Genotype", "Timepoint", "Maturity",
                         "NMJCounting", "Bounds0Red", "Bounds1Red", "Bounds2Red", 
                         "CCDistRed", "CCSize0Red", "CCSize1Red", "CCSize2Red", 
                         "CCRed", "FragmentationRed", "Compactness0Red", "Compactness1Red", 
                         "MeanIntensityRed","ShapeFactor0Red", "ShapeFactor1Red", "ShapeFactor2Red", 
                         "SkeletonLengthRed", "SurfaceVolumeRatio0Red.", "SurfaceVolumeRatio1Red.", "Surface0Red", 
                         "Surface1Red", "Surface2Red", "Surface3Red", "SurfaceDil0Red", 
                         "SurfaceDil1Red", "Volume0Red", "Volume1Red",  "Bounds0Gr", 
                         "Bounds1Gr", "Bounds2Gr", "CCDistGr", "CCSize0Gr", 
                         "CCSize1Gr",  "CCSize2Gr", "CCGr", "Compactness0Gr",  
                         "Compactness1Gr", "FragmentationGr", "ShapeFactor0Gr", "ShapeFactor1Gr", 
                         "ShapeFactor2Gr", "SkeletonLengthGr", "SurfaceVolumeRatio0Gr", "SurfaceVolumeRatio1Gr", 
                         "Surface0Gr", "Surface1Gr", "Surface2Gr", "Surface3Gr", 
                         "SurfaceDil0Gr", "SurfaceDil1Gr", "Volume0Gr", "Volume1Gr", 
                         "IoUGr", "AveDist", "ComDist", "Coverage", 
                         "HausDist", "NIntersection", "NUnion")

MLDataset <- MLDataset[rowSums(is.na(MLDataset)) != ncol(MLDataset),] %>%
  select(-Sample, -Timepoint) %>% select( NMJCounting, everything()) %>%
  filter(Muscle != "MUSCLE", Fenotype != "FENOTYPE")

#Change Variables from Text so caret can parse them
MLDatasetChanged <- MLDataset
MLDatasetChanged$Muscle[MLDatasetChanged$Muscle == "FBP"] <- 0
MLDatasetChanged$Muscle[MLDatasetChanged$Muscle == "Lumbricals"] <- 1
MLDatasetChanged$Fenotype[MLDatasetChanged$Fenotype == "Healthy"] <- 0
MLDatasetChanged$Fenotype[MLDatasetChanged$Fenotype == "Disease"] <- 1
MLDatasetChanged$Genotype[MLDatasetChanged$Genotype == "WT"] <- 0
MLDatasetChanged$Genotype[MLDatasetChanged$Genotype == "FUSD14"] <- 1
MLDatasetChanged$Genotype[MLDatasetChanged$Genotype == "Gars"] <- 2
MLDatasetChanged$Genotype[MLDatasetChanged$Genotype == "GARS"] <- 2
MLDatasetChanged$Genotype[MLDatasetChanged$Genotype == "SOD1 G93A"] <- 3
MLDatasetChanged$Genotype[MLDatasetChanged$Genotype == "TDP43 M323K"] <- 4
MLDatasetChanged$Maturity[MLDatasetChanged$Maturity == "Early-mature"] <- 0
MLDatasetChanged$Maturity[MLDatasetChanged$Maturity == "Mid-mature"] <- 1
MLDatasetChanged$Maturity[MLDatasetChanged$Maturity == "Late-mature"] <- 2

MLDatasetChanged <- drop_na(MLDatasetChanged, NMJCounting)

#Drop average dist column and merge the two types of unhealthy NMJs
MLDatasetChangeNA <- MLDatasetChanged %>% select(-AveDist) %>%
  drop_na()

MLDatasetChangeNAMerge <- MLDatasetChangeNA %>% 
  mutate(NMJCounting = replace(NMJCounting, NMJCounting == "Denervated" | NMJCounting == "Partial" , "Degenerating")) %>%
  mutate(NMJCounting = replace(NMJCounting, NMJCounting == "Fully", "Healthy")) 

MLDatasetTrimmed <- MLDatasetChangeNAMerge

#########################################################################################################################
#########################################Create ML Models################################################################
#########################################################################################################################

#Create a ML dataset without k-fold validation
#Partition the data into training and test data
trainRowNumbers <- createDataPartition(MLDatasetTrimmed$NMJCounting, p=0.8, list=FALSE)

trainData <- MLDatasetTrimmed[trainRowNumbers,]

testData <- MLDatasetTrimmed[-trainRowNumbers,]

#Format data for filtering
x = select(trainData, -NMJCounting)
y = trainData$NMJCounting

trainData$NMJCounting <- y


#Fit the ML model with k-fold validation
RFFitNoKFold <- train(as.factor(NMJCounting) ~ ., 
                      data = trainData, 
                      method = "ranger")


RFPredNoKFold <- predict(RFFitNoKFold, testData) 
# compare predicted outcome and true outcome
confusionMatrix(RFPredNoKFold, as.factor(testData$NMJCounting))

save(RFFitNoKFold, file="NMJRFFit")

#########################################################################################################################

#Create a ML dataset with k-fold validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)


RFFitKFold <- train(as.factor(NMJCounting) ~ ., 
                    data = trainData, 
                    method = "ranger",
                    trControl = fitControl)

RFPredKFold <- predict(RFFitKFold, testData) 
confusionMatrix(RFPredKFold, as.factor(testData$NMJCounting))

save(RFFitKFold, file="NMJRFKfoldFit")

#########################################################################################################################

# Create a ML dataset that saves predictions so AUC can be used 
trainData <- select(trainData, -ShapeFactor1Gr)
ctrl <- trainControl(method="cv", 
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)
rfFitWctrl <- train(NMJCounting ~ ., data=trainData, 
                    method="rf", preProc=c("center", "scale"), 
                    trControl=ctrl)

res <- evalm(rfFitWctrl,gnames='rf')

save(rfFitWctrl, file="NMJRFFitAUC")