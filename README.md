# Prediction of Neuromuscular Junction inervation using machine learning
## What is this?
We have created machine learning models using measurements taken from Neuromuscular Junctions which aim to be able to classify them into healthy or denervating. While we have tried seperation into partially and fully degenerating, the performance decreases.

## How can I use it
Your data must contain the followring columns which will be used for classificatoin
```
Bounds0Red, Bounds1Red, Bounds2Red, CCDistRed, CCSize0Red, CCSize1Red, CCSize2Red, CCRed, FragmentationRed, Compactness0Red, Compactness1Red, MeanIntensityRed, ShapeFactor0Red, ShapeFactor1Red, ShapeFactor2Red, SkeletonLengthRed, SurfaceVolumeRatio0Red, SurfaceVolumeRatio1Red, Surface0Red, Surface1Red, Surface2Red, Surface3Red, SurfaceDil0Red, SurfaceDil1Red, Volume0Red, Volume1Red,  Bounds0Gr, Bounds1Gr, Bounds2Gr, CCDistGr, CCSize0Gr, CCSize1Gr,  CCSize2Gr, CCGr, Compactness0Gr, Compactness1Gr, FragmentationGr, ShapeFactor0Gr, ShapeFactor1Gr, ShapeFactor2Gr, SkeletonLengthGr, SurfaceVolumeRatio0Gr, SurfaceVolumeRatio1Gr, Surface0Gr, Surface1Gr, Surface2Gr, Surface3Gr, SurfaceDil0Gr, SurfaceDil1Gr, Volume0Gr, Volume1Gr, IoUGr, AveDist, ComDist, Coverage, HausDist, NIntersection, NUnion
```
If your data contains only these rows, please use either **NMJRFKfoldFit.rds** or **NMJRFFit.rds** for your predictions.

It may also contain: 
```
Muscle, Fenotype, Genotype, Maturity 
```
If it does contain these rows then please use either **NMJRFKfoldFitNormalParam.rds** or **NMJRFFitNormalParam.rds** for predictions. This provides minor benifits but is by no means necessary.

The models we have created recognise the following options for each of those columns
- Muscle: FDB, or Lumbricals
- Fenotype: Healthy, or Disease
- Genotype: WT, FUSD14, GARS, SOD1 G93A, or TDP43 M323K
- Maturity: Early-mature, Mid-mature, or Late-mature

If your variables don't use these values, please use the other dataset them out.

In order to use a file to predict your data, create a folder containing your data as a csv, and the predictive model you want to use. 

Create an R file in this folder, set the working directory to the file locatoin and run this script
```
MLModel <- readRDS("Model.rds")
YourData <- read_csv("YourData.csv")

PredictionValues <- predict(MLModel, YourData)
YourData$PredictionValues <- PredictionValues
```

Replace the file names with your file names, and set the relevant working directories.

The only packages required are caret and tidyverse
