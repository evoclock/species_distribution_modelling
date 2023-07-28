# Bioclim and Mahalanobis methods

# Recreate the data used prior
library(dismo)
library(maptools)
data(wrld_simpl)

predictors = stack(list.files(file.path(system.file(package="dismo"), 'ex'), pattern='grd$', full.names=TRUE ))
file = file.path(system.file(package="dismo"), "ex/bradypus.csv")

# Presence data
bradypus = read.table(file,  header=TRUE,  sep=',')
bradypus = bradypus[,-1]
presvals = extract(predictors, bradypus)

# Absence data
set.seed(0)
backgr = randomPoints(predictors, 500)
absvals = extract(predictors, backgr)
pb = c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))

sdmdata = data.frame(cbind(pb, rbind(presvals, absvals)))
sdmdata[,'biome'] = as.factor(sdmdata[,'biome'])

# Drop categorical variables if needed for a given model that cannot handle 
# these

pred_nf = dropLayer(predictors, 'biome')

# Split the presence data into training and test sets
set.seed(0)
group = kfold(bradypus, 5)
# Opting to make all groupds but group 1 the training set
pres_train = bradypus[group != 1, ]
pres_test = bradypus[group == 1, ]

# if desired for expediency, restrict predictions to a specific area

ext = extent(-90, -32, -33, 23)

# setting the background data for training and testing sets. 
# The first layer in the raster stack is used as a mask. This ensures that
# 1000 random points only occur within the spatial extent of the rasters, and 
# within  cells that do not report NAs, furthermore, there also should be a 
# single absence point per cell.
# The background points are further restricted to be within 12.5% of the 
# specified ext (extf = 1.25)

set.seed(10)
backg = randomPoints(pred_nf, n=1000, ext=ext, extf = 1.25)
colnames(backg) = c('lon', 'lat')

# kfold the background into 5 partitions
group = kfold(backg, 5)

# as with the presence data, all partitions but group 1 are the training set
backg_train = backg[group != 1, ]
backg_test = backg[group == 1, ]

# Plot the first raster and draw the extent area with absence and presence 
# points from respective test and training sets.
# Background training set = yellow -
# Background test set = black -
# Presence training set = green +
# Presence test set = blue +

r = raster(pred_nf, 1)
plot(!is.na(r), col=c('white', 'light grey'), legend=FALSE)
plot(ext, add=TRUE, col='red', lwd=2)
points(backg_train, pch='0', cex=0.5, col='yellow')
points(backg_test, pch='0',  cex=0.5, col='black')
points(pres_train, pch= '+', col='purple')
points(pres_test, pch='+', col='blue')

# Presence-only methods
# Bioclim and Mahalanobis
# Caveat: these do not perform well when considering climate change

# Use the pred_nf as we need the dataset without categorical variables
bc = bioclim(pred_nf, pres_train)
plot(bc, a=1, b=2, p=0.85)

e = evaluate(pres_test, backg_test, bc, pred_nf)
e
## class          : ModelEvaluation
## n presences    : 23
## n absences     : 200
## AUC            : 0.6890217
## cor            : 0.1765706
## max TPR+TNR at : 0.08592151

# find a threshold
tr = threshold(e, 'spec_sens')
tr
##  0.08592151

# Use the rasterstack with predictor variable to make a prediction to
# a raster layer for the area we set up with ext
# Use the predict function like this so that it ensures compatibility with
# models defined in all packages
pb = predict(pred_nf, bc, ext = ext, progress='')
pb

## class      : RasterLayer
## dimensions : 112, 116, 12992  (nrow, ncol, ncell)
## resolution : 0.5, 0.5  (x, y)
## extent     : -90, -32, -33, 23  (xmin, xmax, ymin, ymax)
## crs        : +proj=longlat +datum=WGS84 +no_defs
## source     : memory
## names      : layer
## values     : 0, 0.7096774  (min, max)

par(mfrow=c(1,2))
plot(pb, main='Bioclim, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(pb > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')

# Now Mahalanobis distance prediction

mm = mahal(pred_nf, pres_train)
e = evaluate(pres_test, backg_test, mm, pred_nf)
e
## class          : ModelEvaluation
## n presences    : 23
## n absences     : 200
## AUC            : 0.7686957
## cor            : 0.1506777
## max TPR+TNR at : 0.1116504

# It only takes into account correlations of the variables in the data set
# and it is independent of the scale of the measurements
pm = predict(pred_nf, mm, ext=ext, progress='')

par(mfrow=c(1,2))
pm[pm < -10] <- -10
plot(pm, main='Mahalanobis distance')
plot(wrld_simpl, add=TRUE, border='dark grey')

# as before find the threshold in the evaluation
tr = threshold(e, 'spec_sens')

plot(pm > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')

# Presence-absence models
# GLM & machine learning methods

# create a dataframe from the presence and background training sets
train = rbind(pres_train, backg_train)

# create a Boolean array matched to presence and background training sets 
pb_train = c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))

# extract environmental predictors from the rasterStack and training dataset
envtrain = extract(predictors, train)

# match the Boolean array to the presence absence training data and store in
# envtrain
envtrain = data.frame(cbind(pa = pb_train, envtrain) )

envtrain[,'biome'] = factor(envtrain[,'biome'], levels=1:14)

# extract environmental predictors from the rasterStack and presence 
# test dataset
testpres = data.frame(extract(predictors, pres_test) )

# extract environmental predictors from the rasterStack and absence test
# dataset. I don't know why Hijmans doesn't just create 
# test = rbind(pres_test, backg_test)
testbackg = data.frame(extract(predictors, backg_test) )

testpres[,'biome'] = factor(testpres[,'biome'], levels=1:14)
testbackg[,'biome'] = factor(testbackg[,'biome'], levels=1:14)

# logistic regression without interaction terms:
gm1 = glm(pa ~ bio1 + bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17,
           family = binomial(link = "logit"), data = envtrain)
summary(gm1)
coef(gm1)


gm2 = glm(pa ~ bio1+bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17,
           family = gaussian(link = "identity"), data = envtrain)
summary(gm2)
coef(gm2)

gm3 = glm(pa ~ bio1+bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17,
          family = poisson(link = "log"), data = envtrain)
summary(gm3)
coef(gm3)

# evaluate all three

ge1 = evaluate(testpres, testbackg, gm1)
ge1
## class          : ModelEvaluation
## n presences    : 23
## n absences     : 200
## AUC            : 0.8156522
## cor            : 0.308183
## max TPR+TNR at : -2.565312

ge2 = evaluate(testpres, testbackg, gm2)
ge2
# class          : ModelEvaluation 
# n presences    : 23 
# n absences     : 200 
# AUC            : 0.7886957 
# cor            : 0.3629842 
# max TPR+TNR at : 0.08711959

ge3 = evaluate(testpres, testbackg, gm3)
ge3
# class          : ModelEvaluation 
# n presences    : 23 
# n absences     : 200 
# AUC            : 0.8163043 
# cor            : 0.3031901 
# max TPR+TNR at : -2.606367


# Predict using all three

pg1 = predict(predictors, gm1, ext=ext)
par(mfrow=c(1,2))
plot(pg1, main='GLM/Binomial, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr = threshold(ge2, 'spec_sens')
plot(pg1 > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='*', cex = 0.75)

pg2 = predict(predictors, gm2, ext=ext)
par(mfrow=c(1,2))
plot(pg2, main='GLM/Gaussian, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr = threshold(ge2, 'spec_sens')
plot(pg2 > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='*', cex = 0.75)

pg3 = predict(predictors, gm3, ext=ext)
par(mfrow=c(1,2))
plot(pg3, main='GLM/Poisson, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr = threshold(ge2, 'spec_sens')
plot(pg3 > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='*', cex = 0.75)

# Machine learning methods
# MaxEnt

maxent()
xm = maxent(predictors, pres_train, factors='biome')
plot(xm)
response(xm)

# Evaluate

e = evaluate(pres_test, backg_test, xm, predictors)
e
## class          : ModelEvaluation
## n presences    : 23
## n absences     : 200
## AUC            : 0.8336957
## cor            : 0.3789954

# Predict
px = predict(predictors, xm, ext=ext, progress='')
par(mfrow=c(1,2))
plot(px, main='Maxent, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr = threshold(e, 'spec_sens')
plot(px > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='*', cex = 0.75)

# Random Forest

# The function randomForest can take a formula or, in two separate arguments, 
# a data.frame with the predictor variables, and a vector with the response. 
# If the response variable is a factor (categorical), randomForest will do 
# classification, otherwise it will do regression. Whereas with species 
# distribution modeling we are often interested in classification 
# (species is present or not), it is my experience that using regression 
# provides better results. rf1 does regression, rf2 and rf3 do classification 
# (they are exactly the same models). See the function tuneRF for optimizing 
# the model fitting procedure.

library(randomForest)
model = pa ~ bio1 + bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17
rf1 = randomForest(model, data=envtrain)

# If using rf2 remember to specify your response variable as a factor
model = factor(pa) ~ bio1 + bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17
rf2 = randomForest(model, data=envtrain)
rf2
# Call:
#   randomForest(formula = model, data = envtrain) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 2
# 
# OOB estimate of  error rate: 10.19%
# Confusion matrix:
#   0  1 class.error
# 0 774 26   0.0325000
# 1  65 28   0.6989247

rf3 = randomForest(envtrain[,1:8], factor(pb_train))
rf3
# Call:
#   randomForest(x = envtrain[, 1:8], y = factor(pb_train)) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 2
# 
# OOB estimate of  error rate: 0%
# Confusion matrix:
#   0  1 class.error
# 0 800  0           0
# 1   0 93           0

erf = evaluate(testpres, testbackg, rf1)
erf
## class          : ModelEvaluation
## n presences    : 23
## n absences     : 200
## AUC            : 0.8580435
## cor            : 0.5010053
## max TPR+TNR at : 0.1060667

# rf2 and rf3 cannot be evaluated like this

pr = predict(predictors, rf1, ext=ext)
par(mfrow=c(1,2))
plot(pr, main='Random Forest, regression')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr = threshold(erf, 'spec_sens')
plot(pr > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='*', cex=0.75)

# SVM

library(kernlab)
svm = ksvm(pa ~ bio1+bio5+bio6+bio7+bio8+bio12+bio16+bio17, data=envtrain)
esv = evaluate(testpres, testbackg, svm)
esv
## class          : ModelEvaluation
## n presences    : 23
## n absences     : 200
## AUC            : 0.7576087
## cor            : 0.3738667
## max TPR+TNR at : 0.02857293

# predict based on the SVM model
ps = predict(predictors, svm, ext=ext)
par(mfrow=c(1,2))
plot(ps, main='Support Vector Machine')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(esv, 'spec_sens')
plot(ps > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='*', cex=0.75)

# Combining model predictions using a rasterStack for each model

models = stack(pb, pm, pg1, pg2, pr, ps)
names(models) = c("bioclim", "mahalanobis", "glm binomial", "glm gaussian",
                  "rf", "svm")
plot(models)

# Do not generate a plot of the mean prediction without first making sure
# all models are between 0 and 1

# We will generate a weighted plot using AUC scores.
# substract 0.5 (the random expectation) and square the result to give further
# weight to higher AUC values

auc = sapply(list(ge1, ge2, erf, esv), function(x) x@auc)
w = (auc-0.5)^2
m2 = weighted.mean( models[[c("glm.binomial", "glm.gaussian","rf", "svm")]], w)
plot(m2, main='weighted mean of four models')
