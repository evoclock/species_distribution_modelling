
# This part of the tutorial is rather basic and only covers glm and how to
# create a suitability score map with the predict function. It then goes into
# model evaluation techniques which is covered in a little more detail.


library(dismo)
sdmdata = readRDS("outputs/sdm.Rds")
presvals = readRDS("outputs/pvals.Rds")

# glm for some of the climatic variables
m1 = glm(pb ~ bio1 + bio5 + bio12 + bio16 + bio17, data = sdmdata)
# m1 = glm(pb ~ bio1 + bio5 + bio12, data = sdmdata)
summary(m1)

# glm for all variables
m2 = glm(pb ~ ., data = sdmdata)
summary(m2)


# Since models in dismo do not use a formula we use the matrix we built from
# the raster layer (presvals)

bc = bioclim(presvals[,c('bio1', 'bio5', 'bio12', 'bio16', 'bio17')])
# bc = bioclim(presvals[,c('bio1', 'bio5', 'bio12')])
class(bc)
bc
pairs(bc)


# Model prediction

bio1 = c(40, 150, 200)
bio5 = c(60, 115, 290)
bio12 = c(600, 1600, 1700)
bio16 = c(400, 800, 1200)
bio17 = c(50, 400, 500)

pd = data.frame(cbind(bio1, bio5, bio12, bio16, bio17))
# pd = data.frame(cbind(bio1, bio5, bio12))
pd

predict(m1, pd)

response(bc)

# Ideally, we always want to create a map of scores, for which we need to 
# provide the predict score with a raster object and model object.

predictors = stack(list.files(file.path(system.file(package="dismo"), 'ex'), 
                              pattern='grd$', full.names=TRUE ))
names(predictors)

p = predict(predictors, m1)

plot(p)


# Evaluation


# Does the model seem sensible, ecologically?
# Do the fitted functions (the shapes of the modeled relationships) make sense?
# Do the predictions seem reasonable? (map them, and think about them)?
# Are there any spatial patterns in model residuals? (see Leathwick and 
# Whitehead 2001 for an interesting example)

# In unbiased data, a high AUC indicates that sites with high predicted 
# suitability values tend to be areas of known presence and locations with 
# lower model prediction values tend to be areas where the species is not 
# known to be present (absent or a random point). An AUC score of 0.5 means 
# that the model is as good as a random guess. 

# Here we illustrate the computation of the correlation coefficient and 
# AUC with two random variables. p (presence) has higher values, and represents 
# the predicted value for 50 known cases (locations) where the species is 
# present, and a (absence) has lower values, and represents the predicted value 
# for 50 known cases (locations) where the species is absent.

p = rnorm(50, mean=0.7, sd=0.3)
a = rnorm(50, mean=0.4, sd=0.4)

par(mfrow=c(1, 2))
plot(sort(p), col='red', pch=21)
points(sort(a), col='blue', pch=24)
legend(1, 0.95 * max(a,p), c('presence', 'absence'),
       pch=c(21,24), col=c('red', 'blue'))
comb = c(p, a)
group = c(rep('presence', length(p)), rep('absence', length(a)))
boxplot(comb~group, col=c('blue', 'red'))

# Calculate the correlation coefficient 

grouped = c(rep(1, length(p)), rep(0, length(a)))
cor.test(comb, grouped)$estimate

# and AUC
mv = wilcox.test(p, a)
auc = as.numeric(mv$statistic) / (length(p) * length(a))
auc

# alternatively with the evaluate function from dismo

e = evaluate(p = p, a = a)
class(e)
e

par(mfrow = c(1, 2))
density(e)
boxplot(e, col = c('blue', 'red'))

# Eval on real data (presence-only). Divide the data into two random sets

# set the sampling criteria for the training dataset
samp = sample(nrow(sdmdata), round(0.75 * nrow(sdmdata)))
traindata = sdmdata[samp,]
traindata = traindata[traindata[,1] == 1, 2:9]

# Now take the complement for the evaluation dataset
testdata = sdmdata[-samp,]

# retrieve the rows that correspond to the training dataset from the bioclim 
# object
bc = bioclim(traindata)

e = evaluate(testdata[testdata==1,], testdata[testdata==0,], bc)
e

plot(e, 'ROC')

# with k-fold testing

pres = sdmdata[sdmdata[,1] == 1, 2:9]

# background data will only be used for model testing so no need to perform
# k-fold on it
back = sdmdata[sdmdata[,1] == 0, 2:9]

k = 5
group = kfold(pres, k)

# Now fit and test the model five times.

e = list()
for (i in 1:k) {
  # train with whichever group does not equal i
  train = pres[group != i,]
  # test with the group that equals i (i.e. fit the model with it)
  test = pres[group == i,]
  bc = bioclim(train)
  e[[i]] = evaluate(p = test, a = back, bc)
}

# retrieve AUC values, max of the sum of sensitivity (i.e. true positive rate),
# and the specificity (i.e. true negative rate). Threshold "spec_sens" is 
# sometimes used as a threshold for setting cells to presence or absence.

auc = sapply(e, function(x){x@auc})
mean(auc)

sapply( e, function(x){ threshold(x)['spec_sens']})

# Removing “spatial sorting bias” (the difference between the distance from 
# testing-presence to training-presence and the distance from testing-absence 
# to training-presence points) through “point-wise distance sampling”. To 
# address variance in AUC values dependent on the spatial extent used to select
# background points.

file = file.path(system.file(package="dismo"), "ex/bradypus.csv")
bradypus = read.table(file,  header=TRUE,  sep=',')
# species name not needed as before
bradypus = bradypus[,-1]
presvals = extract(predictors, bradypus)
set.seed(0)

backgr = randomPoints(predictors, 500)
nr = nrow(bradypus)
s = sample(nr, 0.25 * nr)
pres_train = bradypus[-s, ]
pres_test = bradypus[s, ]
nr = nrow(backgr)
set.seed(9)
s = sample(nr, 0.25 * nr)
back_train = backgr[-s, ]
back_test = backgr[s, ]

sb = ssb(pres_test, back_test, pres_train)
sb[,1] / sb[,2]

# sb[,1] / sb[,2] is an indicator of spatial sorting bias. In the absence of 
# SSB, this value will be 1, however, if the value is close to zero it indicates
# strong SSB.

# point wise sampling to remove SSB

i = pwdSample(pres_test, back_test, pres_train, n=1, tr=0.1)
pres_test_pwd = pres_test[!is.na(i[,1]), ]
back_test_pwd = back_test[na.omit(as.vector(i)), ]
sb2 = ssb(pres_test_pwd, back_test_pwd, pres_train)
sb2[1]/ sb2[2]

bc = bioclim(predictors, pres_train)
evaluate(bc, p=pres_test, a=back_test, x=predictors)
# With reduced SSB, AUC is also reduced
evaluate(bc, p=pres_test_pwd, a=back_test_pwd, x=predictors)


