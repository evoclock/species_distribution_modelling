
# Dealing with predictor variables
# Environmental variables are typically organised as raster (grid) type files.
# Each predictor should have an associated raster that represents the variable
# of interest.
# Examples of predictor variables include:
# - climatic
# - soil
# - terrain
# - vegetation
# - land use
# etc
# Avoid ASCII files if possible (mainly due to performance)
# The layers for each analysis should have the same spatial extent, resolution, 
# origin, and projection for consistency of course.


# Using the dismo package we can create a raster stack from the available files

path = file.path(system.file(package = "dismo"), 'ex')

# Fetch the names of all files ending with the ".grd" extension in that folder

library(dismo)
files = list.files(path, pattern='grd$', full.names=TRUE)
files

# Create the raster stack of predictor variables

predictors = stack(files)
names(predictors)
plot(predictors)

# Alternatively, make a plot of a single raster and plot additional data on top
# of it.

# Will try to stick with the geodata package since maptools is being removed 
# from CRAN

library(geodata)

# get the world boundaries
world = world(resolution = 1, level = 0, path = "data/", version="latest")

# get the bradypus data for this example
file = paste(system.file(package = "dismo"), "/ex/bradypus.csv", sep = "")

bradypus = read.table(file, header = TRUE, sep = ',')

# we do not need the species column

bradypus = bradypus[, -1]

# Plot the first layer of the raster stack for example

plot(predictors, 1)
plot(world, add = TRUE)
points(bradypus, col='blue')

library(maptools)
data(wrld_simpl)

plot(predictors, 1)
plot(wrld_simpl, add = TRUE)
points(bradypus, col='blue')

# maptools still produces a better output

#######

# How to extract values from these raster layers

# extract value of predictor variables for points where occurrence of
# bradypus is noted in the dataset
presvals = extract(predictors, bradypus)

set.seed(0)
# extract value of predictor variables for 500 random background points
backgr = randomPoints(predictors, 500)
absvals = extract(predictors, backgr)

# combine into a dataframe where pb indicates whether
# the point is a background or presence point

pb = c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata = data.frame(cbind(pb, rbind(presvals, absvals)))

# Encode the biome column as a factor (categorical data) to avoid it being 
# treated as a numerical variable
sdmdata[,'biome'] = as.factor(sdmdata[,'biome'])

summary(sdmdata)

# Good suggestions are to extract multiple points in a radius as a way to deal 
# with mismatch between location accuracy and grid cell size. 
# If you want to explore the effect of uncertainty in a given location, 10 
# datasets could be made that represent 10 equally valid samples of the 
# environment in that radius, to which you could then fit 10 models to gauge
# uncertainty at that location.

library(corrplot)
corrplot(cor(sdmdata[,2:5]))

pairs(sdmdata[,2:5], cex = 0.1)

saveRDS(sdmdata, "outputs/sdm.Rds")
saveRDS(presvals, "outputs/pvals.Rds")
