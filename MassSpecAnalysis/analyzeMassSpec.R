## analyzeMassSpec - analyze a Mass Spectrometer Output This script
## (as opposed toa function) analyzes the mass spec. ALl input and
## output file names are formed from a unique ID, assumed passed as
## the first command line argument.
args <- commandArgs(TRUE)
ID <- args[1]

fileName = paste("MassSpec-",ID,".csv",sep="")

## These parameters seemed to work best
templateWidth <- 8
peakWidth <-  3
threshold <- 1.5
maxPeptideLength <- 16
  
## First, load the file, which is presumed to be a 2 column CSV
massSpec <- read.table(fileName,skip=5,sep=',')

## Convert that to a matrix because I hate typing $ all the time
## and I can't rely upon the column header being consistent
massSpecMatrix <- data.matrix(massSpec)

## Now, for some of these input files, there are multiple types of
## data in there, as in they give us some sort of smoothed or
## processed data after the stuff we really care about. We also
## don't care about masses below a certain level, so to reduce the
## size, we do some editing. First, get rid of those that have
## masses less than 500, which is to say less than a few amino acids
massSpecMatrix = massSpecMatrix[which(massSpecMatrix[,1] > 500),]

## Now find any that go backwards in mass
Temp = which(diff(massSpecMatrix[,1]) < 0)
if (length(Temp) > 1) {
    massSpecMatrix = massSpecMatrix[1:Temp[1],]
}

## Extract the dimensions thereof
matrixSize <- dim(massSpecMatrix)
numPoints <- matrixSize[1]

## We are now going to generate a normalizer and a limit. The
## normalizer is a split window probabilistic normalizer, which is
## basically the maximum intensity around a point excluding the
## area near the point, defined as the peakWidth. The limit is the
## maximum value in the peakWIdth area, and is used later to
## return only 1 peak within the peak width
normalizer <- mat.or.vec(numPoints,1)
peaks <- mat.or.vec(numPoints,1)

## extract the masses and intensities individually
masses <- massSpecMatrix[,1]
intensities <- massSpecMatrix[,2]

## Now, we have to go through data and compute those things
## Initialize the indices, including the length of the vector
maxIndex <- length(masses)
peakEnd <- 1
peakStart <- 1
templateEnd <- 1
templateStart <- 1
currentIndex <- 1
    
## Now, let's initialize those indices. First, find the first
## point that is more than the template width away from the
## first point, which is by definition the templateStart point
while(currentIndex < maxIndex &&
      masses[currentIndex] - masses[templateStart] < templateWidth)
    currentIndex = currentIndex + 1

## Now the peak peaks and the end of the template start from
## there.
peakEnd <- currentIndex
peakStart <- currentIndex
templateEnd <- currentIndex

## Now, go down from there to find the start of the peak
while(peakStart > 1 &&
      masses[currentIndex] - masses[peakStart] < peakWidth)
    peakStart = peakStart - 1

## and the end of the peak window and template window above that.
while(peakEnd <  maxIndex &&
      masses[peakEnd] - masses[currentIndex] < peakWidth)
    peakEnd= peakEnd + 1
while(templateEnd <  maxIndex &&
      masses[templateEnd] - masses[currentIndex] < templateWidth)
    templateEnd = templateEnd + 1

## Now, while there is room to move the template end upwards ...
while(templateEnd < maxIndex) {
    
    ## Compute the normalizer and the central peaks in the
    ## current window
    normalizer[currentIndex] =
        max(c(max(intensities[templateStart:(peakStart-1)]),
              max(intensities[(peakEnd+1):templateEnd])))
    peaks[currentIndex] = max(intensities[peakStart:peakEnd])
    
    ## Move up the current index ...
    currentIndex = currentIndex+1
    
    ## And like above, update the limits of the peak and templates
    while(peakStart < currentIndex &&
          masses[currentIndex] - masses[peakStart+1] > peakWidth)
        peakStart = peakStart + 1
    while(templateStart <  currentIndex &&
          masses[currentIndex]-masses[templateStart+1] > templateWidth)
        templateStart = templateStart + 1
    while(templateEnd <  maxIndex &&
          masses[templateEnd+1] > masses[templateEnd] &&
          masses[templateEnd] - masses[currentIndex] < templateWidth)
        templateEnd = templateEnd + 1
    while(peakEnd <  templateEnd &&
          masses[peakEnd] - masses[currentIndex] < peakWidth)
        peakEnd= peakEnd + 1
    
}
    
## Normalize and find the ones where the level is higher than the
## area around it by the threshold and where the peak is at the
## intensity of the point.
temp = intensities / normalizer
indices = which(temp > threshold & intensities == peaks & masses < 2000)

## Now plot the mass spec. First just where the peaks are and then
## the entire one so that we can see where it found the peaks.
png(paste("MassSpec-",ID,".png",sep=""));
limits = c(0,1.1*max(intensities))
plot(masses[indices],
     intensities[indices],
     main=paste('Mass Spec for ID',ID),
     xlab="Masses (da)",
     ylab="Intensities",ylim=limits)
lines(masses,intensities)

## Now, let's get the compositions by running the code. We had to
## put a link to the executable in a path that I could execute
## from. This is that path.
for (index in indices){
    print(paste("Test Index ",index," mass ",masses[index]))
    command <- paste("./computePeptideComposition ",
                     masses[index],
                     " ",
                     maxPeptideLength,
                     " Composition-",index,".csv",
                     sep='')
    print(paste("Excecute Command: ",command))
    system(command)
    composition <- read.csv(paste("Composition-",index,".csv",sep=''))
    print(paste("Found ",nrow(composition)," Compositions"))
    text(x=masses[peaks[index]],
         y=intensities[peaks[index]],
         paste(nrow(composition)),
         adj=c(1,0.5))
}

## This flushes the plotting to the file
dev.off()

