## analyzeMassSpec - analyze a Mass Spectrometer Output. This script
## (as opposed to a function) analyzes the mass spec and can be called
## from the shell using Rscript to invoke it as follows:
##
## Rscript analyzeMassSpec 123456
##
## Where the first argument (in this case 123456) is a unique ID used
## to keep instantiations from stepping on each other.
##
## The input mass spec is expected to be in a file called
## MassSpec-ID.csv as a two-column .csv with no more than 5 header
## lines.
##
## It will plot the Mass Spec into a ping file called MassSpec-ID.png
## and identify up to 8 peaks to be analyzed. It will then, using the
## program computeParallelComposition, generate a file called
## Compositions-ID-[0-7].csv, each one containing the composition list
## for the identified peak weights.

## Get the command line arguments that matter
args <- commandArgs(TRUE)
ID <- args[1]

## Form the name of the input file
fileName = paste("./mic-data/",ID,"-data.csv",sep="")

## These parameters seemed to work best when it comes to identifying
## peaks inthe Mass Spec using a split window normalizer.
templateWidth <- 8
peakWidth <-  3
threshold <- 1.5
maxPeptideLength <- 16
  
## Load the file, which is presumed to be a 2 column CSV. Use
## read.table since it can deal with comments and take the "skip"
## argument.
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

## Now, we have to go through data and compute the split window
## normalizer. IN order to do this quickly, given that the masses are
## not uniformly separated, we set up the limits of the template and
## the "peak". The peak is the middle part of the normalizer that we
## don't include in the normalization.
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

## Now the peak limits and the end limit and the end of the template
## start from there.
peakEnd <- currentIndex
peakStart <- currentIndex
templateEnd <- currentIndex

## Now, go down from there to find the start of the peak
while(peakStart > 1 &&
      masses[currentIndex] - masses[peakStart] < peakWidth)
    peakStart <- peakStart - 1

## and the end of the peak window and template window above that.
while(peakEnd <  maxIndex &&
      masses[peakEnd] - masses[currentIndex] < peakWidth)
    peakEnd <- peakEnd + 1
while(templateEnd <  maxIndex &&
      masses[templateEnd] - masses[currentIndex] < templateWidth)
    templateEnd <- templateEnd + 1

## Now, while there is room to move the template end upwards ...
while(templateEnd < maxIndex) {
    
    ## Compute the normalizer and the central peaks in the
    ## current window
    normalizer[currentIndex] <-
        max(c(max(intensities[templateStart:(peakStart-1)]),
              max(intensities[(peakEnd+1):templateEnd])))
    peaks[currentIndex] <- max(intensities[peakStart:peakEnd])
    
    ## Move up the current index ...
    currentIndex <- currentIndex+1
    
    ## And like above, update the limits of the peak and templates
    while(peakStart < currentIndex &&
          masses[currentIndex] - masses[peakStart+1] > peakWidth)
        peakStart <- peakStart + 1
    while(templateStart <  currentIndex &&
          masses[currentIndex]-masses[templateStart+1] > templateWidth)
        templateStart <- templateStart + 1
    while(templateEnd <  maxIndex &&
          masses[templateEnd+1] > masses[templateEnd] &&
          masses[templateEnd] - masses[currentIndex] < templateWidth)
        templateEnd <- templateEnd + 1
    while(peakEnd <  templateEnd &&
          masses[peakEnd] - masses[currentIndex] < peakWidth)
        peakEnd<- peakEnd + 1
    
}
    
## Normalize and find the ones where the level is higher than the area
## around it by the threshold and where the peak is at the intensity
## of the point. the variable indices holds the indices of the
## identified peaks.
##
## NOTE: At this point, to make debugging faster, we limit the mass to 2000.
temp <- intensities / normalizer
indices <- which(temp > threshold & intensities == peaks & masses < 4000)

## now let's sort those indices by intensity and pick the top 8
sorted <- sort(intensities[indices],index.return=TRUE,decreasing=TRUE)
indices <- indices[sorted$ix]
if (length(indices) > 8) {
    indices <- indices[1:8]
}

## Now plot the mass spec. First just where the peaks are and then
## the entire one so that we can see where it found the peaks.
png(paste("./mic-output/MassSpec-",ID,".png",sep=""))
limits <- c(0,1.1*max(intensities))

if (length(indices) > 0) {
    mass = masses[indices[1]]
    intensity = intensities[indices[1]]
    titleString = ""
} else {
    titleString = "Mass Spec: NO PEAKS FOUND"
}

plot(masses,
     intensities,
     type="l",
     main=titleString,
     xlab="Masses (da)",
     ylab="Intensities",ylim=limits)

if (length(indices) > 0) {

    ## Put cirlces at the peaks and a solid one where we compute the MassSpec
    points(masses[indices],intensities[indices],type="p")
    points(mass,intensity,type="p",pch=19)

    ## Now, let's get the compositions by running the code. We had to
    ## put a link to the executable in a path that I could execute
    ## from. This is that path.
    command <-  paste("/usr/local/bin/computeParallelPeptideComposition ", ID, mass);
    print(paste("Excecute Command: ",command))
    system(command)

    ## Now read in the generated file and count the number of
    ## rows so we can annotate the display and put a composition in the title
    acidNames = c("G","A","S","P","V","T","C","L","N","D",
                  "Q","K","E","M","H","F","R","Y","W")
    if (length(indices) > 0) {
        fileName <- paste("Compositions-",ID,"-0.csv",sep='')
        if (file.exists(fileName)) {
            composition <- read.csv(fileName)
            print(paste("Found ",nrow(composition)," Compositions"))
            titleString = "Mass Spec With Composition: "
            if (nrow(composition) > 1) {
                titleString = paste(titleString," (1 of ",nrow(composition),")")
            }
            subString = ""
            for (acidIndex in 1:length(acidNames)) {
                subString = paste(subString,
                                  acidNames[acidIndex],
                                  composition[1,acidIndex],
                                  sep="");
            }
            title(main=titleString);
            mtext(subString,side=3,line=0)
        } else {
            title(main="Mass Spec: Inversion Error")
        }
    }
}

## This flushes the plotting to the file
dev.off()

