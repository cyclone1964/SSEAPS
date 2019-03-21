# analyzeMassSpec - analyze a Mass Spectrometer Output

analyzeMassSpec <- function(fileName,
                            templateWidth=16,
                            peakWidth = 6,
                            threshold=1.5,
                            maxPeptideLength=16) {
  
    ## First, load the file, which is presumed to be a 2 column CSV
    massSpec <- read.csv(fileName)

    ## Convert that to a matrix because I hate typing $ all the time
    ## and I can't rely upon the column header being consistent
    massSpecMatrix <- data.matrix(massSpec)

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
    limits <- mat.or.vec(numPoints,1)

    ## extract the masses and intensities individually
    masses <- massSpecMatrix[,1]
    intensities <- massSpecMatrix[,2]

    # WE make a growing list of the indices of the peaks
    peaks = c()

    ## Now, we have to go through data and compute those things
    for (index in 1:matrixSize[1]) {
        window <- intensities[which((masses > masses[index] - templateWidth/2 & 
                                     masses < masses[index] - peakWidth/2) |
                                    (masses > masses[index] + peakWidth/2 &
                                     masses < masses[index] + templateWidth/2))]
        center <- intensities[which(masses >= masses[index] - peakWidth/2 & 
                                    masses <= masses[index] + peakWidth/2)]
        
        normalizer[index] <- max(window)
        limits[index] <- max(center)
    }

    ## Normalize and find the ones where the level is higher than the
    ## area around it by the threshold
    temp = intensities / normalizer
    peaks = which(temp > threshold & intensities == limits)

    ## Now plot the mass spec. First just where the peaks are and then
    ## the entire one so that we can see where it found the peaks.
    plot(masses[peaks],
         intensities[peaks],
         main=fileName,
         xlab="Masses (da)",
         ylab="Intensities")
    lines(masses,intensities)
  
    ## Now, let's get the compositions by running the code. We had to
    ## put a link to the executable in a path that I could execute
    ## from. This is that path.
    path = "/Users/Matt"
    for (index in 1:length(peaks)){
        command <- paste(path,
                         "/computePeptideComposition ",
                         masses[peaks[index]],
                         " ",
                         maxPeptideLength,
                         " Composition.csv",
                         sep='')
        print(paste("Excecute Command: ",command))
        system(command)
        composition <- read.csv("Composition.csv")
        print(paste("Found ",nrow(composition)," Compositions"))
        text(x=masses[peaks[index]],
             y=intensities[peaks[index]],
             paste(nrow(composition)),
             adj=c(1,0.5))
    }
}

