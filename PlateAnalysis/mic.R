require(ggplot2)
require(reshape2)
#require(gridExtra)
#library(dplyr)
#library(tidyr)
#library(tidyverse)
#library(data.table)
library("tools")

#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
args <- commandArgs(TRUE)
id <- args[1]
bacteria <- args[2]
assay <- args[3]
medium <- args[4]
testPep1 <- args[5]
testPep2 <- args[6]
antibiotic <- args[7]
antibiotic_concentration <- args[8]
timepoints <- args[9]

timepoints <- as.list(strsplit(timepoints, ",")[[1]])
timepoints <- as.numeric(timepoints)

# determine file type
file <- paste("./mic-data/",id,"-data.csv",sep="")

# used with csv in correct format
df2 = read.csv(file, header = FALSE)
numify <- function(x) as.numeric(as.character(x))
df2[] <- lapply(df2, numify)

# find number of plates
num_rows <- nrow(df2)
num_plates <- num_rows / 9

# split dataframe into seperate plates
plates <- split(df2,rep(1:num_plates, each=9))

# average blanks for each plate and add to list
blanks = list()
for (i in 1:(num_plates / 2)) {
  blanks1 = mean(plates[[i]][2:7, 1])
  blanks2 = mean(plates[[i+1]][2:7, 1])
  blanks [[i]] <- ((blanks1 + blanks2) / 2)
}

# data frames to hold average concentrations
averageGrowthPl1_1 <- data.frame(matrix(, nrow=num_plates/2, ncol=10)) # plates rows 1-3
averageGrowthPl1_2 <- data.frame(matrix(, nrow=num_plates/2, ncol=10)) # plates rows 4-6
controlRow1 <- data.frame(matrix(, nrow=num_plates/2, ncol=10))
lastRow1 <- data.frame(matrix(, nrow=num_plates/2, ncol=10))

# add concentrations to column names
initial_conc <- as.numeric(antibiotic_concentration)
concentrations = list()
for (i in 1:10) {
  concentrations[i] <- initial_conc
  initial_conc <- initial_conc / 2
}
colnames(averageGrowthPl1_1) <- c(concentrations)
colnames(averageGrowthPl1_2) <- c(concentrations)
colnames(controlRow1) <- c(concentrations)
colnames(lastRow1) <- c(concentrations)

# create lists to hold average concentrations (lists will then be used to create row in data frame)
averageGrowthListPl1_1 = list() # Set of plates including plate 1 rows 1-3
averageGrowthListPl1_2 = list() # Set of plates including plate 1 rows 4-6
controlRowList1 = list()
lastRowList1 = list()

# find averages for all columns including both plates (peptides, second to last row and last row)
for (i in seq(1,num_plates,2)) {
  for (j in 1:10) {
    # Peptide 1
    plate1Avg_1 <- mean(plates[[i]][2:4, j+1])
    plate2Avg_1 <- mean(plates[[i+1]][2:4, j+1])
    averageGrowthListPl1_1[[j]] <- ((plate1Avg_1 + plate2Avg_1) / 2)
    
    # Peptide 2
    plate1Avg_2 <- mean(plates[[i]][5:7, j+1])
    plate2Avg_2 <- mean(plates[[i+1]][5:7, j+1])
    averageGrowthListPl1_2[[j]] <- ((plate1Avg_2 + plate2Avg_2) / 2)
    
    # Control (second to last row)
    control1 <- plates[[i]][8, j+1]
    control2 <- plates[[i+1]][8, j+1]
    controlRowList1[[j]] <- ((control1 + control2) / 2)
    
    # ? (last row)
    lastRowVal1 <- plates[[i]][9, 2]
    lastRowVal2 <- plates[[i+1]][9, 2]
    lastRowList1[[j]] <- ((lastRowVal1 + lastRowVal2) / 2)
  }
  averageGrowthPl1_1[i, ] <- averageGrowthListPl1_1
  averageGrowthPl1_2[i, ] <- averageGrowthListPl1_2
  controlRow1[i, ] <- controlRowList1
  lastRow1[i, ] <- lastRowList1
}

# find average of bacteria (last column)
averageGrowthBacteria <- data.frame(matrix(, nrow=num_plates/2, ncol=1))
averageGrowthListBacteria = list() 
for (i in seq(1,num_plates,2)) {
  for (j in 1:1) {
    averageGrowthLast_1 = mean(plates[[i]][2:7, 12])
    averageGrowthLast_2 = mean(plates[[i+1]][2:7, 12])
    averageGrowthListBacteria[[j]] <- ((averageGrowthLast_1 + averageGrowthLast_2) / 2)
  }
  averageGrowthBacteria[i, ] <- averageGrowthListBacteria
}

# subtract blanks from average concentrations
for (i in 1:(num_plates / 2)) {
  averageGrowthPl1_1[i, ] <- averageGrowthPl1_1[i, ] - blanks[[i]]
  averageGrowthPl1_2[i, ] <- averageGrowthPl1_2[i, ] - blanks[[i]]
  controlRow1[i, ] <- controlRow1[i, ] - blanks[[i]]
  lastRow1[i, ] <- lastRow1[i, ] - blanks[[i]]
  averageGrowthBacteria[i, ] <- averageGrowthBacteria[i, ] - blanks[[i]]
}


# Remove NA's from rows
averageGrowthPl1_1 <- averageGrowthPl1_1[complete.cases(averageGrowthPl1_1), ]
averageGrowthPl1_2 <- averageGrowthPl1_2[complete.cases(averageGrowthPl1_2), ]
averageGrowthBacteria <- averageGrowthBacteria[complete.cases(averageGrowthBacteria), ]
controlRow1 <- controlRow1[complete.cases(controlRow1), ]
lastRow1 <- lastRow1[complete.cases(lastRow1), ]

# ADD TIMES (NEED TO CHANGE THIS TO GRAB INPUT FROM USER)
"
Add timepoints where test1 is
"
# remove following 3 lines and add timepoints where test1 is in next 4 lines
averageGrowthPl1_1$time <- c(timepoints)
averageGrowthPl1_2$time <- c(timepoints)
controlRow1$time <- c(timepoints)
lastRow1$time <- c(timepoints)

# find row that is first time over 1440 minutes
num_columns = ncol(averageGrowthPl1_1)
mic_row = 0
for (i in 1:num_rows) {
  if (averageGrowthPl1_1[i, "time"] > 1440) {
    mic_row = i
    i
    break
  }
}
png(paste("./mic-output/mic-plot-",id,".png",sep=""),width = 800,height=1000)

#### Plot peptide 1 MIC after 1440 minutes ####
growth_row1 = ((averageGrowthPl1_1[mic_row,] - averageGrowthPl1_1[1,]) / (averageGrowthBacteria[mic_row] - averageGrowthBacteria[1])) * 100
growth_row1 <- growth_row1[,1:(num_columns-1)]
print(paste("Type of growth",typeof(growth_row1)))
#setattr(growth_row1, "row.names", c("growth"))
row.names(growth_row1) <- "growth"
growth_row1 <- as.data.frame(t(growth_row1))
growth_row1$concentration <- rownames(growth_row1)
growth_row1$concentration <- factor(growth_row1$concentration, levels = growth_row1$concentration)
# reverse order of df so low concentrations are first
growth_row1 <- growth_row1[seq(dim(growth_row1)[1],1),]
growth_row1$concentration <- factor(growth_row1$concentration, levels = growth_row1$concentration)
plot1 <- ggplot(growth_row1, aes(x=concentration, y=growth)) + ggtitle(testPep1) + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_point()

#### Plot peptide 2 MIC after 1440 minutes ####
growth_row2 = ((averageGrowthPl1_2[mic_row,] - averageGrowthPl1_2[1,]) / (averageGrowthBacteria[mic_row] - averageGrowthBacteria[1])) * 100
growth_row2 <- growth_row2[,1:(num_columns-1)]
#setattr(growth_row2, "row.names", c("growth"))
row.names(growth_row2) <- "growth"
growth_row2 <- as.data.frame(t(growth_row2))
growth_row2$concentration <- rownames(growth_row2)
growth_row2$concentration <- factor(growth_row2$concentration, levels = growth_row2$concentration)
# reverse order of df so low concentrations are first
growth_row2 <- growth_row2[seq(dim(growth_row2)[1],1),]
growth_row2$concentration <- factor(growth_row2$concentration, levels = growth_row2$concentration)
plot2 <- ggplot(growth_row2, aes(x=concentration, y=growth)) + ggtitle(testPep2) + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_point() 

#### Plot control MIC after 1440 minutes ####
growth_row3 = ((controlRow1[mic_row,] - controlRow1[1,]) / (averageGrowthBacteria[mic_row] - averageGrowthBacteria[1])) * 100
growth_row3 <- growth_row3[,1:(num_columns-1)]
#setattr(growth_row3, "row.names", c("growth"))
rownames(growth_row3) <- "growth"
growth_row3 <- as.data.frame(t(growth_row3))
growth_row3$concentration <- rownames(growth_row3)
growth_row3$concentration <- factor(growth_row3$concentration, levels = growth_row3$concentration)
# reverse order of df so low concentrations are first
growth_row3 <- growth_row3[seq(dim(growth_row2)[1],1),]
growth_row3$concentration <- factor(growth_row3$concentration, levels = growth_row3$concentration)
plot3 <- ggplot(growth_row3, aes(x=concentration, y=growth)) + ggtitle(antibiotic) + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_point() 

micPeptide1 = "16"
micPeptide2 = "4"
mlcPeptide1 = "32"
mlcPeptide2 = "16"

# Plot peptide 1 mic, mlc and control
pep1_df <- data.frame(time = averageGrowthPl1_1[,'time'],
                 'Peptide 1 MIC' = averageGrowthPl1_1[,micPeptide1],
                 'Peptide 1 MLC' = averageGrowthPl1_1[,mlcPeptide1],
                 'Control Antibiotic' = controlRow1[,micPeptide1],
                 'Broth?' = lastRow1[,micPeptide1])
pep1_df <- melt(pep1_df ,  id.vars = 'time', variable.name = 'key')
pep1_df
plot4 <- ggplot(pep1_df, aes(time, value)) + geom_line(aes(colour = key)) + ggtitle(testPep1) +
  theme(plot.title = element_text(hjust = 0.5))

# Plot peptide 2 mic, mlc and control
pep2_df <- data.frame(time = averageGrowthPl1_2[,'time'],
                      'Peptide 2 MIC' = averageGrowthPl1_2[,micPeptide2],
                      'Peptide 2 MLC' = averageGrowthPl1_2[,mlcPeptide2],
                      'Control Antibiotic' = controlRow1[,micPeptide2],
                      'Broth?' = lastRow1[,micPeptide2])
pep2_df <- melt(pep2_df ,  id.vars = 'time', variable.name = 'key')
pep2_df
plot5 <- ggplot(pep2_df, aes(time, value)) + geom_line(aes(colour = key)) + ggtitle(testPep2) +
  theme(plot.title = element_text(hjust = 0.5))

#grid.arrange(plot1, plot2, plot3, plot4, plot5, ncol=3)
multiplot(plot1,plot2,plot3,plot4,plot5,cols=2)
dev.off()







