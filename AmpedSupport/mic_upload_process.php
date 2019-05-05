<?php session_start(); ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
  <head>
    <meta content="text/html; charset=windows-1252" http-equiv="content-type">
    <link href="css/bootstrap.min.css" rel="stylesheet">
    <link href="css/main.css" rel="stylesheet">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	<link rel="shortcut icon" href="/favicon.ico" type="image/x-icon">
	<link rel="icon" href="/favicon.ico" type="image/x-icon">
	</head>
	<body>
	<?php
   $bacteriaName	= $_POST["Bacteria_Name"];
   $assay			= $_POST["Assay_Name"];
   $medium			= $_POST["Medium_Name"];
   $peptide1		= $_POST["Peptide1_Name"];
   $peptide2		= $_POST["Peptide2_Name"];
   $antibiotic		= $_POST["Antibiotic_Name"];
   $con				= $_POST["Con_Name"];
   $timepoint		= $_POST["Time_Point"];
	?>

	<?php
	//function to open .r from php and write rscript everytime a file is uploaded with dyanmic/unique name.
	function writeFile($idcode) {
	global $fileName;
	$text1 = "#m_m\ndataFile <- read.delim(\"mic-data/".$fileName. "\", header=TRUE, sep=\",\", colClasses = \"character\")\n";
	$constantTXT = <<<END
            library(dplyr)
            library(tidyr)

            dataFile = read.delim(file.choose(), header= FALSE);

            data = dataFile[,c(3,4,5,6,7,8,9,10,11,12,13,14)];
            df2 <- data[complete.cases(data), ]
            numify <- function(x) as.numeric(as.character(x))
            df2[] <- lapply(df2, numify)
            d
            # find number of plates
            num_rows <- nrow(df2)
            num_plates <- num_rows / 9

            # split dataframe into seperate plates
            plates <- split(df2,rep(1:num_plates, each=9))

            blanks = list()
            for (i in 1:num_plates) {
              blanks [[i]] <- mean(plates[[i]][2:7, 1])
            }

            # data frames to hold average concentrations
            averageGrowthPl1_1 <- data.frame(matrix(, nrow=num_plates/2, ncol=10)) # odd plates rows 1-3
            averageGrowthPl1_2 <- data.frame(matrix(, nrow=num_plates/2, ncol=10)) # odd plates rows 4-6
            averageGrowthPl2_1 <- data.frame(matrix(, nrow=num_plates/2, ncol=10)) # even plates rows 1-3
            averageGrowthPl2_2 <- data.frame(matrix(, nrow=num_plates/2, ncol=10)) # even plates rows 4-6

            # create lists to hold average concentrations (lists will then be used to create row in data frame)
            averageGrowthListPl1_1 = list() # Set of plates including plate 1 rows 1-3
            averageGrowthListPl1_2 = list() # Set of plates including plate 1 rows 4-6
            averageGrowthListPl2_1 = list() # Set of plates including plate 2 rows 1-3
            averageGrowthListPl2_2 = list() # Set of plates including plate 2 rows 4-6

            # find averages for plates including first and add them to list then add list to data frame
            for (i in seq(1,num_plates,2)) {
              for (j in 1:10) {
                averageGrowthListPl1_1[[j]] <- mean(plates[[i]][2:4,j+1])
                averageGrowthListPl1_2[[j]] <- mean(plates[[i]][5:7,j+1])
              }
              averageGrowthPl1_1[i, ] <- averageGrowthListPl1_1
              averageGrowthPl1_2[i, ] <- averageGrowthListPl1_2
            }

            # find averages for plates including second and add them to list then add list to data frame
            for (i in seq(2, num_plates, 2)) {
              for (j in 1:10) {
                averageGrowthListPl2_1[[j]] <- mean(plates[[i]][2:4,j+1])
                averageGrowthListPl2_2[[j]] <- mean(plates[[i]][5:7,j+1])
              }
              averageGrowthPl2_1[i, ] <- averageGrowthListPl2_1
              averageGrowthPl2_2[i, ] <- averageGrowthListPl2_2
            }

            averageGrowthLastCol <- data.frame(matrix(, nrow=num_plates, ncol=1))
            averageGrowthListLastCol = list() # Set of plates including plate 1 rows 2-4
            # find averages for plates including second and add them to list then add list to data frame
            for (i in 1:num_plates) {
              for (j in 1:1) {
                averageGrowthListLastCol[[j]] <- mean(plates[[i]][2:7, 12])
              }
              averageGrowthLastCol[i, ] <- averageGrowthListLastCol
            }

            # ADD ALL OF THESE TO LIST TO MAKE CODE CLEANER
            # subtract blanks from average concentrations
            for (i in 1:num_plates) {
              averageGrowthPl1_1[i, ] <- averageGrowthPl1_1[i, ] - blanks[[i]]
              averageGrowthPl1_2[i, ] <- averageGrowthPl1_2[i, ] - blanks[[i]]
              averageGrowthPl2_1[i, ] <- averageGrowthPl2_1[i, ] - blanks[[i]]
              averageGrowthPl2_2[i, ] <- averageGrowthPl2_2[i, ] - blanks[[i]]
              averageGrowthLastCol[i, ] <- averageGrowthLastCol[i, ] - blanks[[i]]
            }

            # Remove NA's from rows
            averageGrowthPl1_1 <- averageGrowthPl1_1[complete.cases(averageGrowthPl1_1), ]
            averageGrowthPl1_2 <- averageGrowthPl1_2[complete.cases(averageGrowthPl1_2), ]
            averageGrowthPl2_1 <- averageGrowthPl2_1[complete.cases(averageGrowthPl2_1), ]
            averageGrowthPl2_2 <- averageGrowthPl2_2[complete.cases(averageGrowthPl2_2), ]
            averageGrowthLastCol <- averageGrowthLastCol[complete.cases(averageGrowthLastCol), ]


            is.num <- sapply(averageGrowthPl1_1, is.numeric)
            averageGrowthPl1_1[is.num] <- lapply(averageGrowthPl1_1[is.num], round, 5)
            is.num <- sapply(averageGrowthPl1_2, is.numeric)
            averageGrowthPl1_2[is.num] <- lapply(averageGrowthPl1_2[is.num], round, 5)
            is.num <- sapply(averageGrowthPl2_1, is.numeric)
            averageGrowthPl2_1[is.num] <- lapply(averageGrowthPl2_1[is.num], round, 5)
            is.num <- sapply(averageGrowthPl2_2, is.numeric)
            averageGrowthPl2_2[is.num] <- lapply(averageGrowthPl2_2[is.num], round, 5)
            is.num <- sapply(averageGrowthLastCol, is.numeric)
            averageGrowthLastCol[is.num] <- lapply(averageGrowthLastCol[is.num], round, 5)

            # need to remember why I did this
            averageGrowthPl1_1[nrow(averageGrowthPl1_1)+1,] <- NA
            averageGrowthPl1_2[nrow(averageGrowthPl1_2)+1,] <- NA
            averageGrowthPl2_1[nrow(averageGrowthPl2_1)+1,] <- NA
            averageGrowthPl2_2[nrow(averageGrowthPl2_2)+1,] <- NA

            # add to table to then write to csv?
            write.table(averageGrowthPl1_1, "filename.csv", col.names=FALSE, sep=",")
            write.table(averageGrowthPl1_2, "filename.csv", col.names=FALSE, sep=",", append=TRUE)
            write.table(averageGrowthPl2_1, "filename.csv", col.names=FALSE, sep=",", append=TRUE)
            write.table(averageGrowthPl2_2, "filename.csv", col.names=FALSE, sep=",", append=TRUE)
            write.table(averageGrowthLastCol, "filename.csv", col.names=FALSE, sep=",", append=TRUE)

            # write to csv
            write.csv(rbind(averageGrowthPl1_1, averageGrowthPl1_2,
                            averageGrowthPl2_1, averageGrowthPl2_2, averageGrowthLastCol), "filename.csv")


END;
	$text2 = "png(file = \"mic-output/".$fileName. "-pep1.png\")\n
	matplot(Pep1\$Time.point,Pep1[,c(1:12)],type = c(\"b\"),lty = 1,lwd =2,  xlab = \"Time Points\",ylab = \"OD\",main = \"MIC_MLC Peptide1\",pch=c(15,16,17),col = 1:12)\n
	legend(\"topleft\", legend = colnames(Pep1), col=1:12, box.lwd =par(\"lwd\") ,pch=c(15,16,17)) # optional legend\n
	dev.off()\n";
	$text3 = "png(file = \"mic-output/".$fileName. "-pep2.png\")\n
	matplot(Pep2\$Time.point,Pep2[,c(1:12)],type = c(\"b\"),lty = 1,lwd =2,  xlab = \"Time Points\",ylab = \"OD\",main = \"MIC_MLC Peptide2\",pch=c(15,16,17),col = 1:12)\n
	legend(\"topleft\", legend = colnames(Pep2), col=1:12, box.lwd =par(\"lwd\") ,pch=c(15,16,17)) # optional legend\n
	dev.off()\n";
	$text4 = "TotalPep <- rbind(Pep1,Pep2)\n
	write.csv(TotalPep,\"mic-output/".$fileName. "-result.csv\")";

    $myfile= fopen(getcwd()."/R-scripts/" . $idcode . "-mic.R", 'w') or die("Could not create the file.");
    fwrite($myfile, $text1);
    fwrite($myfile, $constantTXT);
	fwrite($myfile, $text2);
	fwrite($myfile, $text3);
	fwrite($myfile, $text4);
    fclose($myfile);
	}

	include("top_header_start.php"); include("top_header_logo.php"); include("top_header_menu.php");
	echo "<div class=\"container\"><div class=\"row\"><div class=\"col-sm-12\"><br/><br/><table>\n\n";
		//uploaded file parameters
	$_FILES["upload_file"]["name"];     // file's original filename
	$_FILES["upload_file"]["type"];     // file's mimetype
	$_FILES["upload_file"]["tmp_name"]; // temporary filename
	$_FILES["upload_file"]["size"];     // file's size
	$_FILES["upload_file"]["error"];    // error code


	if ($_FILES["upload_file"]["error"] > 0) {
	$thisMsg = "<p>Error transfer";
	}
	$idcode= time();
	$fileName = $idcode. "-" . "data.csv";
	$finalLoc = getcwd()."/mic-data/".$fileName;
	$URLloc = "/mic-data/".$fileName;
	$result = move_uploaded_file($_FILES["upload_file"]["tmp_name"],$finalLoc);

	if (!$result) {
	$thisMsg = "<p>Unfortunately the file could not be uploaded: $fileName.</p>";
	} else {
	$thisMsg = "<p>Download excel file: <a href=\"$URLloc\" target=\"_blank\">$fileName</a>.</p>\n";
	}

	writeFile($idcode); // create R script
	$Rcommand = "Rscript R-scripts/" . $idcode . "-mic.R";
	exec($Rcommand); // execute R script

	$f = fopen("mic-output/" .$fileName. "-result.csv", "r");
	while (($line = fgetcsv($f)) !== false) {
        echo "<tr>"."&nbsp;&nbsp;&nbsp";
        foreach ($line as $cell) {
                echo "<td>"."&nbsp;&nbsp;&nbsp" . htmlspecialchars($cell) . "</td>";
        }
        echo "</tr>\n";
	}
	fclose($f);
	$pep1png = "mic-output/".$fileName."-pep1.png";
	$pep2png = "mic-output/".$fileName."-pep2.png";
	echo "&nbsp;&nbsp;&nbsp&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<br/>";
	echo "<br/><br/>";
	echo "\n</table>";
	echo "<br/>";
	echo $thisMsg;
	echo "<br/><br/>";
	echo '<img src="'.$pep1png.'">';
	echo '<img src="'.$pep2png.'">';
	echo "<br/>";
	echo '<a href="'.$pep1png.'" download>Download-PEP1</a><br/>';
	echo '<a href="'.$pep2png.'" download>Download-PEP2</a>';
	##echo user input
	echo "<p>Bacteria Name :$bacteriaName"."&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"."Peptide2 Name : $peptide2</p>";
	echo "<p>Assay Name :$assay"."&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"."Antibiotic Name :$antibiotic</p>";
	echo "<p>Medium :$medium"."&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"."Antibiotic concentration :$con</p>";
	echo "<p>Peptide1 Name :$peptide1"."&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"."Timepoint :$timepoint</p>";
?>
<?php include("footer.php"); ?>
</body>
</html>
