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
	
	<?php include("top_header_start.php"); include("top_header_logo.php"); include("top_header_menu.php");
	
	echo "<div class=\"container\"><div class=\"row\"><div class=\"col-sm-12\"><br/><br/><table>\n\n";
		//uploaded file parameters
	$_FILES["upload_file"]["name"];     // file's original filename
	$_FILES["upload_file"]["type"];     // file's mimetype
	$_FILES["upload_file"]["tmp_name"]; // temporary filename
	$_FILES["upload_file"]["size"];     // file's size
	$_FILES["upload_file"]["error"];    // error code
	
        $bacteria = $_POST["Bacteria_Name"];
        $assay = $_POST["Assay_Name"];
        $medium = $_POST["Medium_Name"];
        $peptide1 = $_POST["Peptide1_Name"];
        $peptide2 = $_POST["Peptide2_Name"];
        $antibiotic = $_POST["Antibiotic_Name"];
        $timepoint = $_POST["Time_Point"];
        $con = $_POST["Con_Name"];
	
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
	?>

	<?php
	$Rcommand = "Rscript R-scripts/mic.R $idcode $bacteria $assay $medium $peptide1 $peptide2 $antibiotic $con $timepoint";
	//$Rcommand = "Rscript R-scripts/mic.R $idcode $timepoint";
	exec($Rcommand); // execute R script	
	$micpng = "./mic-output/mic-plot-".$idcode.".png";
	echo "&nbsp;&nbsp;&nbsp&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<br/>";
	echo "<br/><br/>";
	echo "\n</table>";
	echo "<br/>";
	echo "<br/><br/>";
	echo '<img src="'.$micpng.'">';
	?>
<?php include("footer.php"); ?>
</body>
</html>
