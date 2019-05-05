<?php session_start(); 
if ($_SESSION['is_auth'] != true) {
	$_SESSION['redirect'] = "2toolspage.php";
	header ("location: ../login.php"); exit();
}
?>
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
		include("top_header_start.php"); include("top_header_logo.php"); include("top_header_menu.php");
		?>
		<div class="container">
		<div class="row">
			<div class="col-sm-4" style="padding-left: 250px;">
              <h4 style="text-align: left;"><span style="font-weight: normal;"><a href="https://dev.amped.uri.edu/mic.php"><img src ="/images/mictool.png" alt="" style="width:200px;height:200px;" /> <br/><a href="https://dev.amped.uri.edu/mic.php">MIC/MLC Calculator</a></span><br/>
			</div>
			<div class="col-sm-6" style="padding-left: 200px;">
				<h4 style="text-align: left;"><span style="font-weight: normal;"><a href="https://dev.amped.uri.edu/mic.php"><img src ="/images/ms.png" alt="" style="width:200px;height:200px;" /> <br/><a href="https://dev.amped.uri.edu/mic.php">Mass Spectrometry Calculator</a></span><br/>
			</div>
		</div>
		</div>
	
	
		<?php include("footer.php"); ?>	
	</body>
</html>	