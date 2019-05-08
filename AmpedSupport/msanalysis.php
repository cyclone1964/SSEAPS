<?Php session_start();
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
    <script type="text/javascript">
	
	function isFilled(elm) {
		if (elm == "" ||
			elm == null) {
			return false;
		}
		else {
			return true;
		}
	}
	</script>
	</head>
	<body>

     <?php include("top_header_start.php"); include("top_header_logo.php"); include("top_header_menu.php"); 
         ?>
	<div class="container" >
	 <form action="msanalysis_upload_process.php" method="post" enctype="multipart/form-data" name="form1" id="form1" onSubmit="return valForm(this);"> 	 
		<div class="row">
            <div class="col-sm-6">
              &nbsp;
            </div>
            <div align="center";>
              <h3>Mass Spec Analysis</h3>
            </div>
        </div>
		<br>
           <div class="row">
            <div class="col-sm-4">
			<img src="images/msintensity.png" style="width:200px;height:140px;">
              <h5>Upload your .csv file in format shown in the template</br></h5>
            </div>
			<div class="col-sm-6">
				<div style="text-align: left; padding-left: 150px"><input style="height: 35px;" name="upload_file" type="file" accept= " .csv" id="upload_file" value=""/>
				</div>
			    <div style="text-align: left; padding-left: 150px"><input style="height: 35px; width: 150px; font-size: 14px;"
                  value="Submit" id="Submit" name="Submit" type="Submit">
			    </div>
		        <br>
		    </div>
		   </div>			
		</form>
    </div>
	<?php include("footer.php"); ?>
  </body>
</html>
