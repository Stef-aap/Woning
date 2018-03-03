
<?php
$Indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;' ;

//list all relevant files
$files = glob ( '*.{htm,html,jpg,png,gif}', GLOB_BRACE);
foreach($files as $file){
   echo $Indent.'<a href="'.$file.'">'.$file.'</a><br />';
}

echo '<br>';
echo $Indent.$Indent.'<a href=..>'.'Pagina naar boven'.'</a><br />';
echo $Indent.$Indent.'<a href=/>'.'Top Pagina'.'</a><br />';

// list all subdirectories
$files = array_filter ( glob('*'), 'is_dir' );
foreach($files as $file){
   echo $Indent.'<a href="'.$file.'">'.$file.'</a><br />';
}
?>

