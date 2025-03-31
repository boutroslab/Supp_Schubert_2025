#!/usr/bin/perl -s -w
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Capture::Tiny ':all';
use File::Find;

my $barcode=$ARGV[0];
my $output_folder=$ARGV[1];
my $tilesize=$ARGV[2];
my $format=$ARGV[3];

my %Q=();
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){		
		foreach my $number (1..24){
			foreach my $field (1..4){
				my @temp=($letter,$number,$field);
				$Q{$letter.$number."_".$field} = \@temp;
			}
		}
	}
	
while(%Q){
    foreach my $key (sort keys(%Q)){
            my $letter=@{$Q{$key}}[0];
            my $number=@{$Q{$key}}[1];
            my $field=@{$Q{$key}}[2];
            if(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab"){
            #convert the R result
                    delete $Q{$key};
                    $base64_string = capture_stdout {system('convert '.$output_folder."/".$barcode."_".$letter.$number."_".$field.'_segmented.png -thumbnail 50x50 png:- | openssl enc -base64  ')}; #-out '.$output_folder.'/'.$letter.$number.'.b64'
                    #system("chmod -R a+r ".$output_folder."/".$barcode."_".$letter.$number."_".$field.'_segmented.tif ');
            #read in the tab-delim line for this image
                    my $count=0;
                    my @names=();
                    my @result=();
                    open $result_tab , $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab";
                            while(<$result_tab>){
                                    chomp $_;
                                    if($count<=0){
                                            @names=split("\t",$_);
                                            $count++;
                                    }elsif($count>=1){
                                            @result=split("\t",$_);
                                    }		
                            }
                    close $result_tab;								
                 if(!(-e $output_folder."/plate_image.svg")){
                            make_new_plate_svg($output_folder."/plate_image",$format,$tilesize,"plate_image");	#start new plate svg file
                            add_image_to_svg($letter.$number."_".$field,$base64_string,$output_folder,$tilesize,$format,$barcode);
                            close_svg($output_folder."/plate_image.svg");
                    }else{
                            open_svg($output_folder."/plate_image.svg");
                            add_image_to_svg($letter.$number."_".$field,$base64_string,$output_folder,$tilesize,$format,$barcode);
                            close_svg($output_folder."/plate_image.svg");
                    }
                    my $element=0;
                    foreach my $key (@names){
                            my $value=$result[$element];							
                            $element++;
                            if(!(-e $output_folder."/".$key.".svg")){
                                    make_new_plate_svg($output_folder."/".$key,$format,$tilesize,$key);	#start new plate svg file
                                    add_rect_to_svg($letter.$number."_".$field,$value,$output_folder,$key,$tilesize,$format,$barcode);
                                    adjust_colors_in_svg($output_folder,$key);
                                    close_svg($output_folder."/".$key.".svg");
                            }else{
                                    open_svg($output_folder."/".$key.".svg");
                                    add_rect_to_svg($letter.$number."_".$field,$value,$output_folder,$key,$tilesize,$format,$barcode);
                                    adjust_colors_in_svg($output_folder,$key);
                                    close_svg($output_folder."/".$key.".svg");
                            }
                    }
                                       
            }
    }
    combine_to_html($output_folder);    
    sleep 1;
}



#########################################################################################
#name:      remove all image files in folder
#function:  renames any tif file in a given folder in a way that each non word character 
#			or white space gets removed
#input:     (path/to/folder string) 
#output:    renamed files if open was possible
#########################################################################################
sub rename_tifs_in_folder {
	my $indir=$_[0];
	opendir $dir, $indir or die print "folder ".$indir." not found";
	foreach $file (sort(readdir($dir))){
		if($file=~m/(.+)\.tif/ && -e $indir.'/'.$file && !(-z $indir.'/'.$file)){
			$name=$newname=$1;
			$name=~s/\s/\\ /g;
			$name=~s/([\(\)])/\\$1/g;
			$newname=~s/\W//g;		
			if($name ne $newname){
				system('mv '.$indir.'/'.$name.'.tif '.$indir.'/'.$newname.'.tif');
			}
		}
	}
	closedir $dir;
	return(1);      
}
sub wanted {
	if($_=~m/(.+)\.tif/){
		my $newname=$1;
		$newname=~s/\W//g;
		$newname=~s/^_//g;
		$newname=~s/fld/_/g;
		$newname=~s/wv/_/g;
		$newname=~s/DAPIDAPI/DAPI/g;
		$newname=~s/Cy3Cy3/Cy3/g;
		$newname=~s/Cy5Cy5/Cy5/g;
		$newname=~s/FITCFITC/FITC/g;
		$newname=~s/0(\d)/$1/g;
		$File::Find::dir=~s/.+\/(\S+)_.+$/$1/g;
		rename($_,$File::Find::dir."_$newname.tif");
		
	}
}
#########################################################################################
#name:      make a new svg and write svg header
#function:  writes an browser readable svg header in a given file
#input:     (path/to/filename format tilesize [string] [6|12|24|96|384] [integer] ) 
#output:    a new svg file with header and platedesign
#########################################################################################
sub make_new_plate_svg {
	my $filename=$_[0];		#complete path to file including filename prefix "*.svg"
	my $format=$_[1];		#may be 6, 12, 24, 96, 384 well plate
	my $tilesize=$_[2];		#may be whatever
	my $id=$_[3];
	open my $output ,'>'.$filename.'.svg' or die $!;
	print $output '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';
	print $output '<svg version="1.1" onload="init(evt)" baseProfile="full" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:ev="http://www.w3.org/2001/xml-events">
	<script type="text/ecmascript">
		<![CDATA[
			function init(evt)
			{
				if ( window.svgDocument == null )
				{
					svgDocument = evt.target.ownerDocument;
				}
				tooltip = svgDocument.getElementById(\'tooltip\');
				tooltip_bg = svgDocument.getElementById(\'tooltip_bg\');
			}
			function ShowTooltip(evt, mouseovertext)
			{
				tooltip.setAttributeNS(null,"x",evt.clientX+11);
				tooltip.setAttributeNS(null,"y",evt.clientY+35);
 				tooltip.firstChild.data = mouseovertext;
 				
 				tooltip_bg.setAttributeNS(null,"x",evt.clientX+8);
				tooltip_bg.setAttributeNS(null,"y",evt.clientY+16);
				length = tooltip.getComputedTextLength();
				tooltip_bg.setAttributeNS(null,"width",length+8);			
				
 				tooltip_bg.setAttributeNS(null,"visibility","visible");
 				tooltip.setAttributeNS(null,"visibility","visible");
			}
			function HideTooltip()
			{
				tooltip.setAttributeNS(null,"visibility","hidden");
				tooltip_bg.setAttributeNS(null,"visibility","hidden");
			}
			function putToTop(evt, whom) {
				//get node reference
				var element = svgDocument.getElementById(\'whom\');
				//appendChild after the last child
				element.parentNode.appendChild(element);
			}
		]]>
	</script>
	<style>
		.tooltip{
			font-size: 20px;
		}
		.tooltip_bg{
			fill: white;
			stroke: black;
			stroke-width: 1;
			opacity: 0.85;
		}
	</style>
	'."\n";
	my $x=int($tilesize/2);
	my $y=$x+$tilesize;
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){
		print $output	'<text x="'.$x.'" y="'.$y.'" font-family="Verdana" font-size="12" fill="black" >'.$letter.'</text>'."\n";
		$y+=$tilesize;		
	}
	$y=int($tilesize/2);
	$x=$y+$tilesize;		
	foreach my $number (1..24){
		print $output	'<text x="'.$x.'" y="'.$y.'" font-family="Verdana" font-size="12" fill="black" >'.$number.'</text>'."\n";
		$x+=$tilesize;
	}
	print $output '<g id="'.$id.'">'."\n";
	close $output;
	return(1);      
}
#########################################################################################
#name:      add image to svg
#function:  adds an base64 encoded image to the prepared svg
#input:     (path/to/filename base64image /path/to/folder [string] [string] [string]) 
#output:    svg
#########################################################################################
sub add_image_to_svg {
	my $filename=$_[0];		#complete path to file including filename 
	my $b64=$_[1];
	my $output_folder=$_[2];
	my $tilesize=$_[3];
	my $format=$_[4];
	my $barcode=$_[5];
	my %coords=();
	my $x=int($tilesize/2);
	my $y=int($tilesize/2);
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){
		$y+=int($tilesize/2);
		$x=int($tilesize/2);			
		foreach my $number (1..24){
			foreach my $field (1..2){
				$x+=int($tilesize/2);
				my @temp=($x,$y);
				$coords{$letter.$number."_".$field}=\@temp;
			}
		}
		$y+=int($tilesize/2);
		$x=int($tilesize/2);
		foreach my $number (1..24){
			foreach my $field (3..4){
				$x+=int($tilesize/2);
				my @temp=($x,$y);
				$coords{$letter.$number."_".$field}=\@temp;
			}
		}
	}	
	open ($output, ">>".$output_folder."/plate_image.svg") or die $!;
		my $file=$barcode."_".$filename."_segmented.png";
		$x=@{$coords{$filename}}[0];
		$y=@{$coords{$filename}}[1];
		print $output '<image x="'.$x.'" y="'.$y.'" id="'.$filename.'" onmousemove="ShowTooltip(evt,\''.$filename.'\')" onmouseout="HideTooltip()" onclick="window.open(\''.$file.'\',\'_blank\');" onmouseover="this.style.cursor=\'pointer\'" width="'.(int($tilesize/2)-4).'px" height="'.(int($tilesize/2)-4).'px" xlink:href="data:image/png;base64,'.$b64.'"></image>'."\n";
	close $output;
	return(1);      
}
#########################################################################################
#name:      add rect to svg
#function:  adds an base64 encoded image to the prepared svg
#input:     (path/to/filename digital_value /path/to/folder [string] [double] [string]) 
#output:    svg
#########################################################################################
sub add_rect_to_svg {
	my $filename=$_[0];		#complete path to file including filename 
	my $value=$_[1];
	my $output_folder=$_[2];
	my $name=$_[3];
	my $tilesize=$_[4];
	my $format=$_[5];
	my $barcode=$_[6];
	my %coords=();
	my $x=int($tilesize/2);
	my $y=int($tilesize/2);
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){
		$y+=int($tilesize/2);
		$x=int($tilesize/2);			
		foreach my $number (1..24){
			foreach my $field (1..2){
				$x+=int($tilesize/2);
				my @temp=($x,$y);
				$coords{$letter.$number."_".$field}=\@temp;
			}
		}
		$y+=int($tilesize/2);
		$x=int($tilesize/2);
		foreach my $number (1..24){
			foreach my $field (3..4){
				$x+=int($tilesize/2);
				my @temp=($x,$y);
				$coords{$letter.$number."_".$field}=\@temp;
			}
		}
	}
	open ($output, ">>".$output_folder."/".$name.".svg") or die $!;
		my $file=$barcode."_".$filename."_segmented.png";
		$x=@{$coords{$filename}}[0];
		$y=@{$coords{$filename}}[1];
		print $output '<rect x="'.$x.'" y="'.$y.'" id="'.$filename.'" width="'.(int($tilesize/2)-4).'px" height="'.(int($tilesize/2)-4).'px" style="fill:rgb(0,0,0);" onmousemove="ShowTooltip(evt,\''.$value.'\')" onmouseout="HideTooltip()" onclick="window.open(\''.$file.'\',\'_blank\');" onmouseover="this.style.cursor=\'pointer\'" />'."\n";
	close $output;
	return(1);      
}
#########################################################################################
#name:      close svg
#function:  closes an browser readable svg in a given file if closed already than unclose
#input:     (path/to/filename string) 
#output:    a closed or open svg
#########################################################################################
sub open_svg {
	open ($FH, "+< $_[0]") or die "can't update $file: $!";
	while ( <$FH> ) {
		$addr = tell($FH) unless eof($FH);
		if(eof($FH)){
			truncate($FH, $addr) or die "can't truncate $file: $!";
			return (1);
		} 
	}
	close $FH;
	return(1);      
}
#########################################################################################
#name:      close svg
#function:  closes an browser readable svg in a given file if closed already than unclose
#input:     (path/to/filename string) 
#output:    a closed or open svg
#########################################################################################
sub close_svg {
	open my $output ,'>>'.$_[0] or die $!;
		print $output '</g><rect class="tooltip_bg" id="tooltip_bg" x="0" y="0" rx="4" ry="4" width="30" height="30" visibility="hidden"/><text class="tooltip" id="tooltip" x="0" y="0" visibility="hidden">Tooltip</text></svg>'."\n";
	close $output;
	return(1);      
}
#########################################################################################
#name:      create total svg 
#function:  collects all other svgs as grouped elements and combines them into one single svg
#input:     (/path/to/folder [string]) 
#output:    svg
#########################################################################################
sub combine_to_svg {
	my $output_folder=$_[0];	
	open ($output, ">".$output_folder."/total.svg") or die $!;
	print $output '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';
	print $output '<svg version="1.1" onload="init(evt)" baseProfile="full" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:ev="http://www.w3.org/2001/xml-events">
	<script type="text/ecmascript">
		<![CDATA[
			function init(evt)
			{
				if ( window.svgDocument == null )
				{
					svgDocument = evt.target.ownerDocument;
				}
				tooltip = document.getElementById(\'tooltip\');
				tooltip_bg = document.getElementById(\'tooltip_bg\');
			}
			function ShowTooltip(evt, mouseovertext)
			{
				tooltip.setAttributeNS(null,"x",evt.clientX+11);
				tooltip.setAttributeNS(null,"y",evt.clientY+35);
 				tooltip.firstChild.data = mouseovertext;
 				
 				tooltip_bg.setAttributeNS(null,"x",evt.clientX+8);
				tooltip_bg.setAttributeNS(null,"y",evt.clientY+16);
				length = tooltip.getComputedTextLength();
				tooltip_bg.setAttributeNS(null,"width",length+8);			
				
 				tooltip_bg.setAttributeNS(null,"visibility","visible");
 				tooltip.setAttributeNS(null,"visibility","visible");
			}
			function HideTooltip()
			{
				tooltip.setAttributeNS(null,"visibility","hidden");
				tooltip_bg.setAttributeNS(null,"visibility","hidden");
			}
			function putToTop(evt, whom) {
				//get node reference
				var element = svgDocument.getElementById(whom);
				//appendChild after the last child
				element.parentNode.appendChild(element);
				element.parentNode.appendChild(tooltip_bg);
				element.parentNode.appendChild(tooltip);
								
				others = svgDocument.getElementsByClassName("clickable");	
					for(var i=0 ; i<others.length; i++){
						others[i].setAttributeNS(null,"fill","black");	
					}
				var that = evt.target;
				that.setAttributeNS(null,"fill","blue");			
			}
		]]>
	</script>
	<style>
		.tooltip{
			font-size: 20px;
		}
		.tooltip_bg{
			fill: white;
			stroke: black;
			stroke-width: 1;
			opacity: 0.85;
		}
		.clickable { cursor: pointer;  }
		.clickable:hover { color: blue;  }
		.clickable:visited {
			color: red;
			font-family: Verdana;
			text-decoration:none;
		}
	</style>'."\n";
	opendir $outdir, $output_folder;
	my $y=25;
	my $firstfile=0;
	foreach my $file (readdir($outdir)){		
		my $bool=0;
		if($file=~m/(.*)\.svg/ && !($file=~m/total.svg/ )){
			my $name=$1;
			print $output '<text x="1380" y="'.$y.'"  class="clickable" font-family="Verdana" font-size="12" fill="black" onclick="putToTop(evt, \''.$name.'\')" >'.$name.'</text>'."\n";
			$y+=14;
			open $FH, "<".$output_folder."/".$file;
				while (<$FH>){
					my $line =$_;
					if($line=~m/\<g id.*/){
						print $output $line;
						$bool=1;
					}elsif($bool==1){
						if($line=~m/^(.*)<\/svg>/){
							print $output $1;
						}else{
							if($firstfile==0 && $line=~m/\<text .*/){
								print $output $line;
							}elsif(!($line=~m/\<text .*/)){
								print $output $line;
							}
						}
					}	
				}
			close $FH;
			$firstfile++;
		}
	}
	closedir $outdir;
		print $output '<rect class="tooltip_bg" id="tooltip_bg" x="0" y="0" rx="4" ry="4" width="30" height="30" visibility="hidden"/><text class="tooltip" id="tooltip" x="0" y="0" visibility="hidden">Tooltip</text></svg>'."\n";
	close $output;
	return(1);      
}
#########################################################################################
#name:      create total svg 
#function:  collects all other svgs as grouped elements and combines them into one single svg
#input:     (/path/to/folder [string]) 
#output:    svg
#########################################################################################
sub combine_to_html {
	my $output_folder=$_[0];	
	open ($output, ">".$output_folder."/total.html") or die $!;
	print $output '<!DOCTYPE html>
			<html lang="enc" onload="toggle_well_field(\'stuff\')">
			<head>
			<meta http-equiv="Content-Type" content="text/html;charset=utf-8" autocomplete="off">
			<meta name="Author" content="Florian Heigwer"><title>Image Surveillance</title>
			<script type="text/javascript">
			    
			    function init(evt)
			    {
				    if ( window.svgDocument == null )
				    {
					    svgDocument = evt.target.ownerDocument;
				    }
				    tooltip = svgDocument.getElementById(\'tooltip\');
				    tooltip_bg = svgDocument.getElementById(\'tooltip_bg\');
			    }
			    function ShowTooltip(evt, mouseovertext)
			    {
				    tooltip.setAttributeNS(null,"x",evt.clientX+11);
				    tooltip.setAttributeNS(null,"y",evt.clientY+35);
				    tooltip.firstChild.data = mouseovertext;
				    
				    tooltip_bg.setAttributeNS(null,"x",evt.clientX+8);
				    tooltip_bg.setAttributeNS(null,"y",evt.clientY+16);
				    length = tooltip.getComputedTextLength();
				    tooltip_bg.setAttributeNS(null,"width",length+8);			
				    
				    tooltip_bg.setAttributeNS(null,"visibility","visible");
				    tooltip.setAttributeNS(null,"visibility","visible");
			    }
			    function HideTooltip()
			    {
				    tooltip.setAttributeNS(null,"visibility","hidden");
				    tooltip_bg.setAttributeNS(null,"visibility","hidden");
			    }
			    function putToTop(evt, whom) {
				    //get node reference
				    var element = svgDocument.getElementById(whom);
				    //appendChild after the last child
				    element.parentNode.appendChild(element);
				    element.parentNode.appendChild(tooltip_bg);
				    element.parentNode.appendChild(tooltip);
								    
				    others = svgDocument.getElementsByClassName("clickable");	
					    for(var i=0 ; i<others.length; i++){
						    others[i].setAttributeNS(null,"fill","black");	
					    }
				    var that = evt.target;
				    that.setAttributeNS(null,"fill","blue");			
			    }
			    function change_plot(feature){
					var el = document.getElementById(\'plate_container\');
					el.innerHTML = \'<object id="plate_plot" type="image/svg+xml" width="100%" height="1000" data="\'+feature+\'.svg"></object>\'
			    }
			    Element.prototype.remove = function() {
				this.parentElement.removeChild(this);
			    }
			    NodeList.prototype.remove = HTMLCollection.prototype.remove = function() {
				for(var i = 0, len = this.length; i < len; i++) {
				    if(this[i] && this[i].parentElement) {
					this[i].parentElement.removeChild(this[i]);
				    }
				}
			    }
			    function adjust_colors(){
				var a= document.getElementById(\'plate_plot\')
				var color_vector=[
				    ["#67001F","#690423","#6C0927","#6F0E2C","#721330","#751834","#781D39","#7B223D","#7E2741","#812C46","#84314A","#87364F","#893B53","#8C4057","#8F455C","#924A60","#954F64","#985469","#9B596D","#9E5E71","#A16376","#A4687A","#A76D7F","#A97283","#AC7787","#AF7C8C","#B28190","#B58694","#B88B99","#BB909D","#BE95A1","#C19AA6","#C49FAA","#C7A4AF","#C9A9B3","#CCAEB7","#CFB3BC","#D2B8C0","#D5BDC4","#D8C2C9","#DBC7CD","#DECCD1","#E1D1D6","#E4D6DA","#E7DBDF","#E9E0E3","#ECE5E7","#EFEAEC","#F2EFF0","#F5F4F4","#F4F4F5","#EFF0F2","#EAECEF","#E5E8EC","#E0E4E9","#DCE0E6","#D7DCE3","#D2D8E0","#CDD4DD","#C8D0DA","#C3CCD7","#BEC8D4","#B9C4D1","#B4C0CE","#B0BCCB","#ABB8C8","#A6B4C4","#A1B0C1","#9CACBE","#97A8BB","#92A4B8","#8DA0B5","#899CB2","#8498AF","#7F94AC","#7A90A9","#758CA6","#7088A3","#6B84A0","#66809D","#617C9A","#5C7897","#587494","#537091","#4E6C8E","#49688B","#446488","#3F6085","#3A5C82","#35587F","#30547C","#2C5079","#274C76","#224873","#1D4470","#18406D","#133C6A","#0E3867","#093464","#053061","#053061"],
				    ["#F7F7F7","#F5F4F4","#F4F2F2","#F2EFF0","#F1EDEE","#EFEAEC","#EEE8E9","#ECE5E7","#EBE3E5","#E9E0E3","#E8DEE1","#E7DBDF","#E5D9DC","#E4D6DA","#E2D4D8","#E1D1D6","#DFCFD4","#DECCD1","#DCCACF","#DBC7CD","#D9C5CB","#D8C2C9","#D7C0C7","#D5BDC4","#D4BBC2","#D2B8C0","#D1B6BE","#CFB3BC","#CEB1B9","#CCAEB7","#CBACB5","#C9A9B3","#C8A7B1","#C7A4AF","#C5A2AC","#C49FAA","#C29DA8","#C19AA6","#BF98A4","#BE95A1","#BC939F","#BB909D","#B98E9B","#B88B99","#B78997","#B58694","#B48492","#B28190","#B17F8E","#AF7C8C","#AE7A89","#AC7787","#AB7585","#A97283","#A87081","#A76D7F","#A56B7C","#A4687A","#A26678","#A16376","#9F6174","#9E5E71","#9C5C6F","#9B596D","#99576B","#985469","#965266","#954F64","#944D62","#924A60","#91485E","#8F455C","#8E4359","#8C4057","#8B3E55","#893B53","#883951","#87364F","#85344C","#84314A","#822F48","#812C46","#7F2A44","#7E2741","#7C253F","#7B223D","#79203B","#781D39","#771B37","#751834","#741632","#721330","#71112E","#6F0E2C","#6E0C29","#6C0927","#6B0725","#690423","#680221","#67001F","#67001F"],
				    ["#F7F7F7","#F4F4F5","#F2F2F3","#EFF0F2","#EDEEF0","#EAECEF","#E8EAED","#E5E8EC","#E3E6EA","#E1E4E9","#DEE2E7","#DCE0E6","#D9DEE4","#D7DCE3","#D4DAE1","#D2D8E0","#CFD6DE","#CDD4DD","#CBD2DB","#C8D0DA","#C6CED8","#C3CCD7","#C1CAD5","#BEC8D4","#BCC6D2","#B9C4D1","#B7C2CF","#B4C0CE","#B2BECC","#B0BCCB","#ADBAC9","#ABB8C8","#A8B6C6","#A6B4C4","#A3B2C3","#A1B0C1","#9FAEC0","#9CACBE","#9AAABD","#97A8BB","#95A6BA","#92A4B8","#90A2B7","#8DA0B5","#8B9EB4","#899CB2","#869AB1","#8498AF","#8196AE","#7F94AC","#7C92AB","#7A90A9","#778EA8","#758CA6","#728AA5","#7088A3","#6E86A2","#6B84A0","#69829F","#66809D","#647E9C","#617C9A","#5F7A99","#5C7897","#5A7696","#587494","#557293","#537091","#506E8F","#4E6C8E","#4B6A8C","#49688B","#476689","#446488","#426286","#3F6085","#3D5E83","#3A5C82","#385A80","#35587F","#33567D","#30547C","#2E527A","#2C5079","#294E77","#274C76","#244A74","#224873","#1F4671","#1D4470","#1A426E","#18406D","#163E6B","#133C6A","#113A68","#0E3867","#0C3665","#093464","#073262","#053061","#053061"]
				]
				var letter = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"];
				var number = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]
				var field = ["1","2","3","4"]
				var value_array=[];
				letter.forEach(function(l) {
				    number.forEach(function(n) {
					field.forEach(function(f) {
					    if (a.contentDocument.getElementById(l+n+\'_\'+f)) {
						rx=/ShowTooltip\(evt\,\'(.+?)\'/g;
						value_array.push(rx.exec(a.contentDocument.getElementById(l+n+\'_\'+f).outerHTML)[1])
					    }
					});
				    });
				});
				if (document.getElementById(\'z-norm\').checked) {
				    var i=0;
				    var avg=mean(value_array);
				    var sd=variance(value_array);
					value_array.forEach(function(val) {
					    value_array[i]=(value_array[i]-avg)/sd;
					    i+=1;
					});
				}
				var i=0;
				letter.forEach(function(l) {
				    number.forEach(function(n) {
					field.forEach(function(f) {
					    if (a.contentDocument.getElementById(l+n+\'_\'+f)) {			
						a.contentDocument.getElementById(l+n+\'_\'+f).setAttribute("style","fill:"+color_vector[document.getElementById(\'color\').selectedIndex][Math.floor(100*((value_array[i]-Math.min.apply( Math, value_array ))/(Math.max.apply( Math, value_array )-Math.min.apply( Math, value_array ))))]+";");
						//console.log(Math.floor(100*((value_array[i]-Math.min.apply( Math, value_array ))/(Math.max.apply( Math, value_array )-Math.min.apply( Math, value_array )))))
						i+=1;
					    }
					});
				    });
				});	
			    }
			    function toggle_well_field(well){
				var a= document.getElementById(\'plate_plot\')
				if (document.getElementById(\'well_view\').value == \'well\') {
					if (a.contentDocument.activeElement.getElementById(\'wells\')) {
					    a.contentDocument.activeElement.getElementById(\'wells\').remove()
					}
					    var color_vector=[
						["#67001F","#690423","#6C0927","#6F0E2C","#721330","#751834","#781D39","#7B223D","#7E2741","#812C46","#84314A","#87364F","#893B53","#8C4057","#8F455C","#924A60","#954F64","#985469","#9B596D","#9E5E71","#A16376","#A4687A","#A76D7F","#A97283","#AC7787","#AF7C8C","#B28190","#B58694","#B88B99","#BB909D","#BE95A1","#C19AA6","#C49FAA","#C7A4AF","#C9A9B3","#CCAEB7","#CFB3BC","#D2B8C0","#D5BDC4","#D8C2C9","#DBC7CD","#DECCD1","#E1D1D6","#E4D6DA","#E7DBDF","#E9E0E3","#ECE5E7","#EFEAEC","#F2EFF0","#F5F4F4","#F4F4F5","#EFF0F2","#EAECEF","#E5E8EC","#E0E4E9","#DCE0E6","#D7DCE3","#D2D8E0","#CDD4DD","#C8D0DA","#C3CCD7","#BEC8D4","#B9C4D1","#B4C0CE","#B0BCCB","#ABB8C8","#A6B4C4","#A1B0C1","#9CACBE","#97A8BB","#92A4B8","#8DA0B5","#899CB2","#8498AF","#7F94AC","#7A90A9","#758CA6","#7088A3","#6B84A0","#66809D","#617C9A","#5C7897","#587494","#537091","#4E6C8E","#49688B","#446488","#3F6085","#3A5C82","#35587F","#30547C","#2C5079","#274C76","#224873","#1D4470","#18406D","#133C6A","#0E3867","#093464","#053061","#053061"],
						["#F7F7F7","#F5F4F4","#F4F2F2","#F2EFF0","#F1EDEE","#EFEAEC","#EEE8E9","#ECE5E7","#EBE3E5","#E9E0E3","#E8DEE1","#E7DBDF","#E5D9DC","#E4D6DA","#E2D4D8","#E1D1D6","#DFCFD4","#DECCD1","#DCCACF","#DBC7CD","#D9C5CB","#D8C2C9","#D7C0C7","#D5BDC4","#D4BBC2","#D2B8C0","#D1B6BE","#CFB3BC","#CEB1B9","#CCAEB7","#CBACB5","#C9A9B3","#C8A7B1","#C7A4AF","#C5A2AC","#C49FAA","#C29DA8","#C19AA6","#BF98A4","#BE95A1","#BC939F","#BB909D","#B98E9B","#B88B99","#B78997","#B58694","#B48492","#B28190","#B17F8E","#AF7C8C","#AE7A89","#AC7787","#AB7585","#A97283","#A87081","#A76D7F","#A56B7C","#A4687A","#A26678","#A16376","#9F6174","#9E5E71","#9C5C6F","#9B596D","#99576B","#985469","#965266","#954F64","#944D62","#924A60","#91485E","#8F455C","#8E4359","#8C4057","#8B3E55","#893B53","#883951","#87364F","#85344C","#84314A","#822F48","#812C46","#7F2A44","#7E2741","#7C253F","#7B223D","#79203B","#781D39","#771B37","#751834","#741632","#721330","#71112E","#6F0E2C","#6E0C29","#6C0927","#6B0725","#690423","#680221","#67001F","#67001F"],
						["#F7F7F7","#F4F4F5","#F2F2F3","#EFF0F2","#EDEEF0","#EAECEF","#E8EAED","#E5E8EC","#E3E6EA","#E1E4E9","#DEE2E7","#DCE0E6","#D9DEE4","#D7DCE3","#D4DAE1","#D2D8E0","#CFD6DE","#CDD4DD","#CBD2DB","#C8D0DA","#C6CED8","#C3CCD7","#C1CAD5","#BEC8D4","#BCC6D2","#B9C4D1","#B7C2CF","#B4C0CE","#B2BECC","#B0BCCB","#ADBAC9","#ABB8C8","#A8B6C6","#A6B4C4","#A3B2C3","#A1B0C1","#9FAEC0","#9CACBE","#9AAABD","#97A8BB","#95A6BA","#92A4B8","#90A2B7","#8DA0B5","#8B9EB4","#899CB2","#869AB1","#8498AF","#8196AE","#7F94AC","#7C92AB","#7A90A9","#778EA8","#758CA6","#728AA5","#7088A3","#6E86A2","#6B84A0","#69829F","#66809D","#647E9C","#617C9A","#5F7A99","#5C7897","#5A7696","#587494","#557293","#537091","#506E8F","#4E6C8E","#4B6A8C","#49688B","#476689","#446488","#426286","#3F6085","#3D5E83","#3A5C82","#385A80","#35587F","#33567D","#30547C","#2E527A","#2C5079","#294E77","#274C76","#244A74","#224873","#1F4671","#1D4470","#1A426E","#18406D","#163E6B","#133C6A","#113A68","#0E3867","#0C3665","#093464","#073262","#053061","#053061"]
					    ]
					    var letter = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"];
					    var number = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]
					    var field = ["1","2","3","4"]
					    var value_array=[];
					    var no_nan_array=[];
					    letter.forEach(function(l) {
						number.forEach(function(n) {
						    var field_array=[]
						    field.forEach(function(f) {
							if (a.contentDocument.getElementById(l+n+\'_\'+f)) {
							    rx=/ShowTooltip\(evt\,\'(.+?)\'/g;
							    field_array.push(rx.exec(a.contentDocument.getElementById(l+n+\'_\'+f).outerHTML)[1])
							}
						    });
						    value_array.push(mean(field_array))
						    if (!isNaN(mean(field_array))) {
							no_nan_array.push(mean(field_array))
						    }			    
						});
					    });
					    
					    var i=0;
					    var newgElement = document.createElementNS("http://www.w3.org/2000/svg", \'g\');
					    newgElement.setAttribute("id","wells");
					    letter.forEach(function(l) {
						number.forEach(function(n) {
						    if (!isNaN(value_array[i])) {
							var x=a.contentDocument.getElementById(l+n+\'_\'+1).x.baseVal.value
							var y=a.contentDocument.getElementById(l+n+\'_\'+1).y.baseVal.value
							var newElement = document.createElementNS("http://www.w3.org/2000/svg", \'rect\'); //Create a path in SVG\'s namespace
							newElement.setAttribute("x",x); //Set path\'s data
							newElement.setAttribute("y",y); //Set path\'s data
							newElement.setAttribute("width","36px");
							newElement.setAttribute("height","36px");
						       
							newElement.setAttribute("onmousemove","ShowTooltip(evt,\'"+value_array[i]+"\')");
							newElement.setAttribute("onmouseout","HideTooltip()");
							newElement.setAttribute("onmouseover","this.style.cursor=\'pointer\'");
							
							newElement.setAttribute("style","fill:"+color_vector[document.getElementById(\'color\').selectedIndex][Math.floor(100*((value_array[i]-Math.min.apply( Math, no_nan_array ))/(Math.max.apply( Math, no_nan_array )-Math.min.apply( Math, no_nan_array ))))]+";");
			
							newgElement.appendChild(newElement);
							a.contentDocument.activeElement.appendChild(newgElement);
							
							tooltip = a.contentDocument.getElementById(\'tooltip\');
							tooltip_bg = a.contentDocument.getElementById(\'tooltip_bg\');
							a.contentDocument.activeElement.appendChild(tooltip_bg);
							a.contentDocument.activeElement.appendChild(tooltip);
							
						    }
						i+=1;
						});    
					    });
					    a.contentDocument.activeElement.getElementById(document.getElementById(\'feature\').value).setAttributeNS(null,"visibility","hidden");
				    }else{
					if (a.contentDocument.activeElement.getElementById(\'wells\')) {
					    a.contentDocument.activeElement.getElementById(\'wells\').setAttributeNS(null,"visibility","hidden");
					    a.contentDocument.activeElement.getElementById(document.getElementById(\'feature\').value).setAttributeNS(null,"visibility","visible");
					}		
				    }
			    }
			    function mean(array) {
				    var sum = 0;
				    for (i = 0; i < array.length; i++)
					{
					    sum += parseFloat(array[i]); //don\'t forget to add the base
					}
				    var avg = sum/array.length;
				    return (avg)
			    }
			    function variance(array) {
				    var avg=mean(array);
				    var sum = 0;
				    for (i = 0; i < array.length; i++)
					{
					    sum += Math.pow( (parseFloat(array[i]) - avg ) , 2 ) ; //don\'t forget to add the base
					}
				    var vari = Math.pow((sum/array.length), 1/2 );
				    return (vari)
			    }
			</script>
			<style>
				.tooltip{
					font-size: 20px;
				}
				.tooltip_bg{
					fill: white;
					stroke: black;
					stroke-width: 1;
					opacity: 0.85;
				}
				.clickable { cursor: pointer;  }
				.clickable:hover { color: blue;  }
				.clickable:visited {
					color: red;
					font-family: Verdana;
					text-decoration:none;
				}
			</style>
			</head>
			    <body>
				<select name="feature" id="feature" onchange="change_plot(document.getElementById(\'feature\').value);">'."\n";
	opendir $outdir, $output_folder;
	foreach my $file (sort readdir($outdir)){
		if($file=~m/(.*)\.svg/ && !($file=~m/total.svg/ )){
			my $name=$1;
			print $output '<option value=\''.$name.'\' >'.$name.'</option>'."\n";			
		}
	}
	closedir $outdir;
	print $output '</select>
			<select name="well_view" id="well_view" onchange="toggle_well_field(document.getElementById(\'well_view\').value) ; adjust_colors();">
			    <option value=\'fields\' selected="selected" >per field</option>
			    <option value=\'well\' >per well</option>
			    </select>
			<select name="color" id="color" onchange="toggle_well_field(document.getElementById(\'well_view\').value) ; adjust_colors();">
			    <option value=\'divergent\' >divergent (red,white,blue)</option>
			    <option value=\'red\' >white to red</option>
			    <option value=\'blue\' selected="selected" >white to blue</option>	    
			</select>
			<input id="z-norm" style="vertical-align: top; cursor: pointer;" onMouseDown="this.__chk = this.checked" onClick="if (this.__chk) this.checked = false ; toggle_well_field(document.getElementById(\'well_view\').value) ; adjust_colors();" name="z-norm" value="z" type="radio"> normalize colors by z score 
		    <br>
		    <div id="plate_container">
			<object id="plate_plot" type="image/svg+xml" width="100%" height="1000" data="cells.svg" onload="toggle_well_field(document.getElementById(\'well_view\').value) ; adjust_colors();" ></object>
		    </div>
		    </body>
		</html>'."\n";
	close $output;
	return(1);      
}

#########################################################################################
#name:      update tile colors on plate
#function:  parses an svg value plate plot for the values and adjusts the colors 
#			according to the values
#input:     (/path/to/folder name [string] [string]) 
#output:    svg
#########################################################################################
sub adjust_colors_in_svg {
	my $output_folder=$_[0];
	my $name=$_[1];
	my @temp=(0,"inf");
	my @color_vector=('rgb(103,0,31)','rgb(178,24,43)','rgb(214,96,77)','rgb(244,165,130)','rgb(253,219,199)','rgb(247,247,247)','rgb(209,229,240)','rgb(146,197,222)','rgb(67,147,195)','rgb(33,102,172)','rgb(5,48,97)');
	open ($input, "<".$output_folder."/".$name.".svg") or die $!;		
		while(<$input>){
			if($_=~m/ShowTooltip\(evt\,\'(\-*\d+[\.\d+]*).*/){
				my $tempval=$1;				
				if($tempval>$temp[0]){					
					$temp[0]=$tempval;
				}elsif($tempval<$temp[1]){
					$temp[1]=$tempval;
				}
			}
		}
	close $input;
	open ($input, "<".$output_folder."/".$name.".svg") or die $!;		
	open $tempf ,">".$output_folder."/".$name."temp";
		while(<$input>){
			my $line=$_;
			my $color_val=0;
			my $color="rgb(0,0,0)";
			if($line=~m/ShowTooltip\(evt\,\'(\-*\d+[\.\d+]*).*/){
				my $val=$1;	
				if($val>$temp[0]){					
					$temp[0]=$val;
				}elsif($val<$temp[1]){
					$temp[1]=$val;
				}
				if(($temp[0]-$temp[1])!=0){			
					$color_val=10*(($val-$temp[1])/($temp[0]-$temp[1]));
				}else{
					$color_val=0
				}			
				$color=$color_vector[int($color_val)];
			}
			$line=~s/(style\=\"fill\:)rgb\(\d+\,\d+\,\d+\)/$1$color/;
			print $tempf $line;
		}		
		close $tempf;
	close $input;
	rename $output_folder."/".$name."temp",$output_folder."/".$name.".svg";
	return(1);      
}