#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Capture::Tiny ':all';
use File::Find;

my $infolder="/data/raw_data/HC1098";

print "Watcher startet \@\: ".$infolder."\n";

my %platesdone;

my $totalplatenum=1;
my $curplatenumber=0;
while ($curplatenumber < $totalplatenum) {
    opendir(my $bigindir, $infolder);
	FORITER: foreach my $subfolder (sort readdir($bigindir)){
	    if (-d $infolder."/".$subfolder 
		  #&& -r $infolder."/".$subfolder
			&& !($subfolder=~m/^\./)
			&& !($subfolder=~m/Plate/)
		  && ($subfolder=~m/HC10984f/)
			&& !exists($platesdone{$infolder."/".$subfolder})  # uncomment this line to allow REMATCH to iterate only once over each plate
		    ) 
		{	
			my $outfolder="/data/results/HC1098/".$subfolder;
			open my $logfile ,">", "/data/scripts/HC1098/run-dir/log_file".$subfolder.".txt";

			if (!(-d $outfolder)) {
			    system('mkdir '.$outfolder);
			    print $logfile $outfolder."\n";
				start_off_image_analysis ("/data/scripts/HC1098/R-pipelines/HC1098_imageAnalysis.R",$infolder."/".$subfolder,$outfolder,384,"DAPI.tif","Cy5.tif","FITC.tif",40,$logfile);
			    $platesdone{$infolder."/".$subfolder}++;
			    print $logfile $infolder."/".$subfolder."\tdone\n";
			}else{
			    print $logfile $outfolder."\n";
				start_off_image_analysis ("/data/scripts/HC1098/R-pipelines/HC1098_imageAnalysis.R",$infolder."/".$subfolder,$outfolder,384,"DAPI.tif","Cy5.tif","FITC.tif",40,$logfile);
				$platesdone{$infolder."/".$subfolder}++;
			    print $logfile $infolder."/".$subfolder."\tdone\n";
			}
			close $logfile;
		}else{
			next FORITER;
		}	
	}
	sleep 1;
	closedir($bigindir);
}


#########################################################################################
#name:      start image analysis
#function:  closes an browser readable svg in a given file if closed already than unclose
#input:     (path/to/filename string) 
#output:    a closed or open svg
#########################################################################################
sub start_off_image_analysis {
	my $path_to_script=$_[0];
	my $input_folder=my $barcode= $_[1];
	$barcode=~s/.+\/(\S+)_\d{4}.+$/$1/g;
	my $output_folder=$_[2];
	my $format=$_[3];
	my $nuclei_scheme=$_[4];	#fld1wvDAPIDAPI.tif
	my $body_scheme=$_[5];		#fld1wvCy3Cy3.tif
	my $extra_scheme=$_[6];	
	my $tilesize=$_[7];		
	my $logfile_int=$_[8];
	my %Q=();
	my $start = time;
	print $logfile_int "hello world 1\n";
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){		
		foreach my $number ("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"){
			foreach my $field (1..4){
				my @temp=($letter,$number,$field);
				$Q{$letter.$number."_".$field} = \@temp;
			}
		}
	}
	print $logfile_int "hello world\n";
        print $logfile_int $input_folder."\n";
        
	if (! -e $output_folder."/total.html") {
		system('echo \'perl /data/scripts/HC1098/perl-schedulers/waiter.quick.pl '.$barcode.' '.$output_folder.' '.$tilesize.' '.$format.'\' | qsub -l walltime=5:00:00 -N '.$output_folder.' &'); #.'\'  '
	}
	print $logfile_int "I am before Queing\n";
	my %smQ;
	my $test=0;
	my $counter=0;
	my $test_injob=0;
	my $queue_length=0;
      QITER: while(%Q){
	    print $logfile_int keys(%Q)."\n";
	          find({wanted => sub {
                                        if(-d $_ && $_=~m/thumbs/){
                                            $File::Find::prune = 1;
                                            return;
                                        }
                                        if($_=~m/(.+\/(\w+)_\S+?\/)(.+_t1)\.tif/ && -r $_){
                                        my $path=$1;
                                        my $bc=$2;
                                        my $newname=$3;
                                        $newname=~s/^(\w)_(\d+)_/$1$2_/g;
                                        $newname=~s/_f(\d)_/_$1_/g;
                                        $newname=~s/_c\d_/_/g;
                                        $newname=~s/xDAPI_mDAPI/DAPI/g;
                                        $newname=~s/xFITC_mFITC/FITC/g;
                                        $newname=~s/xCy5_mCy5/Cy5/g;
                                        $newname=~s/xCy3_mCy3/Cy3/g;
                                        $newname=~s/_z1_t1//g;
                                        rename($_,$path."/".$bc."_$newname.tif");
                                        }
                                            },
              no_chdir => 1},$input_folder);
		    $counter=0;
				CHECKITER: foreach my $key (keys(%Q)){
					my $letter=@{$Q{$key}}[0];
					my $number=@{$Q{$key}}[1];
					my $field=@{$Q{$key}}[2];
					if(!(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab") &&
					   !(-e $output_folder."/".$barcode."_".$letter.$number."_".$field."_segmented.png")) {
					print $logfile_int "searching ".$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme."\n";
                        if(	 (-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
                            && (-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme) 
                            && (-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme)
                            && (-r $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme)
                             && (-r $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme)
                             && (-r $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme)
                            && !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
                            && !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme) 
                            && !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme)
                        ){
                            print $logfile_int "found ".$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme."\n";
                            #sleep 1;
                                $queue_length=capture_stdout {system("qstat -B | grep vm")};					
                                $queue_length=~s/\S+\s+\S+\s+\S+\s+(\d+).*$/$1/;
                                #print $logfile $output_folder."/".$barcode."_".$letter.$number."_".$field."_single_cell.RData\n";
                                print $logfile_int "to analyze ".$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme."\n";
                                if ($queue_length< 300) {
                                    $smQ{$key}=$Q{$key};						
                                    $test++;
                                    $counter++;
                                    if (scalar(keys(%smQ)) >= 16) {								
                                        my $cmd_start='R -f '.$path_to_script.' --slave --args';
                                        my $filename="";
                                        foreach my $smkey (keys(%smQ)){
                                            delete $Q{$smkey};
                                            $filename=$filename." ".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$nuclei_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$body_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$extra_scheme;
                                        $test_injob++;
                                        }
                                        my $cmd_tail=$output_folder;
                                        system('echo "'.$cmd_start.$filename.' '.$cmd_tail.'" | qsub -l walltime=140:00:00 -N HC1098');# 
                                        sleep 1;
                                        %smQ=();
                                        next QITER;										
                                    }elsif(scalar(keys(%Q)) == scalar(keys(%smQ))){
                                        my $cmd_start='R -f '.$path_to_script.' --slave --args';
                                        my $filename="";
                                        foreach my $smkey (keys(%smQ)){
                                            delete $Q{$smkey};
                                            $filename=$filename." ".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$nuclei_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$body_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$extra_scheme;
                                        $test_injob++;
                                        }
                                        my $cmd_tail=$output_folder;
                                        system('echo "'.$cmd_start.$filename.' '.$cmd_tail.'" | qsub -l walltime=140:00:00 -N HC1098 '); #
                                        sleep 1;
                                        %smQ=();
                                        next QITER;
                                    }else{
                                        next CHECKITER;
                                    }
                                    
                                }else{
                                     print $logfile_int "have not found ".$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme."\n";
                                    next CHECKITER;
                                }
                            }elsif((!-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
                            && (!-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme) 
                            && (!-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme)){
                               print $logfile_int " cant find ".$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme."\n"; 
                            }				
						
					}else{
						delete $Q{$key};
					}
					#sleep 1;
				}
				#sleep 1;
			}
		print "I am after Queing\n";
		if ( scalar(keys(%smQ)) > 0) {
			#do the raw analysis either with R ,CP, or Hcell
			my $cmd_start='R -f '.$path_to_script.' --slave --args';
			my $filename="";
			foreach my $smkey (keys(%smQ)){
				$filename=$filename." ".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$nuclei_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$body_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$extra_scheme;
			$test_injob++;
			}
			my $cmd_tail=$output_folder;
			system('echo "'.$cmd_start.$filename.' '.$cmd_tail.'" | qsub -l walltime=140:00:00 -N HC1098 '); #
			sleep 1; 
		}
}
		
