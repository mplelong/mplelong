#! /usr/bin/perl                                                                                                                                          
use strict;
use warnings;

{ # main    SEE END OF FILE FOR COMMENTS & EXAMPLES                                                                                                                                      

 my $numArgs = $#ARGV + 1;
 if($numArgs != 10 )
 {
     print "usage: ./concat_XY_subplane.pl tslice p1 jproc data_dir data_root out_dir out_root myid iproc_start iproc_end \n";
     last;
 }
 
    my $tslice=$ARGV[0];
    my $p1=$ARGV[1];
    my $jproc=$ARGV[2];
    my $data_dir=$ARGV[3];     # e.g. ../output/3D
    my $data_root=$ARGV[4];    # e.g. XYZ
    my $out_dir=$ARGV[5];      # e.g. ../output/slices/3D
    my $out_root=$ARGV[6];     # e.g. slice
    my $myid=$ARGV[7];         # mpi id of processor running the this script
    
    my $iproc_start=$ARGV[8];     # 1st iproc value in subplane
    my $iproc_end=$ARGV[9];       # last iproc value in subplane
	  
    my $NCKS='ncks';
    my $NCPDQ='ncpdq';
    my $NCRCAT='ncrcat';
    
    my $char_myid=sprintf("%03d",$myid);
    my $char_j=sprintf("%03d",$jproc);     # fixed "j" processor id  
    my $istep=sprintf("%06d",$tslice);
    
    my ($pid,$char_i,$ncfile,$tmpfile1,$tmpfile2,$newfile,$outfile,$cmd,$tmp_dir,$tmpfile);
    
    #----------------------------------------------------------------------------------------------
    #   add trailing underscores so input is consistent w/ filename roots in io_params
    #----------------------------------------------------------------------------------------------
    $data_root = $data_root . "_";
    $out_root  = $out_root . "_";
    
    #----------------------------------------------------------------------------------------------
    #   create a tmp directory to work in, unique for myid
    #----------------------------------------------------------------------------------------------
    $tmp_dir = $out_dir . "/" . $char_myid . "_" . $istep;
    mkdir $tmp_dir,0755 if ! -d $tmp_dir;
    
    #-----------------------------------------                                                                                                                                                          
    #  1d decomposition w/ p1=1
    #  entire plane in 1 file
    #  just extract the requested time slice
    #----------------------------------------- 
    if($p1 == 1){
      $pid=0;
      $char_i = sprintf("%03d", $pid);	 
	  $ncfile  = $data_dir . "/" . $data_root . $char_i . "-" . $char_j . ".nc";
	  $outfile = $out_dir . "/" . $out_root . $istep . ".nc";
	  $cmd="$NCKS -O -d timedimension,$istep $ncfile $outfile >& /dev/null ";
	  system($cmd);
	  $cmd="rm -rf " . $tmp_dir;
      system($cmd);
	  last;	
    } 
    
    #--------------------------------------------------                                                                                                                                                          
    #  Loop over pid from $iproc_start to $iproc_end 
    #  swap record dimension from k to i 
    #--------------------------------------------------                                                                                                                                                            
    for($pid=$iproc_start; $pid<$iproc_end+1; $pid++){
	   $char_i = sprintf("%03d", $pid);	 
	   $ncfile  = $data_dir . "/" . $data_root . $char_i . "-" . $char_j . ".nc";
	   $tmpfile = $tmp_dir . "/" . $myid . "_" . $char_i . "-" . $char_j . ".nc";
	   #-----------------------------------------                                                                                                                                                          
       #  extract the requested time slice
       #  store result in $tmpfile1
       #-----------------------------------------
       $cmd="$NCKS -O -d timedimension,$istep $ncfile $tmpfile >& /dev/null";
	   system($cmd);
	   
	   #-----------------------------------------                                                                                                                                                          
       #  swap record dimension from time to i
       #  overwrite $tmpfile
       #-----------------------------------------
	   $cmd="$NCPDQ -O -a idimension,timedimension $tmpfile $tmpfile >& /dev/null";
	   system($cmd);
    }

    #------------------------------------                                                                                                                                                          
    #  concatenate the iproc files 
    #------------------------------------
    $outfile = $out_dir . "/" . $out_root . $istep . ".nc";
    $cmd="$NCRCAT -O " . $tmp_dir . "/" . $myid . '_*' . "-" . $char_j . ".nc" . " $outfile";
    #print " $cmd \n";
    system($cmd);
    
    #-----------------------------------------                                                                                                                                                          
    #  swap dimensions around, t=record
    #  overwrite $outfile
    #-----------------------------------------
	$cmd="$NCPDQ -O -a timedimension,jdimension,idimension $outfile $outfile >& /dev/null";
	system($cmd);
  
    #------------------------------------                                                                                                                                                          
    #  clean up the temp files 
    #------------------------------------
    $cmd="rm -f " . $tmp_dir . "/" . $myid . '_*' . "-" . $char_j . "_*.nc";
    system($cmd);

    #---------------------------------------------------------
    # get rid of the tmp directory
    #---------------------------------------------------------
    $cmd="rm -rf " . $tmp_dir;
    system($cmd);    
}   #  end 



