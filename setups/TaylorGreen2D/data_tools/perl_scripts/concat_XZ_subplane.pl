#! /usr/bin/perl                                                                                                                                          
use strict;
use warnings;

{ # main                                                                                                                                          

 my $numArgs = $#ARGV + 1;
 if($numArgs != 12 )
 {
     print "usage: ./concat_XZ_subplane.pl tslice p1 p2 data_dir data_root out_dir out_root myid iproc_start iproc_end jproc_start jproc_end \n";
     last;
 }

    my $tslice=$ARGV[0];
    my $p1=$ARGV[1];
    my $p2=$ARGV[2];
    my $data_dir=$ARGV[3];     # e.g. ../output/3D
    my $data_root=$ARGV[4];    # e.g. XYZ
    my $out_dir=$ARGV[5];      # e.g. ../output/slices/3D
    my $out_root=$ARGV[6];     # e.g. slice
    my $myid=$ARGV[7];         # mpi id of processor running the this script
    
    my $iproc_start=$ARGV[8];      # 1st jproc value in subplane
    my $iproc_end=$ARGV[9];        # last jproc value in subplane
    my $jproc_start=$ARGV[10];     # 1st jproc value in subplane
    my $jproc_end=$ARGV[11];       # last jproc value in subplane
     
    my $NCKS='ncks';
    my $NCPDQ='ncpdq';
    my $NCRCAT='ncrcat';
    
    my $istep=sprintf("%06d",$tslice);
    my ($ipid,$jpid,$char_i,$char_j,$ncfile,$jfile,$tmpfile,$newfile,$cmd);
    
    #----------------------------------------------------------------------------------------------
    #   add trailing underscores so input is consistent w/ filename roots in io_params
    #----------------------------------------------------------------------------------------------
    $data_root = $data_root . "_";
    $out_root  = $out_root . "_";
    
    
    #----------------------------------------------------
    #  Loop over jprocs from $jproc_start to $jproc_end 
    #----------------------------------------------------    
    for($jpid=$jproc_start; $jpid<$jproc_end+1; $jpid++){ 
     $char_j = sprintf("%03d", $jpid);   #  "j" processor id
     $jfile= $out_dir . "/" . "stitched_" . $myid ."-" . $char_j . ".nc";
     
     #----------------------------------------------------
     #  Loop over iprocs from $iproc_start to $iproc_end 
     #----------------------------------------------------
     for($ipid=$iproc_start; $ipid<$iproc_end+1; $ipid++){
      $char_i = sprintf("%03d", $ipid);   #  "i" processor id
     
      # raw data file with many slices
      $ncfile = $data_dir . "/" . $data_root . $char_i . "-" . $char_j . ".nc"; 
      # tmp file with only the desired time slice     
      $tmpfile= $out_dir . "/" . "tmp_" . $myid ."-" . $char_i . "-" . $char_j . ".nc";
    
      #------------------------------------------------------------
      # extract the requested time slice and write into a tmp file 
      #------------------------------------------------------------
      $cmd="$NCKS -O -d timedimension,$istep $ncfile $tmpfile";
      system($cmd);
 
      #---------------------------------------
      # swap record dimension from time to i
      #---------------------------------------
      $cmd="$NCPDQ -O -a idimension,kdimension,timedimension $tmpfile $tmpfile >& /dev/null";    
      system($cmd);
     }   # end loop over iproc files at fixed jproc
     
     #------------------------------------------------------------
     #  concatenate across the subset of iprocs, store in $jfile
     #------------------------------------------------------------
     $cmd="$NCRCAT -O " . $out_dir . "/tmp_" . $myid . '-*-' . $char_j . ".nc $jfile";
     system($cmd);
   
     #-----------------------------------------------------------------
     #  swap dimensions, make k be record dimension, overwrite $jfile
     #-----------------------------------------------------------------
     $cmd="$NCPDQ -O -a kdimension,idimension,timedimension $jfile $jfile >& /dev/null";    
     system($cmd);
   
     #---------------------------------------------------------
     # clean up the single time slice files at this j 
     #---------------------------------------------------------
     for($ipid=$iproc_start; $ipid<$iproc_end+1; $ipid++){  
       $char_i = sprintf("%03d", $ipid);   #  "i" processor id    
       $tmpfile= $out_dir . "/" . "tmp_" . $myid ."-" . $char_i . "-" . $char_j . ".nc";
    
       $cmd="rm -f $tmpfile";
       system($cmd);            
     }  # end cleanup loop over jproc files
            
  }    # end loop over jproc files
  
  #-------------------------------------------------------------
  # create the final output file name with the time index in it
  #-------------------------------------------------------------
  $newfile = $out_dir . "/" . $out_root .  $istep . ".nc";
  
  #---------------------------------------------------------
  # cat the subset of jproc files together
  #---------------------------------------------------------
  $cmd="$NCRCAT -O " . $out_dir . "/stitched_" . $myid . '-*' . " $newfile";
  system($cmd);
  
  #------------------------------------------------
  # now swap record dimension from k back to time
  #------------------------------------------------
  $cmd="$NCPDQ -O -a timedimension,kdimension,idimension $newfile $newfile";
  system($cmd);
  
  #---------------------------------------------------------
  # clean up the tmp files
  #---------------------------------------------------------
  $cmd="rm -f " . $out_dir . "/stitched_" . $myid . '-*' . "nc";
  system($cmd);
      
}   #  end 


