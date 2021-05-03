#! /usr/bin/perl                                                                                                                                          
use strict;
use warnings; 

{ # main                                                                                                                                          

 my $numArgs = $#ARGV + 1;
 if($numArgs != 7 )
 {
     print "usage: ./concat_profiles.pl iproc p2 data_dir data_root out_dir out_root myid \n";
     print "example: ./concat_profiles.pl 0 144 '../../output/1D' 'C3' '../../output/slices/1D' 'C3' 0 \n";
     last;
 }

    my $iproc=$ARGV[0];
    my $p2=$ARGV[1];
    my $data_dir=$ARGV[2];     # e.g. ../output/3D
    my $data_root=$ARGV[3];    # e.g. XYZ
    my $out_dir=$ARGV[4];      # e.g. ../output/slices/3D
    my $out_root=$ARGV[5];     # e.g. slice
    my $myid=$ARGV[6];         # mpi id of processor running the this script
     
    my $NCKS='ncks';
    my $NCPDQ='ncpdq';
    my $NCRCAT='ncrcat';
    
    my $char_i=sprintf("%03d",$iproc);     # fixed "i" processor id   (0 when p1=1)
    my ($pid,$char_j,$ncfile,$tmpfile,$newfile,$cmd);
    
    #----------------------------------------------------------------------------------------------
    #   add trailing underscores so input is consistent w/ filename roots in io_params
    #----------------------------------------------------------------------------------------------
    $data_root = $data_root . "_";
    $out_root  = $out_root ;
    
    
    #---------------------------------------
    #  Loop over pid from 0 to p2-1 
    #---------------------------------------    
    for($pid=0; $pid<$p2; $pid++){
  
     $char_j = sprintf("%03d", $pid);   #  "j" processor id
     
     $ncfile = $data_dir . "/" . $data_root . $char_i . "-" . $char_j . ".nc";
     $tmpfile= $out_dir . "/" . "tmp_" . $myid ."-" . $char_j . ".nc";
     
     #---------------------------------------
     # swap record dimension from time to k
     #---------------------------------------
      $cmd="$NCPDQ -O -a kdimension,timedimension $ncfile $tmpfile >& /dev/null"; 
      #print "$cmd \n";   
      system($cmd);
            
  }    # end loop over p2 files
  
  #---------------------------------------------------------
  # create an output file name 
  #---------------------------------------------------------
  $newfile = $out_dir . "/" . $out_root . ".nc";
  
  #---------------------------------------------------------
  # cat all p2 files at this i together, store in newfile
  #---------------------------------------------------------
  $tmpfile = $out_dir . "/" . "tmp.nc";
  $cmd="$NCRCAT -O " . $out_dir . "/tmp_" . $myid . '-*' . " $newfile";
  #print "$cmd \n";
  system($cmd); 
  
  
  
  #---------------------------------------------------------
  # clean up the p2 tmp files
  #---------------------------------------------------------
  for($pid=0; $pid<$p2; $pid++){  
     $char_j = sprintf("%03d", $pid);    
     $tmpfile= $out_dir . "/" . "tmp_" . $myid ."-" . $char_j . ".nc";
    
     $cmd="rm -f $tmpfile";
     #print "$cmd \n";
     system($cmd);            
  }  # end cleanup loop over p2 files
    
}   #  end 



