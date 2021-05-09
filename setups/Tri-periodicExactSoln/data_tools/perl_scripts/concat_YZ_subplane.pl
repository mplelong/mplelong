#! /usr/bin/perl                                                                                                                                          
use strict;
use warnings;

{ # main                                                                                                                                          

 my $numArgs = $#ARGV + 1;
 if($numArgs != 10 )
 {
     print "usage: ./concat_YZ_subplane.pl tslice iproc p2 data_dir data_root out_dir out_root myid jproc_start jproc_end \n";
     last;
 }

    my $tslice=$ARGV[0];
    my $iproc=$ARGV[1];
    my $p2=$ARGV[2];
    my $data_dir=$ARGV[3];     # e.g. ../output/3D
    my $data_root=$ARGV[4];    # e.g. XYZ
    my $out_dir=$ARGV[5];      # e.g. ../output/slices/3D
    my $out_root=$ARGV[6];     # e.g. slice
    my $myid=$ARGV[7];         # mpi id of processor running the this script
    
    my $jproc_start=$ARGV[8];     # 1st jproc value in subplane
    my $jproc_end=$ARGV[9];       # last jproc value in subplane
     
    my $NCKS='ncks';
    my $NCPDQ='ncpdq';
    my $NCRCAT='ncrcat';
    
    my $char_i=sprintf("%03d",$iproc);     # fixed "i" processor id   (0 when p1=1)
    my $istep=sprintf("%06d",$tslice);
    my ($pid,$char_j,$ncfile,$tmpfile,$newfile,$cmd);
    
    #----------------------------------------------------------------------------------------------
    #   add trailing underscores so input is consistent w/ filename roots in io_params
    #----------------------------------------------------------------------------------------------
    $data_root = $data_root . "_";
    $out_root  = $out_root . "_";
    
    
    #------------------------------------------------
    #  Loop over pid from $jproc_start to $jproc_end 
    #------------------------------------------------    
    for($pid=$jproc_start; $pid<$jproc_end+1; $pid++){
  
     $char_j = sprintf("%03d", $pid);   #  "j" processor id
     
     $ncfile = $data_dir . "/" . $data_root . $char_i . "-" . $char_j . ".nc";
     $tmpfile= $out_dir . "/" . "tmp_" . $myid ."-" . $char_j . ".nc";
    
     #---------------------------------------
     # extract the requested time slice
     # and write into a tmp file with both
     # the j index and myid in the filename
     #---------------------------------------
     $cmd="$NCKS -O -d timedimension,$istep $ncfile $tmpfile";
     system($cmd);
 
     #---------------------------------------
     # swap record dimension from time to k
     #---------------------------------------
      $cmd="$NCPDQ -O -a kdimension,timedimension $tmpfile $tmpfile >& /dev/null";    
      system($cmd);
            
  }    # end loop over jproc files
  
  #---------------------------------------------------------
  # create an output file name with the i index in it
  #---------------------------------------------------------
  $newfile = $out_dir . "/" . $out_root .  $istep . ".nc";
  
  #---------------------------------------------------------
  # cat all p2 files at this i together
  #---------------------------------------------------------
  $cmd="$NCRCAT -O " . $out_dir . "/tmp_" . $myid . '-*' . " $newfile";
  system($cmd);
  
  #------------------------------------------------
  # now swap record dimension from k back to time
  #------------------------------------------------
  $cmd="$NCPDQ -O -a timedimension,kdimension $newfile $newfile";
  system($cmd);
  
  #---------------------------------------------------------
  # clean up the tmp files
  #---------------------------------------------------------
  for($pid=0; $pid<$p2; $pid++){  
     $char_j = sprintf("%03d", $pid);    
     $tmpfile= $out_dir . "/" . "tmp_" . $myid ."-" . $char_j . ".nc";
    
     $cmd="rm -f $tmpfile";
     system($cmd);            
  }  # end cleanup loop over jproc files
    
}   #  end 


#---------------------------------------------------------------------------------------
#  (1)  ./create_global_snapshot.pl tslice p1 p2 data_dir data_root out_dir out_root myid
#        use for any concatenating needed for 1d decomposition w/ p1=1
#        use for XZplanes and XYZ blocks for 2d decomposition p1>1 and p2>1
#
#  (2)  ./concat_YZ.pl
#        use for YZ planes for 2d decomposition
#
#  (3)  ./concat_XY.pl
#        use for XY planes for 2d decomposition
#
#   Scripts assume that  output/slices  output/slices/2D and output/slices/3D exist.
#   (I create these directories in the "make outdirs" block in Makefile). Scripts are
#   intended to live in the "input" directory so they can be fine-tuned for different
#   sets of simulations.
#
#   NB the conventions regarding trailing slashes and underscores 
#      when directories and filename roots: i.e. don't add either in calling cmd.
#
#   The XYZ global files get bad names, e.g.  XYZ_30_000000.nc, 
#    rename the result in calling script if desired
#    e.g. 30 is the time slice but 000000 is the unimportant slice number in "new" files
#
#   All tmp files have myid in the filenames so these scripts can be called in parallel
#   from an MPI4PY script, distributing the time steps across processors, i.e. different
#   processors won't try to use the same tmp file names.
#---------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------
#     1D DATA DECOMPOSITION EXAMPLES & TESTS:
#
#        p1=1   p2=4  ==> np=4         test case  (nx,ny,nz) = (5,513,129)
#------------------------------------------------------------------------------------------------------------

 #  NB  2D files are "append", w/ many time slices (starting with 0)  
   #==> ./create_global_snapshot.pl 2 1 4 '../../../output/2D' 'XZplane' '../../../output/slices/2D' 'XZ' 0
 
 #  NB  2D files are "append", w/ many time slices (starting with 0)  
   #==> ./create_global_snapshot.pl 0 1 4 '../../../output/2D' 'YZplane' '../../../output/slices/2D' 'YZ' 0

 #  NB  3D files are "new", w/ only one slice, request time slice 0, specify which slice in filename
   #==> ./create_global_snapshot.pl 0 1 4 '../../../output/3D' 'XYZ_000030' '../../../output/slices/3D' 'XYZ_30' 0

 #  XY planes entirely contained on single  "j" processor when p1=1
 #  e.g.  output/2D/XYplane_iii-003.nc
 #
 # ./concat_XY.pl tslice p1 jproc data_dir data_root out_dir out_root myid
   #==> ./concat_XY.pl 1 2 3 '../../../output/2D' 'XYplane' '../../../output/slices/2D' 'XY' 0
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------
#     2D DATA DECOMPOSITION TESTS:
#
#        p1=2   p2=4  ==> np=8        test case  (nx,ny,nz) = (5,513,129)
#------------------------------------------------------------------------------------------------------------
 #  NB  2D files are "append", w/ many time slices (starting with 0) 
 #        XZ planes involve all p1 x p2 processors 
   #==> ./create_global_snapshot.pl 2 2 4 '../../../output/2D' 'XZplane' '../../../output/slices/2D' 'XZ' 0


 #  NB  3D files are "new", w/ only one slice, request time slice 0, specify which slice in filename
 #        XYZ blocks involve all p1 x p2 processors
   #==> ./create_global_snapshot.pl 0 2 4 '../../../output/3D' 'XYZ_000030' '../../../output/slices/3D' 'XYZ_30' 0
   #==> mv ../../../output/slices/3D/XYZ_30_000000.nc ../../../output/slices/3D/XYZ_30.nc


 #  YZ planes have only a single "i" processor and all "j" processors
 #  e.g.  output/2D/YZplane_001-00*.nc
 #
 # ./concat_YZ.pl tslice iproc p2 data_dir data_root out_dir out_root myid
   #==> ./concat_YZ.pl 2 1 4 '../../../output/2D' 'YZplane' '../../../output/slices/2D' 'YZ' 0

 
 #  XY planes have only a single "j" processor and all "i" processors
 #  e.g.  output/2D/XYplane_iii-003.nc
 #
 # ./concat_XY.pl tslice p1 jproc data_dir data_root out_dir out_root myid
   #==> ./concat_XY.pl 2 2 3 '../../../output/2D' 'XYplane' '../../../output/slices/2D' 'XY' 0
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

