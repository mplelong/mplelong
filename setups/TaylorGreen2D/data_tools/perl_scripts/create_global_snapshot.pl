#! /usr/bin/perl
use strict;
use warnings;

{ # main     SEE END OF FILE FOR COMMENTS & EXAMPLES

 my $numArgs = $#ARGV + 1;
 if($numArgs != 8 )
 {
     print "usage: ./create_global_snapshot.pl tslice p1 p2 data_dir data_root out_dir out_root myid \n";
     last;
 }

    my $tslice=$ARGV[0];       # integer time slice, first slice in a file is "0"
    my $p1=$ARGV[1];           # processors used to split x dimension
    my $p2=$ARGV[2];           # processors used to split y & z dimensions
    my $data_dir=$ARGV[3];     # where the simulation output data lives, e.g. ../output/3D
    my $data_root=$ARGV[4];    # e.g. XYZ
    my $out_dir=$ARGV[5];      # where to store the new files, e.g. ../output/slices/3D
    my $out_root=$ARGV[6];     # e.g. slice
    my $myid=$ARGV[7];         # mpi id of processor running the this script, set to 0 if no MPI used
    my $char_myid=sprintf("%03d",$myid);
    my ($tmp_dir,$cmd);

    #------------------------------------------------------------------------------------
    #   add trailing underscores so input is consistent w/ filename roots in io_params
    #------------------------------------------------------------------------------------
    $data_root = $data_root . "_";
    $out_root  = $out_root . "_";


    #----------------------------------------------------------------------------------------------
    #   create a tmp directory to work in, unique for myid
    #----------------------------------------------------------------------------------------------
    $tmp_dir = $out_dir . "/" . $char_myid;
    mkdir $tmp_dir,0755 if ! -d $tmp_dir;


    #----------------------------------------------------------------------------------------------
    #  for each (x) i processor 0,1,...p1-1 , concatenate across the z (p2) direction
    #  this concatenation requires swapping of record indices
    #   ==> this creates a set of p1 temp files
    #----------------------------------------------------------------------------------------------
    &cat_slices_2d_decomp( $tslice, $p1, $p2, $data_dir, $data_root, $tmp_dir, $out_root, $myid );

    #----------------------------------------------------------------------------------------------
    #  concatenate across the p1 temp files
    #   ==> this creates a set of p1 temp files
    #----------------------------------------------------------------------------------------------
    &concat_across_x( $tslice, $p1, $tmp_dir, $out_root, $myid );

    #----------------------------------------------------------------------------------------------
    #  mv the concatenated file from myid's tmp directory to the specified output directory
    #----------------------------------------------------------------------------------------------
    $cmd = "mv $tmp_dir/XYZ* $out_dir/. ";
    system($cmd);

    #---------------------------------------------------------
    # now safely get rid of myid's tmp directory
    #---------------------------------------------------------
    $cmd="rm -rf $tmp_dir";
    system($cmd);

}  # end of main



sub cat_slices_2d_decomp
{
 use strict;
 use warnings;
 my $numArgs = @_;
 if($numArgs != 8 )
 {
     print "usage: cat_slices_2d_decomp tslice p1 p2 data_dir data_root out_dir out_root myid \n";
     last;
 }

    my $tslice=$_[0];
    my $p1=$_[1];
    my $p2=$_[2];
    my $data_dir=$_[3];     # e.g. ../output/3D
    my $data_root=$_[4];    # e.g. XYZ_
    my $out_dir=$_[5];      # e.g. ../output/slices/3D
    my $out_root=$_[6];     # e.g. slice_
    my $myid=$_[7];         # mpi id of processor running the this script
    my ($i);


    #------------------------------------------
    #  loop through each of the p1 x slices
    #------------------------------------------
    for($i=0;$i<$p1;$i++){
      &concat_across_z($tslice,$p1,$p2,$i,$data_dir,$data_root,$out_dir,$out_root,$myid);
    }

}   #  end subroutine cat_slices_2d_decomp



sub concat_across_z
{
 use strict;
 use warnings;
 my $numArgs = @_;
 if($numArgs != 9 )
 {
     print "usage: concat_across_z tslice p1 p2 i data_dir data_root out_dir out_root myid \n";
     print " $_[0] \n";
     print " $_[1] \n";
     print " $_[2] \n";
     print " $_[3] \n";
     print " $_[4] \n";
     print " $_[5] \n";
     print " $_[6] \n";
     print " $_[7] \n";
     print " $_[8] \n";
     last;
 }

    my $tslice=$_[0];
    my $p1=$_[1];
    my $p2=$_[2];
    my $i=$_[3];
    my $data_dir=$_[4];     # e.g. ../output/3D
    my $data_root=$_[5];    # e.g. XYZ_
    my $out_dir=$_[6];      # e.g. ../output/slices/3D
    my $out_root=$_[7];     # e.g. slice_
    my $myid=$_[8];         # mpi id of processor running the this script

    my $NCKS='ncks --64bit_offset';
    my $NCPDQ='ncpdq --64bit_offset';
    my $NCRCAT='ncrcat --64bit_offset';

    my $char_i=sprintf("%03d",$i);
    my $istep=sprintf("%06d",$tslice);
    my ($pid,$char_j,$ncfile,$tmpfile,$newfile,$cmd);


    #---------------------------------------
    #  Loop over pid from 0 to p2-1
    #---------------------------------------
    for($pid=0; $pid<$p2; $pid++){

     $char_j = sprintf("%03d", $pid);

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

  }    # end loop over p2 files

  #---------------------------------------------------------
  # create an output file name with the i index in it
  #---------------------------------------------------------
  $newfile = $out_dir . "/" . "slice_$istep" . "_" . $char_i . "_" . $myid . ".nc";

  #---------------------------------------------------------
  # cat all p2 files at this i together
  #---------------------------------------------------------
  $cmd="$NCRCAT -O " . $out_dir . "/tmp_" . $myid . '-*' . " $newfile";
  system($cmd);

  #---------------------------------------------------------
  # clean up the p2 tmp files
  #---------------------------------------------------------
  for($pid=0; $pid<$p2; $pid++){
     $char_j = sprintf("%03d", $pid);
     $tmpfile= $out_dir . "/" . "tmp_" . $myid ."-" . $char_j . ".nc";

     $cmd="rm -f $tmpfile";
     system($cmd);
  }  # end cleanup loop over p2 files


}   #  end subroutine concat_across_z


sub concat_across_x
{
 use strict;
 use warnings;
 my $numArgs = @_;
 if($numArgs != 5 )
 {
     print "usage: concat_across_x tslice p1 out_dir out_root myid \n";
     last;
 }

  my $tslice=$_[0];
  my $p1=$_[1];
  my $out_dir=$_[2];      # e.g. ../output/slices/3D
  my $out_root=$_[3];     # e.g. slice_
  my $myid=$_[4];         # mpi id of processor running the this script

  my $NCKS='ncks --64bit_offset';
  my $NCPDQ='ncpdq --64bit_offset';
  my $NCRCAT='ncrcat --64bit_offset';

  my($pid,$char_i,$ncfile,$outfile,$cmd);
  my $istep=sprintf("%06d",$tslice);


  #------------------------------------
  #  Loop over pid from 0 to p1-1
  #  swap record dimension from k to i
  #  overwriting the file
  #------------------------------------
  for($pid=0; $pid<$p1; $pid++){
	 $char_i = sprintf("%03d", $pid);
	 $ncfile = $out_dir . "/" . "slice_$istep" . "_" . $char_i . "_" . $myid . ".nc";
	 if($p1 > 1){
	  $cmd="$NCPDQ -O -a idimension,kdimension $ncfile $ncfile >& /dev/null";
	  system($cmd);
	 }
   }

  #------------------------------------
  #  concatenate the p1 files
  #------------------------------------
  $outfile = $out_dir . "/" . $out_root . $istep . ".nc";
  if($p1 > 1){
   $cmd="$NCRCAT -O " . $out_dir . "/" . "slice_$istep" . '_*'  . "_" . $myid . ".nc" . " $outfile";
   #print "$cmd \n";
   system($cmd);
  }
   else{
    $ncfile = $out_dir . "/" . "slice_$istep" . "_" . $char_i . "_" . $myid . ".nc";
    $cmd="mv $ncfile" . " $outfile";
    system($cmd);
  }

  #------------------------------------
  #  clean up the p1 temp files
  #------------------------------------
  $cmd="rm -f " . $out_dir . "/" . "slice_$istep" . '_*'  . "_" . $myid . ".nc";
  system($cmd);

}  # end subroutine concat_across_x $tslice $p1 $out_dir $out_root $myid;






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
 #  ./create_global_snapshot.pl 2 1 4 '../../../output/2D' 'XZplane' '../../../output/slices/2D' 'XZ' 0

 #  NB  2D files are "append", w/ many time slices (starting with 0)
 #  ./create_global_snapshot.pl 0 1 4 '../../../output/2D' 'YZplane' '../../../output/slices/2D' 'YZ' 0

 #  NB  3D files are "new", w/ only one slice, request time slice 0, specify which slice in filename
 #  ./create_global_snapshot.pl 0 1 4 '../../../output/3D' 'XYZ_000030' '../../../output/slices/3D' 'XYZ_30' 0

 #  XY planes entirely contained on single  "j" processor when p1=1
 #  e.g.  output/2D/XYplane_iii-003.nc
 #
 # ./concat_XY.pl tslice p1 jproc data_dir data_root out_dir out_root myid
 #  ./concat_XY.pl 1 2 3 '../../../output/2D' 'XYplane' '../../../output/slices/2D' 'XY' 0
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








