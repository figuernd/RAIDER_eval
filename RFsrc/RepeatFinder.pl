#!/usr/local/bin/perl

# Natalia Volfovsky and Steven L. Salzberg April 2000
#Usage :
#RepeatFinder.pl -f <input file in multiple fasta format> -r <length of minimal exact repeat>  -o <overlap (%)> -g <gap> -e <expect value for BLAST >  -i <identification name> 

# Indentification of the  options

use Getopt::Std;
#use IPC::Open2;

#DEFAULT VALUES
$MUM=50;
$GAP=25;
$EXPECT=0.00000001;
$VERB="yes";
$DIR=".";

# Indentification of the  options
getopt('ifregomvd');

if(!$opt_f or !$opt_i)
{
die("\nPlease enter input file name and identification name :\n          
RepeatFinder.pl -f <input file in multiple fasta format> -i  <identification name>  \n                                                                    
DEFAULT VALUES: minimal exact repeat length 50\n 
                gap for merging             25\n
                expect value for BLAST     0.00000001\n
                verbose                    yes \n
                scrap directory            current directory\n
To change this values - the full command line:
   
RepeatFinder.pl -f <input file in multiple fasta format> -r <length of minimal exact repeat> -g <gap>[-o <overlap>] -e <expect value for BLAST >  -i <identification name> -v <verbose[yes/no]> -d <full path of scrap directory > \n\n")
}
if(!$opt_d){
    $opt_d=$DIR;
}else{
unless (-d "$opt_d"){die("Input error: \n There are no directory $opt_d");}
}
if(-e <*.$opt_i>) {die("Task with  id $opt_i is running or solution exists");}
if(-e <*.$opt_i.>) {die("Task with  id $opt_i is running or solution exists");}
unless(-r "$opt_f") {die("can't read file $opt_f \n");}

if(!$opt_e){
    $opt_e=$EXPECT;
}else{
unless ($opt_e=~ /\d+.?\d*/) {die(" Wrong expect value : $opt_e \n");}
}
if(!$opt_r){
    $opt_r=$MUM;
}elsif ($opt_r=~ /\D+/ ){die(" Wrong minimal repeat length : $opt_r \n"); }

if(!length($opt_g)){
    $opt_g=$GAP;
}
if ($opt_o>100 or $opt_o <0  ) {die("Wrong overlap : $opt_o \n");}
if ( !$opt_m ){ 
$maxl=1000000000;
}
else {
$maxl=$opt_m;
}
if(!length($opt_v)){
    $opt_v=1;
}elsif($opt_v=~/[yes|y]/i){
    $opt_v=1;
}
else{
    $opt_v=0;
}

 
$outfilename1 =$opt_d."/real.seq.$opt_i";
$outfilename2=$opt_d."/seq.coord.$opt_i";
$outfilename3=$opt_d."/all.repeats.$opt_i";
$outfilename4=$opt_d."/all.coord.$opt_i";
$outfilename5=$opt_d."/Qresult.$opt_i";
$outfilename6=$opt_d."/R_map.$opt_i";

###INPUT of the file in multiple fasta format and transform it into onesequence
#file

($mfasta,@seqcoord)=data_input($opt_f);

open( outfile1, ">$outfilename1");
print outfile1 ">all.data\n";
for ($i = 1; $i <= length($mfasta) ; $i++) {
		 
			print outfile1 substr($mfasta,$i-1,1);
		
                       if (($i % 60) == 0) {
			    print outfile1 "\n";
			}
		    }
close(outfile1);
###
open( outfile2, ">$outfilename2"); 
for $i ( 0 .. $#seqcoord ) {
      print outfile2 "\t  @{$seqcoord[$i]} \n";
  }
close(outfile2);

####REPUTER CALL######
@allr=Reputer_call($opt_i,$opt_r,$maxl,$outfilename1,@seqcoord);

open( outfile3, ">$outfilename3"); 
open( outfile4, ">$outfilename4"); 
$last=$allr[0][0];
print outfile4 "\t  $last \n";
for $i ( 0 .. $#allr ) {
      if($last!=$allr[$i][0]){
       print outfile4 "\t  $allr[$i][0] \n";
       $last =$allr[$i][0];
   }
      print outfile3 "\t  @{$allr[$i]} \n";
  }
close(outfile3);
close(outfile4);
###
####CLASS DEFINITION  before BLAST
#
if(length($opt_o)>0){
$command="BBclass_search  -r $opt_r -o $opt_o  $outfilename3 $outfilename4 $outfilename2";
}else{
$command="GBBclass_search  -r $opt_r -g $opt_g  $outfilename3 $outfilename4 $outfilename2";
}
$outmapname=$opt_d."/Repeat.mfasta.$opt_i";

$verb=1;
@Map=ClassFinder($opt_i,$opt_r,$opt_o,$command,$outmapname,$mfasta,$verb,$out_name);
##
####BLAST
@Qresult=BlastBlock($opt_i,$opt_e,$opt_d,@Map);


open( outfile5, ">$outfilename5");
$last1=$Qresult[0][0];
$last2=$Qresult[0][1];
print outfile5 "$Qresult[0][0] $Qresult[0][1] \n";
for $i ( 0 .. $#Qresult ) {
    if($last1!=$Qresult[$i][0] or $last2!=$Qresult[$i][1]){
	print outfile5 "\t  @{$Qresult[$i]} \n";
        $last1=$Qresult[$i][0];
	$last2=$Qresult[$i][1];
    }
}
close(outfile5);

open( outfile6, ">$outfilename6");
for $i ( 0 .. $#Map ) {
print outfile6 " $Map[$i][0] $Map[$i][1] $Map[$i][2] $Map[$i][3]\n";
}
close(outfile6);

# REDISTRIBUTION  after BLAST


$command="ABclass_search $outfilename6 $outfilename5 $outfilename2 AB_map.$opt_i";


$outmapname="AB.mfasta.$opt_i.";
$out_name="AB.map.$opt_i.";


@Map2=ClassFinder($opt_i,$opt_r,$opt_o,$command,$outmapname,$mfasta,$opt_v,$out_name);

CleanDir($opt_i,$opt_d);
###########################
###########################
###########################
sub data_input{
 my($opt_f)=@_;
 my($index,$index1,$prev_asmbl,$line, $asmbl_id,$mfasta,@seqcoord);
 unless (open(SF,$opt_f))  {
     die ("can't open file $opt_f.\n");
 }
 $seq="";
 #$mfasta=">all.data\n";
 $mfasta="";
 $index=0;
 $index1=1;
 $indS=0;
 $prev_asmbl=0;
 while  ($line = <SF>)  {
     
     if ($line =~/>/){
	 $line=~ s/>//g;
	 $asmbl_id=$line;
	 if(!$prev_asmbl) {
	     $prev_asmbl=$line;
	 }
	 if (substr($seq,0,1)=~/[ATGCatgc]/){
	     $seqcoord[$indS][0]=$prev_asmbl;
	     $seqcoord[$indS][1]=$index1;
	     $indS++;
	     $prev_asmbl=$asmbl_id;
	     $seq =~ s/\s+//g;
	     $mfasta.=$seq;
	     for ($i = 1; $i <= length($seq) ; $i++) {
		 $index++;
	     }
	     $index1=$index+1;
	 }
	 
	 $seq="";
     }

     else {
	 $seq.=$line;
     }
     
 }
 
 
 $seqcoord[$indS][0]=$asmbl_id;
 $seqcoord[$indS][1]=$index1;
 $indS++;
 $seq =~ s/\s+//g;
 $mfasta.=$seq;
 for ($i = 1; $i <= length($seq) ; $i++) {
     $index++;
 }
 $asmbl_id=-1;
 $seqcoord[$indS][0]=$asmbl_id;
 $seqcoord[$indS][1]=$index;
 
 return($mfasta,@seqcoord);
}

############################
sub Reputer_call{
    my($opt_i,$opt_r,$maxl,$outfilename1,@seqcoord)=@_;
    my($line,$i,@allrsrt);
    $cmd="reputer -f -p -l $opt_r $outfilename1";
    
    open(res,"$cmd |");
    $i=0;$i1=0;
    while($line=<res>){
	if($line!~/\#/){
	    ($one,$two,$three,$four)=split(/\s/,$line);
            #print "$two, $four, $one \n";
            if($one<= $maxl and $two !=$four and $four-$two+1>$one){
	    @cl=cleanRepeat($two+1,$four+1,$one,@seqcoord);
 	      
             while($#cl != $#buffer1 ){
               for($k1=0;$k1<=$#cl ;$k1++){
		   $buffer1[$k1][0]=$cl[$k1][0];
                   $buffer1[$k1][1]=$cl[$k1][1];
                   $buffer1[$k1][2]=$cl[$k1][2];
	       }
            for ($k=0;$k<=$#buffer1 ;$k++){
		@buffer2=cleanRepeat($buffer1[k][0], $buffer1[k][1],$buffer1[k][2],@seqcoord);
               for($k1=0;$k1<=$#buffer2 ;$k1++){
		   $cl[$k+$k1][0]=$buffer2[$k1][0];
                   $cl[$k+$k1][1]=$buffer2[$k1][1];
                   $cl[$k+$k1][2]=$buffer2[$k1][2];
	       }
           
	    }
          }  
            for($k=0;$k<=$#cl;$k++){
#clean
                if($k){$i=$i+2;
		       $i1=$i1+2;}

		$allr[$i][0]=  $cl[$k][0];
		$allr[$i+1][1]=$cl[$k][0];    
	        $allr[$i][1]=  $cl[$k][1];
		$allr[$i+1][0]=$cl[$k][1];
		$allr[$i][2]=  $cl[$k][2];
		$allr[$i+1][2]=$cl[$k][2];
                $allr[$i][3]=$i1;
		$allr[$i+1][3]=$i1;
	    }    
		
	    #print "$allr[$i][0], $allr[$i][1], $allr[$i][2] \n";
	    $i=$i+2;
	    $i1++;
	}
	#print $line;
	}
    }
    @allrsrt=sort{ &my_sort1 ($a,$b)} @allr;
	return(@allrsrt);
}
##########
sub cleanRepeat{
 my($cor1,$cor2,$length,@seqcoord)=@_;
 my(@check);
 $lo=0;$hi=$#seqcoord ;
 #$cor1++;
 #$cor2++;
 $Curcbac2= BinaryBsearch($lo,$hi,$cor2,@seqcoord);
 $Curcbac1= BinaryBsearch($lo,$hi,$cor1,@seqcoord);
 $Curcbacl2=BinaryBsearch($lo,$hi,$cor2+$length-1,@seqcoord);
 $Curcbacl1=BinaryBsearch($lo,$hi,$cor1+$length-1,@seqcoord);
 
  #print "$cor1 $cor2 $Curcbac2  $Curcbacl2 $Curcbac1 $Curcbacl1 \n";

 $right=0;$left=0; $right1=0;$left1=0;

 if($Curcbac2== $Curcbacl2 && $Curcbac1== $Curcbacl1){
     ##fprintf (fo1,"%10d %10d %10d \n",  cr1,cr2,leng2);
     $check[0][0]=$cor1;
     $check[0][1]=$cor2; 
     $check[0][2]=$length;
   #  ($check[0][0],$check[0][1],$check[0][2])=cleanRepeat($cor1,$cor2,$length,@seqcoord);
    
 }
 elsif($Curcbac2!= $Curcbacl2 && $Curcbac1== $Curcbacl1){
     $Bound22= $seqcoord[$Curcbac2+1][1] ; 
     $diff2= $Bound22-$cor2;
     if ($diff2>=$opt_r){
	 # fprintf (fo1,"%10d %10d %10d \n",  cr1,cr2,diff2);
	 $check[0][0]=$cor1;
	 $check[0][1]=$cor2; 
	 $check[0][2]=$diff2;
#  ($check[0][0],$check[0][1],$check[0][2])=cleanRepeat($cor1,$cor2,$length,@seqcoord);
     }
     else {$left=1;}
     
     $diff22=$cor2+$length-$Bound22-1;
     if ($diff22>=$opt_r){
	 #fprintf (fo1,"%10d %10d %10d \n",  cr1+leng2-diff22,Bound22,diff22);
	 if(!$check[0][0]) {$j=0;}else {$j=1;}
	 $check[$j][0]=$cor1+$length-$diff22-1;
	 $check[$j][1]=$Bound22; 
	 $check[$j][2]=$diff22;
     }
     elsif($left){} #/*fprintf (ferr,"%10d %10d %10d\n",  cr1,cr2,leng2)*/;
 }
 elsif($Curcbac2== $Curcbacl2 && $Curcbac1!=$ Curcbacl1){
     $Bound11= $seqcoord[$Curcbac1+1][1] ;
     $diff1= $Bound11-$cor1;
     if ($diff1>=$opt_r){
	 # fprintf (fo1,"%10d %10d %10d \n",  cr1,cr2,diff1);
         $check[0][0]=$cor1;
	 $check[0][1]=$cor2; 
	 $check[0][2]=$diff1;
     }
     else {$right=1; }
     $diff11=$cor1+$length-$Bound11-1;
     if ($diff11>=$opt_r){
	 #   fprintf (fo1,"%10d %10d %10d \n",  Bound11,cr2+leng2-diff11,diff11);
	 if(!$check[0][0]) {$j=0;}else {$j=1;}
	 $check[$j][0]=$Bound11;
	 $check[$j][1]= $cor2+$length-$diff11-1;
	 $check[$j][2]=$diff11;
     }
	 elsif ($right) {}#/*fprintf (ferr,"%10d %10d %10d\n",  cr1,cr2,leng2)*/
 }
 elsif ($Curcbac2!= $Curcbacl2 and $Curcbac1!= $Curcbacl1){
     
     $Bound11= $seqcoord[$Curcbac1+1][1] ;
     $Bound22= $seqcoord[$Curcbac2+1][1] ; 
     
     $diff1= $Bound11-$cor1;
     $diff2= $Bound22-$cor2;
     
     $gendiff=$diff1<$diff2 ? $diff1:$diff2;
     
     if ($gendiff >=$opt_r){
	 #fprintf (fo1,"%10d %10d %10d \n",  cr1,cr2,gen_diff);
         $check[0][0]=$cor1;
	 $check[0][1]=$cor2; 
	 $check[0][2]=$gendiff;
     }
     else{
	 $left=1;}
     
     $diff11=$cor1+$length-$Bound11 -1;
     $diff22=$cor2+$length-$Bound22-1;
     
     $gendiff = $diff11 < $diff22 ? $diff11 : $diff22;
     
     if ($gendiff >=$opt_r){
	   #fprintf (fo1,"%10d %10d %10d \n",cr1+leng2-gen_diff,cr2+leng2-gen_diff,gen_diff);
         if(!$check[0][0]) {$j=0;}else {$j=1;}
	 $check[$j][0]= $cor1+$length-$gendiff-1;
	 $check[$j][1]= $cor2+$length-$gendiff-1;
	 $check[$j][2]= $gendiff;
     }
     else{
	     $right=1;
	 }
     $gendiff = $diff1 < $diff22 ? $diff1 : $diff22;
     if ($gendiff >=$opt_r && $right && $left){
	   #fprintf (fo1,"%10d %10d %10d \n",cr1,cr2+leng2-gen_diff,gen_diff);
         if(!$check[0][0]) {$j=0;}else {$j=1;}
	 $check[$j][0]= $cor1;
	 $check[$j][1]= $cor2+$length-$gendiff-1;
	 $check[$j][2]= $gendiff;
     }
     else{
	     $right1=1;
	 }
     $gendiff = $diff11 < $diff2 ? $diff11 : $diff2;
     
     if ($gendiff >=$opt_r && $right && $right1 && $left){
	 #fprintf (fo1,"%10d %10d %10d \n",cr1+leng2-gen_diff,cr2,gen_diff);
	 if(!$check[0][0]) {$j=0;}else {$j=1;}
	 $check[$j][0]= $cor1+$length-$gendiff-1;
	 $check[$j][1]= $cor2;
	 $check[$j][2]= $gendiff;
     }
     else{
	     $left1=1;
	 }
     if($right && $left && $right1 && $left1){}
     #/*fprintf (ferr,"%10d %10d %10d\n",  cr1,cr2,leng2)*/;
 }
     
return(@check); 
}
######
sub ClassFinder{

    my($opt_i,$opt_r,$opt_o,$command,$outmapname,$reals,$verb,$out_name)=@_;
    my($i,$max_length,$length,$coord,$class1,@MapArray);

    if($verb){
	open( outmap, ">$outmapname");
    }
    if($out_name){
      open( out, ">$out_name");
    }
    $max_length=length($reals);
    $num=0;
    $cnum=0;
    $prev_class=1;
 
    open(map1,"$command |");

    while($line1=<map1>){
     $line1=~ s/[()]+//g;
    
     ($stam,$class1,$coord,$length, $copies,$a_name, $loc_coord)=split(/\s+/,$line1);
     $check=$coord +$length-1;
     #print " $class1 $coord  $length $copies\n ";
     
     if ($check < $max_length){
	 $check_st=substr($reals,$coord-1,$length);
	 if(!$a_name){
	     print  outmap ">Class $class1 Coord $coord Length $length Copies $copies\n";}
         else{
             if ($prev_class eq $class1){
		 $cnum++;
	     }
             else {
		 $prev_class=$class1;
                 $cnum=1;
	     }
             $loc2=$loc_coord +$length-1;
             $name="C".$class1."R".$cnum;
#########
#OUTPUT
########
             print out "$class1 $name  $a_name $loc_coord $loc2 \n";
#####
	     if($verb){
	     print  outmap ">$name  Assembl $a_name Coord $loc_coord Length $length Copies $copies\n";
	 }
            
	 }
	 ($MapArray[$num][0], $MapArray[$num][1],$MapArray[$num][2],$MapArray[$num][3],$MapArray[$num][4])=($class1,$coord,$length,$copies, $check_st);
         
         $num++;
        if($verb){
	 for ($i = 0; $i < $length ; $i++) {
	     
	     print  outmap substr($check_st,$i,1);
	
	     if ((($i+1) % 60) == 0) {
		 print  outmap "\n";
	     }
	 }
	 if ((($i) % 60) != 0) {
	     print  outmap "\n";
	 }
	 #print "$check_st \n";
     }
     }
     
     #print $line1;
 }
    close(outmap);
    close(map1);
    close(out);
    return(@MapArray);
}
#######
sub BlastBlock{
    my($opt_i,$opt_e, $opt_d,@Map)=@_;
    my($i,$ind,@qresult);
    $Rname1=$opt_d."/REPdb.$opt_i";
    $Rname2=$opt_d."/Repeat.mfasta.$opt_i";
    system("pressdb -o $Rname1 $Rname2");
    
    $queryname =$opt_d."/query.$opt_i";
    
    #open( outquery, ">$queryname");
    
    $cmd3= "blastn $Rname1 $queryname E=$opt_e B=10000 V=10000";
    $class_prev=0;
    $ind=0;
    $maplength=$#Map;
    for($i=0;$i<=$maplength;$i++){
    
	    open( outquery, ">$queryname");
	    print  outquery ">Class $Map[$i][0] Coord $Map[$i][1] Length $Map[$i][2] Copies $Map[$i][3] \n";
            $class1=$Map[$i][0];
	    print  outquery "$Map[$i][4] \n";
	  
	    close (outquery);
	    open (result, "$cmd3 |");
	    $new_line1=0;
	    while ($line1 = <result>) {
		#print $line1;
		if($line1 =~/>/){
                    #print $line1;
		    $line1=~/>Class ([\d\w\S\s]{10}) Coord ([\d\w\S\s]{10}) Length ([\d\w\S\s]{10}) Copies ([\d\w\S\s]{10})/;
                 ($stam,$rclass1,$rcoord,$rlength, $rcopies,$a_name, $loc_coord)=split(/\s+/,$line1); 
		    $new_line1=1;
	
		    
		}
				
		if ($class1 != $rclass1 and  $new_line1 >0){
		    if (($rclass_prev != $rclass1 and $class1 == $class_prev) or ($class1 != $class_prev)){
			#print outfile2 " $class1 $rclass1 \n";
			$qresult[$ind][0]=$class1<$rclass1 ? $class1 :$rclass1;
                        $qresult[$ind][1]=$class1>$rclass1 ? $class1 :$rclass1;

		     #($qresult[$ind][0],$qresult[$ind][1])=($class1,$rclass1);
			#print "$qresult[$ind][0],$qresult[$ind][1] \n";
                        $ind++;
			$new_line1=0;
			$rclass_prev=$rclass1;
			$class_prev=$class1;
		    }
		}   
	    }
	    close (result);
	    
	}
     @qresult=sort{ &my_sort1 ($a,$b)} @qresult;
    return(@qresult);   
}
#######

#######
sub BinaryBsearch{
    my($lo,$hi,$val,@seqcoord)=@_;
    my($i);
    while($lo<=$hi){
	#$i=($lo +$hi)/2;
        $i=(($lo+$hi)-(($lo+$hi)%2))/2;
#	print "$i $val  $seqcoord[$i][1] $seqcoord[$i+1][1]\n";
	if($seqcoord[$i][1]<=$val and $val<$seqcoord[$i+1][1]){
	   return($i);
       }
       elsif ($val<$seqcoord[$i][1]){
	   $hi=$i-1;
       }
	else{
	    $lo=$i+1;
	}
       
    }
return(-1);
}

########################

sub my_sort1 
  { 
  my $as=shift; 
  my $bs=shift; 
 
  return 1 if($as->[0] > $bs->[0]); 
  return -1 if($as->[0] < $bs->[0]); 
  return 1 if($as->[2] < $bs->[2]); 
  return -1 if($as->[2] > $bs->[2]);
  return 1 if($as->[1] > $bs->[1]); 
  return -1 if($as->[1] <$bs->[1]); 
  return 0; 
  } 


#################
sub PrintFile{
    my($FileName,@Array)=@_;

open( out, ">$FileName");

for $i ( 0 .. $#Array ) {
      print out "\t  @{$Array[$i]} \n";
  }
    close(out);
}
#####################
sub CleanDir{
    my($opt_i,$opt_d)=@_;
   
    
    @deletefiles=(<$opt_d/*.$opt_i>,<$opt_d/REPdb.$opt_i.*>);
    foreach(@deletefiles){
	#print "$_ \n";
	if($_ !~ /^AB/){
	 unlink $_;
     }
    }
}



