#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Find;
use File::Path;
use Cwd;
#use diagnostics;

 ########################### WELCOME  ##############################
#                                                                   #
#                          sgRNAcas9                                # 
# ---a tool for fast designing CRISPR sgRNA with high specificity   #
#                                                                   #
# AUTHOR  : Xie Shengsong                                           #
# Email   : ssxieinfo\@gmail.com                                    #
# Homepage: www.biootools.com                                       #
#           BiooTools (Biological online tools)                     #  
# Shanghai Institute of Biochemistry and Cell Biology (SIBCB)       #
# School of Veterinary Medicine, Huazhong Agricultural University   #
# Version: sgRNAcas9_3.0.4                                          #
# Begin       : 2013.12.9                                           #
# LAST REVISED: 2015.2.5                                            #
 ###################################################################

my $oldtime = time();

my ($Inputfile_Fasta, $truncat, $GC_l, $GC_m, $Genome, $Option, $Type, $Seqmap_vesion, $Num_mismatch, $offset_s, $offset_e, $path);

GetOptions( "i=s" => \$Inputfile_Fasta,        #Input file
            "x=i" => \$truncat,                #Length of sgRNA[20]
            "l=i" => \$GC_l,                   #The minimum value of GC content [20]
			"m=i" => \$GC_m,                   #The maximum value of GC content [80] 
	        "g=s" => \$Genome,                 #The reference genome sequence
			"o=s" => \$Option,                 #Searching CRISPR target sites using DNA strands based option(s/a/b)
			"t=s" => \$Type,                   #Type of gRNA searching mode(s/p) 
			"v=s" => \$Seqmap_vesion,          #Operation system [w, for windows; l, for linux-64; u, for linux-32;  m, for MacOSX-64; a, for MacOSX-32]
			"n=i" => \$Num_mismatch,           #Maximum number of mismatches [5]
			"s=i" => \$offset_s,               #The minimum value of sgRNA offset [-2]
            "e=i" => \$offset_e,               #The maximum value of sgRNA offset [32]
            "p=s" => \$path,                   #Output path
          );

#default
$truncat ||= "20";  
$GC_l ||= "20";                                #20 % < GC% < 80 %
$GC_m ||= "80";
$Seqmap_vesion ||= "l";                        #linux
$Num_mismatch ||="5";                          #number of mismatches: 5
$offset_s ||="-3";                             #sgRNA offset: -2 to 32 bp
$offset_e ||="33";
my $dir_default = getcwd;                      #default output
$path ||= $dir_default;

my $dir =$path;

mkdir("$dir/sgRNAcas9.report_$truncat.$Option",0755)||die "Can't create directory: Directory exists at $dir. Please delete, move or rename the exist directory before you run this program.$!" ;

open  (LOG, ">>$dir/sgRNAcas9.report_$truncat.$Option/sgRNAcas9.Log.txt") || die "Can't open sgRNAcas9.Log.txt for writing!" ."\n";
print  LOG "################################# Log ###########################################".        "\n\n";
#print "Writing Log information.                                                                      " ."\n";
print  LOG "#                              sgRNAcas9                                                  " ."\n";
print  LOG "#     ---a tool for fast designing CRISPR sgRNA with high specificity                     " ."\n";          
print  LOG "#                                                                                         " ."\n";
print  LOG "#       contact:  Xie Shengsong, Email: ssxieinfo\@gmail.com                                .\n\n";
######


#mkdir("$dir/sgRNAcas9.report_$truncat.$Option/A.Final_report",0755)||die "can't create directory: $!" ;

print "\n\tWelcome to sgRNAcas9\n";
print "\t---a tool for fast designing CRISPR sgRNA with high specificity\n";
print "\t---------------------------------------------------------\n";
print "Version   : 3.0"."\n";
print "Copyright : Free software"."\n";
print "Author    : Shengsong Xie"."\n";
print "Email     : ssxieinfo\@gmail.com"."\n";
print "Homepage  : www.biootools.com"."\n";

my $local_time;
$local_time = localtime();
print "Today     : $local_time\n\n";
print  LOG "# Time, begin at $local_time."."\n";
print  LOG "# Usage: perl $0 -i $Inputfile_Fasta -x $truncat -l $GC_l -m $GC_m -g $Genome -o $Option -t $Type -v $Seqmap_vesion -n $Num_mismatch -s $offset_s -e $offset_e -p $path 2>log.txt\n\n";

################################### format seq ###################################
print "Start sgRNAcas9 program........\n";
print "Step1: Format target sequences.\n";
print  LOG "# Start sgRNAcas9 program........\n";
print  LOG "# Step1: Format target sequences.\n";

open (Inseq, $Inputfile_Fasta) || die "Can't open $Inputfile_Fasta for reading!\n";
open (FASTA, ">$dir/sgRNAcas9.report_$truncat.$Option/TargetSeq.fa") || die "Can't open TargetSeq.fa for writing!\n";

my $TmpTit="";
my $TmpSeq="";
my $TmpTit_S="";
my $TmpTit_A="";
my $TmpSeq_Acomp;
my $TmpSeq_ARevcomp;

if ($Option eq "b") {              #B=both:sense and anti-sense strand

    while (<Inseq>) {

	  chomp $_; 
	 
      if (/^>(\S+)/){
	   
	    if ($TmpTit && $TmpSeq) {
	
			 $TmpTit_S=$TmpTit."_S";
		     $TmpTit_A=$TmpTit."_A";
		
			 $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
		     $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
		     $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
				
		     print FASTA ">$TmpTit_S\n$TmpSeq\n";	
	         print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";
		    }
	
	         $TmpTit=$1;	
	         $TmpSeq="";
	
	   }elsif(/(\w+)/){
	
	         $TmpSeq.=$1;
		
	  }
	 
	}

       if ($TmpTit && $TmpSeq) {

	        my $TmpTit_S=$TmpTit."_S";
            my $TmpTit_A=$TmpTit."_A";

	        my $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
            $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
            my $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
     
            print FASTA ">$TmpTit_S\n$TmpSeq\n";
            print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n"; 
	
  }
}elsif ($Option eq "a") {           #Non-template: A: anti-sense strand

	while (<Inseq>) {

		chomp $_;

	if (/^>(\S+)/){
	
	    if ($TmpTit && $TmpSeq) {
		  	  
		     $TmpTit_A=$TmpTit."_A";
		
			 $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
		     $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
		     $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
				
		     print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";	
		  }
	
	         $TmpTit=$1;	
	         $TmpSeq="";
	
	  }elsif(/(\w+)/){
	
	         $TmpSeq.=$1;
		
	  }
	}

       if ($TmpTit && $TmpSeq) {

            my $TmpTit_A=$TmpTit."_A";

            my $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
            $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/; 
            my $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
   
            print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";
   
  }
}elsif ($Option eq "s" || $Option eq "") {     #Template: S: sense strand

	while (<Inseq>) {

		chomp $_;

	if (/^>(\S+)/){
	
	   if ($TmpTit && $TmpSeq) {
		  	  
		    $TmpTit_S=$TmpTit."_S";
				
		    print FASTA ">$TmpTit_S\n$TmpSeq\n";	
		  
	    }
	
	        $TmpTit=$1;	
	        $TmpSeq="";
	
	  }elsif(/(\w+)/){
	
	    	$TmpSeq.=$1;
		
	  } 
	}

      if ($TmpTit && $TmpSeq) {
	
		   my $TmpTit_S=$TmpTit."_S";

           print FASTA ">$TmpTit_S\n$TmpSeq\n";
    
  }
}
close(Inseq);
close(FASTA);

########################### Find CRISPR targets-single #################################
print "Step2: Find CRISPR targets.\n";
print  LOG "# Step2: Find CRISPR targets.\n";

open(IntPut, "$dir/sgRNAcas9.report_$truncat.$Option/TargetSeq.fa") || die "Can't open TargetSeq.fa for reading!\n";

open(OutPut, ">$dir/sgRNAcas9.report_$truncat.$Option/report_protospacer_single.txt") || die "Can't open report_protospacer_single.txt for writing!\n";
open(OutPut1, ">$dir/sgRNAcas9.report_$truncat.$Option/CRISPR.targets_single.fa") || die "Can't open CRISPR.targets_single.fa for writing!\n";
#open(OutPut2, ">$dir/sgRNAcas9.report_$truncat.$Option/full_length_sgRNA.txt") || die "Can't open full_length_sgRNA.txt for writing!\n";
open(OutPut3, ">$dir/sgRNAcas9.report_$truncat.$Option/CRISPR.targets_S.txt") || die "Can't open CRISPR.targets_single.fa for writing!\n";
open(OutPut4, ">$dir/sgRNAcas9.report_$truncat.$Option/CRISPR.targets_A.txt") || die "Can't open CRISPR.targets_single.fa for writing!\n";

#print OutPut "Candidate sgRNA, Pattern: GGX18NGG, GX19NGG, X20NGG, GC% >=$GC %\n\n";
print OutPut "sgRID\t"."Start\t"."End\t"."CRISPR_target_sequence(5'-3')\t"."Length(nt)\t"."GC%\n";

my $ID=""; 
my $seq="";
my $dCas9handle = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCG";
my $terminator = "UUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU";
my $probe_id="";

while(<IntPut>) {
	chomp $_;

	if (/^>(\S+)/){
		analysis(); $ID=$1;

  }else{
		$seq=$_;
		
  }
}
analysis();
close OutPut;
close OutPut1;
close OutPut3;
close OutPut4;
################################### subroutine ###################################
sub analysis {  
	if($ID eq "" || $seq eq "") {
	  return();
	}
	#my $seq =~ s/\n//g;
	my $SEQRNA  = uc($seq);
	my $SEQDNA  = $SEQRNA;
	$SEQDNA  =~ tr/U/T/;

	my $len = length($seq);

	my $idx_s;
	my $i = 0;

	for($idx_s=-1; $idx_s<length($seq); $idx_s++) {

		my $lentruncat = $truncat+3;
		my $lentruncat2 = $truncat+1;

		if (substr($SEQDNA, $idx_s+1, $lentruncat) =~ /[ATCG]{$lentruncat2}AG/) {    #X?-NAG  
			my $sgRidx_s = $idx_s+2;                                  #For sense strand
			my $sgR = $&;
			#print $sgR."\n";
			my $SgR1 = substr($sgR,0,-3); 
			#print $SgR1."\n";
			my $sgRAT = &AT($sgR);
			#my $sgRGC = &GC($sgR);
			my $sgRGC = &GC($SgR1);

			#my $SgRsd = substr($sgR,8,12);                           #5-X12-NGG-3
			#my $SgRseed ="$SgRsd"."NGG";
			my $SgRNt = substr($sgR,0,$truncat);                      #truncat
			my $sgRlen = length($sgR);
			my $sgRidx_e = $sgRidx_s+$sgRlen-1;
			$i++;
			
			my $sgRidx_A_s =$len-$sgRidx_e+1;                         #For anti-sense strand
			my $sgRidx_A_e =$len-$sgRidx_s+1;
			
			my $sgrna = $SgRNt;
			$sgrna =~ tr/T/U/;
			my $sgRNA = $sgrna.$dCas9handle.$terminator;

				if (!($sgR =~ /T{4,18}/g)){
					#if (!($sgR =~ /A{5,21}|T{4,21}|C{6,21}|G{6,21}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|GC{6,10}|(AAT){5,7}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/g) ){  

					if ($ID=~ /._S$/ and ($sgRGC >= $GC_l and $sgRGC < $GC_m)) {    #For sense strand
           
						print OutPut "$ID\_$i\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						print OutPut3 "$ID\_$i\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						
						print OutPut1 ">$ID\_$i\n"."$sgR\n";
						#print OutPut2 ">$ID\_$i\n"."$sgRNA"."\n";

					}elsif ($ID=~ /._A$/ and ($sgRGC >= $GC_l and $sgRGC < $GC_m) ) {   #For anti-sense strand
           
						print OutPut "$ID\_$i\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						print OutPut4 "$ID\_$i\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						
						print OutPut1 ">$ID\_$i\n"."$sgR\n";
						#print OutPut2 ">$ID\_$i\n"."$sgRNA"."\n";

    			}
				}
		}
	}
}

#######subroutine to calculate GC% content
sub GC { 

    my $seq2 = shift @_; 

    $seq2 = uc($seq2);

    my $seq2DNA;
    $seq2DNA = $seq2;
    $seq2DNA =~ tr/U/T/;

    my $seq2length = length($seq2DNA);

    my $count = 0;
    for (my $i = 0; $i < $seq2length; $i++) {
			my $sub = substr($seq2DNA,$i,1);
			if ($sub =~ /G|C/i) {
	    	$count++;
			}
    }
    my $gc = sprintf("%.1f",$count * 100 /$seq2length);
    return $gc;
}


#######subroutine to calculate AT% content
sub AT { 

    my $seq2 = shift @_; 

    $seq2 = uc($seq2);

    my $seq2DNA;
    $seq2DNA = $seq2;
    $seq2DNA =~ tr/U/T/;

    my $seq2length = length($seq2DNA);

    my $count = 0;
    for (my $i = 0; $i < $seq2length; $i++) {
			my $sub = substr($seq2DNA,$i,1);
			if ($sub =~ /A|T/i) {
	    	$count++;
			}
    }
    my $gc = sprintf("%.1f",$count * 100 /$seq2length);
    return $gc;
}

close IntPut;
close OutPut;
#exit;
