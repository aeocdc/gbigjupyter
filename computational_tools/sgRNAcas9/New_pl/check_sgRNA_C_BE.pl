#!/usr/bin/perl
#$PAM：所运行的PAM区的名称

use Getopt::Long;
# $fi="/home/nixm/run/shm/Project/0-sgRNAcas9/sgRNAcas9_V3.0_GUI/report/sgRNA_rm_repeat.txt";
# $fo="/home/nixm/run/shm/Project/0-sgRNAcas9/sgRNAcas9_V3.0_GUI/report/sgRNA_result-editable.txt";


my ($Inputfile_Fasta, $Outputfile_Fasta);

GetOptions( "i=s" => \$Inputfile_Fasta,        #Input file
            "o=s" => \$Outputfile_Fasta,                #Output file
          );
open FI, "$Inputfile_Fasta" or die "can't open file in\n";
open FO, ">$Outputfile_Fasta" or die "can't open file out\n";

while(<FI>){
	chomp;
	my @str=split/\t/;
	my @length=$str[0];
	my @length=split/_/;
	print FO "$_";
	if ($str[0]=~ m/(_S_)/g){
		my $seq=substr($str[3],3,7);
		while ($seq=~ m/(CAA)/g) {
			my $end=pos($seq);
			my $ss=($str[1]+$end-4)%3;
			if($ss==0){
				my $ll=$end-2;
				my $per=($str[1]+$end-3)*100/$length[0];
				$per=sprintf "%.2f",$per;
				print FO "\tCAA\t1\t$ll\t$per";
			}
			else{
				print FO "\tCAA\t0\t-\t-";
			}
		}
		while ($seq=~ m/(CGA)/g) {
			my $end=pos($seq);
			my $ss=($str[1]+$end-4)%3;
			if($ss==0){
				my $ll=$end-2;
				my $per=($str[1]+$end-3)*100/$length[0];
				$per=sprintf "%.2f",$per;
				print FO "\tCGA\t1\t$ll\t$per";
			}
			else{
				print FO "\tCGA\t0\t-\t-";
			}
		}
		while ($seq=~ m/(CAG)/g) {
			my $end=pos($seq);
			my $ss=($str[1]+$end-4)%3;
			if($ss==0){
				my $ll=$end-2;
				my $per=($str[1]+$end-3)*100/$length[0];
				$per=sprintf "%.2f",$per;
				print FO "\tCAG\t1\t$ll\t$per";
			}
			else{
				print FO "\tCAG\t0\t-\t-";
			}
		}
		print FO "\n";
	}
	else{
		my $seq=substr($str[3],3,7);
		while ($seq=~ m/(CCA)/g) {
			my $end=pos($seq);
			my $ss=($str[2]-$end)%3;
			if($ss==0){
				my $ll=$end-2;
				my $per=($str[2]-$end+3)*100/$length[0];
				$per=sprintf "%.2f",$per;
				print FO "\tCCA\t1\t$ll\t$per";
			}
			else{
				print FO "\tCCA\t0\t-\t-";
			}
		}
			print FO "\n";
	}
		
}
close FI;
close FO;