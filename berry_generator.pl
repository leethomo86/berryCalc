#! /usr/bin/perl

use strict;
use warnings;

use Parallel::ForkManager;
use Cwd;
use File::Copy;
use Math::Trig;

my $rad = 0.01;
my $steps = 5;

my $symmetry = 'time-reversal';
my $method = 'rhf';
my $inputfile = 'input.txt';
open(INPUT,'<',$inputfile) or die $!;
chomp(my $numfiles = <INPUT>);
my @chkfile = ();
while(<INPUT>){
  chomp($_);
  $_ =~ s/^\s+|\s+$//g;
  push(@chkfile,$_);
}
close(INPUT);
die "Error: Number of input files not equal to number specified: $!" unless ($#chkfile+1==$numfiles);  

my $g16root = $ENV{g16root};
my $i = 0;
my $MAX_PROCESSES = 4;
my $nproc=1;
my $mem="1GB";
my $sub = 1;
my $scfsub = 1;
my $proj_state = 0;
my $berry_state = 0;
my $direct = '';

foreach my $temp (@ARGV){
      chomp($temp);
      if($temp =~ /^-nproc=(\d+)$/i){
        $nproc=$1;
      }elsif($temp =~ /^-mem=(\d+[a-zA-Z][a-zA-Z])$/i){
        $mem=$1;
      }elsif($temp =~ /^-nosub/i){
        $sub = 0;
      }elsif($temp =~ /^-noscf/i){
        $scfsub = 0;
      }elsif($temp =~ /^-projstate=(\d+)$/i){
        $proj_state = $1;
      }elsif($temp =~ /^-nocistate=(\d+)$/i){
        $berry_state = $1;
      }elsif($temp =~ /^-direct/i){
        $direct = '--direct';
      }else{
	die("Unrecognized option in Berry phase generator,\n")
      }
}

my @coords = (['O',  0.000000,  0.000000,  0.000000], 
              ['O',  0.000000,  1.200000,  1.150000],
              ['O',  0.000000, -1.200000,  1.150000]);

my @dc = ([ 0.0000000000,  0.0000000000,  0.0000000000], 
          [ 0.0000000000,  1.0000000000,  0.0000000000],
          [ 0.0000000000, -1.0000000000,  0.0000000000]);
  
my @gd = ([ 0.0000000000,  0.0000000000,  0.0000000000],  
          [ 0.0000000000,  0.0000000000,  1.0000000000],
          [ 0.0000000000,  0.0000000000,  1.0000000000]);

my $curdir = getcwd;
die("Projection requested on state that does not exist") if ($proj_state lt -1 or $proj_state ge $numfiles); 
die("NOCI requested on state that does not exist") if ($berry_state lt -1 or $berry_state ge $numfiles); 

if ($scfsub eq 1){
	for (my $state=0; $state<=$numfiles-1; $state++){
		mkdir("STATE_$state") or $!{EEXIST} or die("Can't create directory \"STATE_$state\": $!\n");
		chdir("STATE_$state");
		for (my $theta=0.0; $theta<2*pi; $theta+=2*pi/$steps){
			my $i = $rad*cos($theta);
			my $j = $rad*sin($theta);
		        my $gjf = "ozone_sol"."$state"."_angle_".sprintf("%.3f",$theta);
		        open  (GJF,">$gjf.com") or die "Could not create $gjf.\n";
		        print GJF "%oldchk=$curdir/$chkfile[$state]\n";
		        print GJF "%chk=$curdir/STATE_$state/$gjf.chk\n";
		        print GJF "#P $method chkbas guess=read output=matrixelement\n";
		        if($state eq 0 or $state eq $proj_state){
				print GJF "# IOp(5/194=4,5/10=500) scf=(conven,conver=6,novaracc,maxcycles=250) int=raf3 nosymm\n\n";
			}else{
				print GJF "# IOp(5/194=4,5/10=500) scf=(conven,conver=6,novaracc,maxcycles=250) nosymm\n\n";
			}
		        print GJF "Ozone scan. State = ".$state.", X = ".sprintf("%.5f",$i).", Y = ".sprintf("%.5f",$j)."\n\n";
		        print GJF "0 1\n";
			for(my $line=0; $line<=$#coords; $line++){
			  	printf GJF "%s %.6f %.6f %.6f\n",$coords[$line][0], 
			  		$coords[$line][1]+$i*$dc[$line][0]+$j*$gd[$line][0],
			  		$coords[$line][2]+$i*$dc[$line][1]+$j*$gd[$line][1],
			  		$coords[$line][3]+$i*$dc[$line][2]+$j*$gd[$line][2];
			}
		        print GJF "\n$curdir/STATE_$state/$gjf.mat\n\n";
		        close GJF;
		}
		chdir("..");
	}
}

#build NOCI directory with correct files pointing to SCF directories
if ($proj_state ne -1){
	mkdir("PROJ") or $!{EEXIST} or die("Can't create directory \"NOCI\": $!\n");
	chdir("PROJ");
	for (my $theta=0.0; $theta<2*pi; $theta+=2*pi/$steps){
		my $i = $rad*cos($theta);
		my $j = $rad*sin($theta);
	        my $gjf = "ozone_PROJ_angle_".sprintf("%.3f",$theta);
	        open  (GJF,">$gjf.input") or die "Could not create $gjf.\n";
		print GJF "$curdir/STATE_$proj_state/ozone_sol"."$proj_state"."_angle_".sprintf("%.3f",$theta).".mat\n";
	        close GJF;
		open (OUT,">$gjf.out") or die "Could not open $gjf.out.\n";
		print OUT " Ozone scan. PROJ, X = ".sprintf("%.5f",$i).", Y = ".sprintf("%.5f",$j)."\n\n";
	        close OUT;
	}
	chdir("..");
}

#build NOCI directory with correct files pointing to SCF directories
if ($berry_state ne -1){
	mkdir("NOCI") or $!{EEXIST} or die("Can't create directory \"NOCI\": $!\n");
	chdir("NOCI");
	for (my $theta=0.0; $theta<2*pi; $theta+=2*pi/$steps){
		my $i = $rad*cos($theta);
		my $j = $rad*sin($theta);
	        my $gjf = "ozone_NOCI_angle_".sprintf("%.3f",$theta);
	        open  (GJF,">$gjf.input") or die "Could not create $gjf.\n";
		print GJF "$numfiles\n";
		for (my $state=0; $state<=$numfiles-1; $state++){
			print GJF "$curdir/STATE_$state/ozone_sol"."$state"."_angle_".sprintf("%.3f",$theta).".mat\n";
		}
	        close GJF;
		open (OUT,">$gjf.out") or die "Could not open $gjf.out.\n";
		print OUT " Ozone scan. NOCI, X = ".sprintf("%.5f",$i).", Y = ".sprintf("%.5f",$j)."\n\n";
	        close OUT;
	}
	chdir("..");
}

#go into each directory and run each input before running NOCI
if ($sub eq 1){
	if ($scfsub eq 1){
		for (my $state=0; $state<=$numfiles-1; $state++){
			chdir("STATE_$state");
        		my @files = <*.com>;

        		my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

        		foreach my $file (@files){
        		        $pm->start and next;

        		        my $log_file = $file;
        		        unless($log_file =~ s/\.(?:gjf|com)/.log/){
        		                $log_file .= ".log";
        		        }
        		        my $sys_cmd = "$g16root/g16/g16 -m=$mem -p=$nproc < $file > $log_file";
        		        system("$sys_cmd");
        		        print "Job $file completed.\n";

        		        $pm->finish;
        		}
        		$pm->wait_all_children;
			chdir("..");
		}
	}

	if ($proj_state ne -1){
		chdir("PROJ");
        	my @files = <*.input>;

        	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

        	foreach my $file (@files){
        	        $pm->start and next;

        		my $out_file = $file;
        		unless($out_file =~ s/\.(?:input)/.out/){
        		        $out_file .= ".out";
        		}
			open (INPUT,$file) or die "Can't open input file $file\n";
			chomp(my $projfile = <INPUT>);
			close INPUT;
        	        my $sys_cmd = "~/Documents/Wheeler-HIll-PUHF/WH_PUHF.exe -f $projfile $direct --scan-output --symmetry $symmetry >> $out_file";
        	        system("$sys_cmd");
			open (OUT,">>$out_file") or die "Could not open $out_file.\n";
        	        if ($? eq 0){
        	                print OUT "JOB TERMINATION SUCCESS\n";
        	        }else{
        	                print OUT "JOB TERMINATED WITH ERROR\n";
        	        }
        	        close OUT;
        	        print "Job $file completed.\n";
        	        $pm->finish;
        	}
        	$pm->wait_all_children;
		chdir("..");
	}

	if ($berry_state ne -1){
		chdir("NOCI");
        	my @files = <*.input>;
        	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

        	foreach my $file (@files){
			$pm->start and next;
        		my $out_file = $file;
        		unless($out_file =~ s/\.(?:input)/.out/){
        		        $out_file .= ".out";
        		}
        		my $sys_cmd = "~/Documents/ResHF/ResHF.exe -f $file $direct --scan-output $berry_state >> $out_file";
			system("$sys_cmd");
			open (OUT,">>$out_file") or die "Could not open $out_file.\n";
        	        if ($? eq 0){
        	                print OUT "JOB TERMINATION SUCCESS\n";
        	        }else{
        	                print OUT "JOB TERMINATED WITH ERROR\n";
        	        }
        	        close OUT;
			print "Job $file completed.\n";
        		$pm->finish;
		}
        	$pm->wait_all_children;
		chdir("..");
	}
}
