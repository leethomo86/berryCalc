#! /usr/bin/perl

use strict;
use warnings;

#use lib "/home/lee/PROGRAMS/PERL/gregap/";
#use ESPT::Glog 0.07;

my %result = ();

my $dir = $ARGV[0];
die "Incorrect number of arguments passed to noci_analysis: ".$#ARGV+1 ."\n" if ($#ARGV ne 0);
my @files = <$dir/*.out>;
foreach my $file (@files){
        open(INPUT,$file) or die "Can't open file $file\n";

        my $finish = '';
        my @dist = ();
        my $freqs = '';
        my $Energy = '';

        while(<INPUT>){
		if (/^JOB TERMINATION SUCCESS/){
                        $finish = 1;
                }
                if (/^\s+NOCI Energies/../^\s+Eigenvectors/){
                        chomp($_);
                        unless (/^\s+NOCI Energies/ or /^\s+Eigenvectors/){
                                my @arr = split(/\s+/,$_);
                                $Energy = $Energy . "$arr[2]\t";
                        }
                }
                if (/^\s+.*\s+scan. NOCI, X =\s+(-?\d+\.?\d*), Y =\s+(-?\d+\.?\d*)/){
                        $dist[0] = sprintf("%.5f",$1);
			$dist[1] = sprintf("%.5f",$2);
                }
        }

        if ($finish ne 1){
                print "$file did not finish\n";
		$Energy = '';
        }

        $result{"$dist[0]"}{"$dist[1]"} = "$Energy";

        close INPUT;
}

open(OUTPUT,">$dir.dat") or die "Can't open $dir.dat\n";

#foreach my $name (sort keys %result) {
#    printf OUTPUT "%4.2f\t%28s\n", $name, $result{$name};
#}
foreach my $name (sort keys %result) {
	foreach my $name2 (sort keys %{ $result{$name} }) {
		my $string = "$name\t$name2\t$result{$name}{$name2}";
		my @fields = split(/\s+/,$string);
		my @sfields = sort { $a <=> $b } @fields[2..$#fields];
		$string = join("\t",@sfields);
        	print OUTPUT "$fields[0]\t$fields[1]\t$string\n";
	}
	print OUTPUT "\n";
}

#&plotopt('plot.dat', 'Li-H distance', 'Energy (Har)', 'ENERGY', '', '', 'lines', '2', '');

    #######################################################################
    # SUBROUTINE PLOTOPT                                                  #
    # CREATED 10/05/12                                                    #
    # LAST MODIFIED 10/05/12                                              #
    # LEE THOMPSON                                                        #
    # This subroutine takes data from plot.dat, constructs a gnuplot file #
    # plot.plt and displays the gnuplot graphic. Series and columns input #
    # use "|" as a delimiter. First column is plotted as x data and each  #
    # subsequent column corresponds to y data for each series.            #
    # TODO: Extend by inserting option for surface plot.                  #
    #######################################################################

sub plotopt {
    my($input, $xlabel, $ylabel, $title, $xrange, $yrange, $style, $column, $series) = @_;

    my $command = '';

    open FHDL, ">plot.plt" or die "Cannot open plot input\n";

    my @names = split(/\|/, $series);
    my @columns = split(/\|/, $column);
    if (scalar @columns > 1) {
        $command = sprintf("plot \"%s\" using 1:%s title \"%s\"\n", $input,$columns[0], $names[0]);
        for (my $i = 1; $i <= (scalar @columns - 1); $i++) {
            $command = sprintf("%sreplot \"%s\" using 1:%s title \"%s\"\n",$command, $input, $columns[$i], $names[$i]);
        }
    }
    elsif (scalar @names == $column) {
        $command = sprintf("plot \"%s\" using 1 title \"%s\"\n", $input,$names[0]);
	for (my $i = 1; $i <= $column; $i++) {
            $command = sprintf("%sreplot \"%s\" using %s title \"%s\"\n",$command, $input, $i, $names[$i]);
        }
    }
    else {
        $command= sprintf("plot \"%s\" using 1:%s title \"%s\"\n", $input,$column, $series);
    }

    printf FHDL "set xlabel \"%s\"\n" . "set ylabel \"%s\"\n" . "set xrange [%s]\n" . "set yrange [%s]\n" . "set grid\n" .  "set style data %s\n" . "set title \"%s\"\n" . '%s' . 'pause -1', $xlabel, $ylabel, $xrange, $yrange, $style, $title, $command;
    close FHDL;

    system('gnuplot plot.plt');
}


