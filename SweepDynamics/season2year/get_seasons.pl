#!/bin/perl
use strict;
use warnings;

# set filenames
my $aafilename = $ARGV[0];

# open fasta file with amino acid sequences
open (my $aafile, $aafilename) or die "Could not open file '$aafilename' $!";

my $date = "";

my $minyear = 10000;
my $maxyear = 0;

# get min and max year
while( my $line = <$aafile>)  {
        if ($line =~ /^>/) {

                #extract date
                $line =~ m/\s([12][0-9]{3}[\/|-][01][0-9][\/|-][0123][0-9])\s/;
                $date = $1;

                if (not defined $date) {
                        print "Date not found in amino acid sample $line! Provide date in format YYYY/MM/DD or YYYY-MM-DD\n";
                        next;
                } 

                #year from date
                my @splitdate = split /[\/|-]/, $date;
                my $year = $splitdate[0];

		if ( $year >= $maxyear ) {
			$maxyear = $year;
		}
		if ( $year <= $minyear ) {
			$minyear = $year;
		}
        } 
}

my $minmonth = 12;
seek $aafile, 0, 0;

# get minimum month in minimum year
while( my $line = <$aafile>)  {
        if ($line =~ /^>/) {

                #extract date
                $line =~ m/\s([12][0-9]{3}[\/|-][01][0-9][\/|-][0123][0-9])\s/;
                $date = $1;

                if (not defined $date) {
                        print "Date not found in amino acid sample $line! Provide date in format YYYY/MM/DD or YYYY-MM-DD\n";
                        next;
                } 

                #get month and year from date
                my @splitdate = split /[\/|-]/, $date;
                my $year = $splitdate[0];
                my $month = $splitdate[1];
		
		if ( $year == $minyear ) {
			if ( $month < $minmonth ) {
				$minmonth = $month;
			}
		}
        } 
}

my $startseason = "";
# get starting season
if ( $minmonth > 0 && $minmonth <=3 ) {
	$startseason = "N";
} elsif ( $minmonth > 3 && $minmonth < 10 ) {
	$startseason = "S";
} elsif ( $minmonth >= 10 && $minmonth <=12) {
	$startseason = "N";
	$minyear = $minyear + 1;
}

open(my $newfile, '>', "minmaxseason.txt");
print $newfile "$minyear\t$maxyear\t$startseason";
