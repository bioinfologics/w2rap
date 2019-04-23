#!/usr/bin/perl -w

# Script:  nextclip_index_reference.pl
# Purpose: Produce simple index file of a reference
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

use warnings;
use strict;
use Getopt::Long;

my $in_filename=$ARGV[0];
my $out_filename=$in_filename.".nextclip";

die "Error: You must specify an input filename\n" if not defined $in_filename;

print "Indexing $in_filename to $out_filename\n";

open(INPUTFILE, $in_filename) or die "Can't open $in_filename\n";
open(OUTPUTFILE, ">".$out_filename) or die "Can't open $out_filename\n";

my $id = "";
my $contig_length = 0;
while(<INPUTFILE>) {
    chomp(my $line = $_);
    
    if ($line =~ /^>(\S+)/) {
        if ($contig_length > 0) {
            print OUTPUTFILE $id, "\t", $contig_length, "\n";
        }
        
        $contig_length = 0;
        $id = $1;
    } else {
        $contig_length += length($line);
    }
}

if ($contig_length > 0) {
    print OUTPUTFILE $id, "\t", $contig_length, "\n";
}

close(OUTPUTFILE);
close(INPUTFILE);

print "DONE\n";