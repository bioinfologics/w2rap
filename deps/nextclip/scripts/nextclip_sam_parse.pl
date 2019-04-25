#!/usr/bin/perl -w

# Script:  nextclip_sam_parse.pl
# Purpose: Parse SAM files for NextClip
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/


use warnings;
use strict;
use Getopt::Long;

my $read_one;
my $read_two;
my $out_prefix;

my $unmapped = 0;
my $n_single_mapped = 0;
my $n_map_differently = 0;
my $n_poor_mapping = 0;
my $n_mp = 0;
my $n_pe = 0;
my $n_tandem = 0;
my $n_mp_out_of_range = 0;
my $n_pe_out_of_range = 0;
my $n_tandem_out_of_range = 0;
my $n_dupes = 0;
my $sum_mp_one = 0;
my $sum_mp_two = 0;
my $sum_pe_one = 0;
my $sum_pe_two = 0;
my $sum_tandem_one = 0;
my $sum_tandem_two = 0;
my $counter = 0;
my $mp_limit = 25000;
my $pe_limit = 1000;
my $tandem_limit = 10000;
my $mp_target_too_small=0;
my $pe_target_too_small=0;
my $tandem_target_too_small=0;
my $min_map_q = 0;
my $log_filename = 0;
my %read_one_alignments;
my $reference_filename = "";
my $reference_min_size = 0;
my %reference_lengths;

&GetOptions(
'one:s'        => \$read_one,
'two:s'        => \$read_two,
'out:s'        => \$out_prefix,
'reference:s'  => \$reference_filename,
'mplimit:i'    => \$mp_limit,
'pelimit:i'    => \$pe_limit,
'refminsize:i' => \$reference_min_size,
'minmapq:i'    => \$min_map_q,
'log:s'        => \$log_filename
);

log_print("\nIn nextclip_sam_parse.pl\n");

die if not defined $read_one;
die if not defined $read_two;
die if not defined $out_prefix;

log_print("R1: $read_one\nR2: $read_one\nOut prefix: $out_prefix\nMP limit: $mp_limit\nPE limit: $pe_limit\nRefminsize: $reference_min_size\nMinmapq: $min_map_q\n");
log_print("Reference: $reference_filename\n") if (defined $reference_filename);

# Check for files
unless (-e $read_one) {
    log_and_screen("ERROR: nextclip_sam_parse.pl: Can't find R1 file $read_one\n");
    die;
}

unless (-e $read_two) {
    log_and_screen("ERROR: nextclip_sam_parse.pl: Can't find R2 file $read_two\n");
    die;
}

my $mp_out_filename = $out_prefix."_mp.txt";
my $pe_out_filename = $out_prefix."_pe.txt";
my $tandem_out_filename = $out_prefix."_tandem.txt";
my $unmapped_out_filename = $out_prefix."_unmapped.txt";

log_print("MP out filename: $mp_out_filename\nPE out filename: $pe_out_filename\nTandem out filename: $tandem_out_filename\nUnmapped out filename: $unmapped_out_filename\n");

if ($reference_min_size > 0) {
    if ($reference_filename eq "") {
        log_and_screen("WARNING: nextclip_sam_parse.pl: You've specified a reference minimum size, so you must specify a reference filename. Will reset refminsize to 0.");
        $reference_min_size=0;
    } else {
        load_reference_index();
    }
}

open(my $fh_mp, ">".$mp_out_filename);
open(my $fh_pe, ">".$pe_out_filename);
open(my $fh_tandem, ">".$tandem_out_filename);
open(my $fh_unmapped, ">".$unmapped_out_filename);

# Load R2 sam file and store
print "Opening $read_one\n";
open(my $fh_one, $read_one) or die "Can't open $read_two\n";
<$fh_one>; <$fh_one>;
while(<$fh_one>) {
    chomp(my $sam_line = $_);
    if ($sam_line !~ /^\@/) {
        my @fields = split(/\t/, $sam_line);
        my $id = $fields[0];
        my $flags = $fields[1];
 
        if (($flags & 256) == 0) {
            if (defined $read_one_alignments{$id}) {
                print "Warning: already seen ID $id\n";
                print "   Previous: ", $read_one_alignments{$id}, "\n";
                print "       This: ", $sam_line, "\n";
            } else {
                $read_one_alignments{$id} = $sam_line;
            }
        }
    }
}
close($fh_one);


# Load R1 sam file
print "Opening $read_two\n";
open(my $fh_two, $read_two) or die "Can't open $read_one\n";
<$fh_two>; <$fh_two>;

while(<$fh_two>) {
    chomp(my $sam_line_two = $_);
    if ($sam_line_two !~ /^\@/) {
        my @fields_two = split(/\t/, $sam_line_two);
        my $flags_two = $fields_two[1];
        my $id = $fields_two[0];
    
        if (($flags_two & 256) == 0) {
            if (defined $read_one_alignments{$id}) {
                $counter++;
                my $sam_line_one = $read_one_alignments{$id};
                my @fields_one = split(/\t/, $sam_line_one);
                my $flags_one = $fields_one[1];
                my $ref_one = $fields_one[2];
                my $ref_two = $fields_two[2];
                my $pos_one = $fields_one[3];
                my $pos_two = $fields_two[3];
                my $mapq_one = $fields_one[4];
                my $mapq_two = $fields_two[4];
                my $reverse_one = $flags_one & 16 ? 1:0;
                my $reverse_two = $flags_two & 16 ? 1:0;
                my $length_read_one = length($fields_one[9]);
                my $length_read_two = length($fields_two[9]);
                
                # Are both unmapped?
                if ((($flags_one & 4) == 4) && (($flags_two & 4) == 4)) {
                    $unmapped++;
                    print $fh_unmapped $id, "\t", $length_read_one, "\t", $length_read_two, "\t", $flags_one, "\t", $flags_two, "\t", $ref_one, "\t", $ref_two, "\n";
                }
                # Is just one read unmapped?
                elsif ((($flags_one & 4) == 4) || (($flags_two & 4) == 4)) {
                    $n_single_mapped++;
                    print $fh_unmapped $id, "\t", $length_read_one, "\t", $length_read_two, "\t", $flags_one, "\t", $flags_two, "\t", $ref_one, "\t", $ref_two, "\n";
                }
                # Is mapping quality too low
                elsif (($mapq_one < $min_map_q) || ($mapq_two < $min_map_q)) {
                    $n_poor_mapping++;
                    print $fh_unmapped $id, "\t", $length_read_one, "\t", $length_read_two, "\t", $flags_one, "\t", $flags_two, "\t", $ref_one, "\t", $ref_two, "\n";
                }
                # Do R1 and R2 map to a different reference?
                elsif ($ref_one ne $ref_two) {
                    $n_map_differently++;
                    print $fh_unmapped $id, "\t", $length_read_one, "\t", $length_read_two, "\t", $flags_one, "\t", $flags_two, "\t", $ref_one, "\t", $ref_two, "\n";
                }
                # Work out orientation - MP, PE or tandem
                else {
                    if ($pos_one <= $pos_two) {
                        my $insert = $pos_two - $pos_one;
                        if (($reverse_one == 1) && ($reverse_two == 0)) {
                            $n_mp++;
                            if (($insert >=0) && ($insert <= $mp_limit)) {
                                $sum_mp_one += $length_read_one;
                                $sum_mp_two += $length_read_two;
                                write_to_graph_file($fh_mp, $ref_one, $id, $insert, $length_read_one, $length_read_two);
                            } else {
                                $n_mp_out_of_range++;
                            }
                        } elsif (($reverse_one == 0) && ($reverse_two == 1)) {
                            $n_pe++;
                            if (($insert >= 0) && ($insert <= $pe_limit)) {
                                $sum_pe_one += $length_read_one;
                                $sum_pe_two += $length_read_two;
                                write_to_graph_file($fh_pe, $ref_one, $id, $insert, $length_read_one, $length_read_two);
                            } else {
                                $n_pe_out_of_range++;
                            }
                        } elsif ($reverse_one == $reverse_two) {
                            $n_tandem++;
                            if (($insert >= 0) && ($insert <= $tandem_limit)) {
                                $sum_tandem_one += $length_read_one;
                                $sum_tandem_two += $length_read_two;
                                write_to_graph_file($fh_tandem, $ref_one, $id, $insert, $length_read_one, $length_read_two);
                            } else {
                                $n_tandem_out_of_range++;
                            }
                        }
                    } elsif ($pos_one > $pos_two) {
                        my $insert = $pos_one - $pos_two;
                        if (($reverse_one == 0) && ($reverse_two == 1)) {
                            $n_mp++;
                            if (($insert >= 0) && ($insert <= $mp_limit)) {
                                $sum_mp_one += $length_read_one;
                                $sum_mp_two += $length_read_two;
                                write_to_graph_file($fh_mp, $ref_one, $id, $insert, $length_read_one, $length_read_two);
                            } else {
                                $n_mp_out_of_range++;
                            }
                        } elsif (($reverse_one == 1) && ($reverse_two == 0)) {
                            $n_pe++;
                            if (($insert >= 0) && ($insert <= $pe_limit)) {
                                $sum_pe_one += $length_read_one;
                                $sum_pe_two += $length_read_two;
                                write_to_graph_file($fh_pe, $ref_one, $id, $insert, $length_read_one, $length_read_two);
                           } else {
                               $n_pe_out_of_range++;
                           }
                        } elsif ($reverse_one == $reverse_two) {
                            $n_tandem++;
                            if (($insert >= 0) && ($insert <= $tandem_limit)) {
                                $sum_tandem_one += $length_read_one;
                                $sum_tandem_two += $length_read_two;
                                write_to_graph_file($fh_tandem, $ref_one, $id, $insert, $length_read_one, $length_read_two);
                            } else {
                                $n_tandem_out_of_range++;
                            }
                        }
                    }
                }
            }
        } else {
            $n_dupes++;
        }
    }
}

close($fh_one);
close($fh_mp);
close($fh_pe);
close($fh_tandem);
close($fh_unmapped);

my $pc_mp = 0;
my $pc_mp_in_range = 0;
my $pc_mp_out_of_range = 0;
my $av_mp_one = 0;
my $av_mp_two = 0;
if ($n_mp > 0) {
    $pc_mp = 100.0*$n_mp/$counter;
    $pc_mp_out_of_range = 100.0*$n_mp_out_of_range/$counter;
    $pc_mp_in_range = 100.0*($n_mp - $n_mp_out_of_range)/$counter;
    $av_mp_one = $sum_mp_one/$n_mp;
    $av_mp_two = $sum_mp_two/$n_mp;
};

my $pc_pe = 0;
my $pc_pe_in_range = 0;
my $pc_pe_out_of_range = 0;
my $av_pe_one = 0;
my $av_pe_two = 0;
if ($n_pe > 0) {
    $pc_pe = 100.0*$n_pe/$counter;
    $pc_pe_out_of_range = 100.0*$n_pe_out_of_range/$counter;
    $pc_pe_in_range = 100.0*($n_pe - $n_pe_out_of_range)/$counter;
    $av_pe_one = $sum_pe_one/$n_pe;
    $av_pe_two = $sum_pe_two/$n_pe;
};

my $pc_tandem = 0;
my $pc_tandem_in_range = 0;
my $pc_tandem_out_of_range = 0;
my $av_tandem_one = 0;
my $av_tandem_two = 0;
if ($n_tandem > 0) {
    $pc_tandem = 100.0*$n_tandem/$counter;
    $pc_tandem_out_of_range = 100.0*$n_tandem_out_of_range/$counter;
    $pc_tandem_in_range = 100.0*($n_tandem - $n_tandem_out_of_range)/$counter;
    $av_tandem_one = $sum_tandem_one/$n_tandem;
    $av_tandem_two = $sum_tandem_two/$n_tandem;
};

my $total_good_maps = $n_mp + $n_pe + $n_tandem;
my $total_good_maps_out_of_range = $n_mp_out_of_range + $n_pe_out_of_range + $n_tandem_out_of_range;
my $total_good_maps_in_range = $total_good_maps - $total_good_maps_out_of_range;
my $total_bad_maps = $unmapped + $n_single_mapped + $n_map_differently + $n_poor_mapping;

print "\n";
print  "                           MP limit:\t", $mp_limit, "\n";
print  "                           PE limit:\t", $pe_limit, "\n";
print  "                       Tandem limit:\t", $tandem_limit, "\n";
print  "      Reference minimum contig size:\t", $reference_min_size, "\n";
print  "                     MAPQ threshold:\t", $min_map_q, "\n";
print  "               Number of alignments:\t", $counter, "\n";

if ($counter > 0) {
    printf "\n";
    printf "     Total number with good mapping:\t%i\t%.2f\n", $total_good_maps, 100.0*$total_good_maps/$counter;
    printf "           Total good maps in limit:\t%i\t%.2f\n", $total_good_maps_in_range, 100.0*$total_good_maps_in_range/$counter;
    printf "       Total good maps out of limit:\t%i\t%.2f\n", $total_good_maps_out_of_range, 100.0*$total_good_maps_out_of_range/$counter;
    printf "\n";
    printf "          Number with both unmapped:\t%i\t%.2f\n", $unmapped, 100.0*$unmapped/$counter;
    printf "      Number with one read unmapped:\t%i\t%.2f\n", $n_single_mapped, 100.0*$n_single_mapped/$counter;
    printf "        Number that map differently:\t%i\t%.2f\n", $n_map_differently, 100.0*$n_map_differently/$counter;
    printf "           Number with poor mapping:\t%i\t%.2f\n", $n_poor_mapping, 100.0*$n_poor_mapping/$counter;
    printf "      Total number with bad mapping:\t%i\t%.2f\n", $total_bad_maps, 100.0*$total_bad_maps/$counter;
    printf "\n";
    printf "                Number of mate pair:\t%i\t%.2f\t%.2f\t%.2f\n", $n_mp, $pc_mp, $av_mp_one, $av_mp_two;
    printf "       Number of mate pair in limit:\t%i\t%.2f\n", ($n_mp - $n_mp_out_of_range), $pc_mp_in_range;
    printf "    (MPs aligning to small contigs):\t%i\n", $mp_target_too_small if ($reference_min_size > 0);
    printf "   Number of mate pair out of limit:\t%i\t%.2f\n", $n_mp_out_of_range, $pc_mp_out_of_range;
    printf "\n";
    printf "                 Number of pair end:\t%i\t%.2f\t%.2f\t%.2f\n", $n_pe, $pc_pe, $av_pe_one, $av_pe_two;
    printf "        Number of pair end in limit:\t%i\t%.2f\n", ($n_pe - $n_pe_out_of_range), $pc_pe_in_range;
    printf "    (PEs aligning to small contigs):\t%i\n", $pe_target_too_small if ($reference_min_size > 0);
    printf "    Number of pair end out of limit:\t%i\t%.2f\n", $n_pe_out_of_range, $pc_pe_out_of_range;
    printf "\n";
    printf "                   Number of tandem:\t%i\t%.2f\t%.2f\t%.2f\n", $n_tandem, $pc_tandem, $av_tandem_one, $av_tandem_two;
    printf "          Number of tandem in limit:\t%i\t%.2f\n", ($n_tandem - $n_tandem_out_of_range), $pc_tandem_in_range;
    printf " (Tandem aligning to small contigs):\t%i\n", $tandem_target_too_small if ($reference_min_size > 0);
    printf "      Number of tandem out of limit:\t%i\t%.2f\n", $n_tandem_out_of_range, $pc_tandem_out_of_range;
    printf "\n";
    print  "                    Number of dupes:\t", $n_dupes, "\n";
    print  "                 Recalculated total:\t", $unmapped+$n_single_mapped+$n_map_differently+$n_mp+$n_pe+$n_tandem, "\n";
} else {
    log_and_screen("\nERROR: No alignments - did the BWA stage work correctly?\n");
}

if ($reference_min_size > 0) {
    print "\nNote: Pairs aligning to small (<",$reference_min_size,") contigs have not been written to output.\n";
}

if ($n_mp < 10) {
   log_and_screen("\nERROR: Number of MP alignments (".$n_mp.") too small\n");
} elsif ($n_mp < 1000) {
   log_and_screen("\nWARNING: Number of MP alignments (".$n_mp.") too small\n");
}

printf "\nDONE\n";


sub write_to_graph_file
{
    my $fh=$_[0];
    my $ref=$_[1];
    my $id=$_[2];
    my $insert=$_[3];
    my $length_one=$_[4];
    my $length_two=$_[5];
    my $length = 1;
    
    if ($reference_min_size > 0) {
        if (defined $reference_lengths{$ref}) {
            $length = $reference_lengths{$ref};
        } else {
            die "ERROR: nextclip_sam_parse.pl: No length defined for ID [$ref]\n";
        }
    }

    if ($length >= $reference_min_size) {
        print $fh $id, "\t", $insert, "\t", $length_one, "\t", $length_two, "\n";
    } else {
        if ($fh == $fh_mp) {
            $mp_target_too_small++;
        } elsif ($fh == $fh_pe) {
            $pe_target_too_small++;
        } elsif ($fh == $fh_tandem) {
            $tandem_target_too_small++;
        }
    }
}

sub write_dummy_entries
{
    print $fh_mp "DummyId\t0\t0\t0\n";
    print $fh_pe "DummyId\t0\t0\t0\n";
    print $fh_tandem "DummyId\t0\t0\t0\n";
}

sub load_reference_index
{
    my $reference_index = $reference_filename.".nextclip";

    log_and_screen("Opening index file $reference_index\n");
    log_and_screen("ERROR: nextclip_sam_parse.pl: Can't find $reference_index - have you run nextclip_index_reference.pl as per the manual?\n") unless (-e $reference_index);
    open(INDEXFILE, $reference_index) or die "Can't open reference index $reference_index\n";

    while(<INDEXFILE>) {
        chomp(my $line = $_);
        my @arr = split(/\t/, $line);
        my $id = $arr[0];
        my $length = $arr[1];
        
        if (defined $reference_lengths{$id}) {
            print "WARNING: nextclip_sam_parse.pl: Reference ID $id already used.\n";
        } else {
            $reference_lengths{$id} = $length;
        }
    }
    close(INDEXFILE);
}

# Write to log and screen
sub log_and_screen
{
    log_print(@_);
    print @_;
}

sub log_print
{
    if (defined $log_filename) {
        open(my $log_fh, ">>".$log_filename);
        print $log_fh @_;
        close($log_fh);
    }
}
