#!/usr/bin/perl -w

# Script:  nextclip_make_report.pl
# Purpose: Produce PDF report of NextClip analysis
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

use warnings;
use strict;
use Getopt::Long;

my $library_dir;
my $log_filename;
my $library_name;
my $organism;
my $nextclip_file;
my $reference;
my %attributes;
my %att;

# Get command line options
&GetOptions(
'libdir:s'     => \$library_dir,
'log:s'        => \$log_filename,
'libname:s'    => \$library_name,
'nextclip:s'   => \$nextclip_file,
'organism:s'   => \$organism,
'reference:s'  => \$reference
);

die "ERROR: You must specify a -libdir option.\n" if (not defined $library_dir);
die "ERROR: You must specify a -log option.\n" if (not defined $log_filename);
die "ERROR: You must specify a -libname option.\n" if (not defined $library_name);
die "ERROR: You must specify a -organism option.\n" if (not defined $organism);
die "ERROR: You must specify a -reference option.\n" if (not defined $reference);
die "ERROR: You must specify a -nextclip option.\n" if (not defined $nextclip_file);

my $title=$library_name;
my $output_dir=$library_dir."/latex";
my $latex_file=$library_dir."/latex/".$library_name.".tex";
my $parseable_file=$library_dir."/analysis/".$library_name."_analysis.txt";
my $latex_fh;
my $parseable_fh;

open($latex_fh, ">".$latex_file) or die "ERROR: Can't open latex file $latex_file\n";
open($parseable_fh, ">".$parseable_file) or die "ERROR: Can't open parseable file $parseable_file\n";

$title =~ s/_/\\_/g;

load_reference_index();

if ($organism ne "Unknown") {
    $title="NextClip report for $title ($organism)"
}

log_and_screen("\n");
log_and_screen("Report builder\n");
log_and_screen("\n");
log_and_screen("Library name is ".$library_name."\n");
log_and_screen("Output dir is ".$output_dir."\n");
log_and_screen("NextClip log is ".$nextclip_file."\n");
log_and_screen("Organism is ".$organism."\n");
log_and_screen("Reference is ".$reference."\n");
log_and_screen("LaTeX file is ".$latex_file."\n");

output_latex_header();

print $latex_fh "\\section*{\\large{".$title."}}\n";

get_main_stats();
write_overall_section();

print $latex_fh "\\subsection*{Fragments with junction adaptor in R1 and R2 (A)}\n";
print $latex_fh "\\vspace{-3mm}\n";
output_category_report($library_dir."/logs/parse_A.log", "A");

print $latex_fh "\\subsection*{Fragments with junction adaptor in R2 only (B)}\n";
print $latex_fh "\\vspace{-3mm}\n";
output_category_report($library_dir."/logs/parse_B.log", "B");

print $latex_fh "\\clearpage\n";

print $latex_fh "\\subsection*{Fragments with junction adaptor in R1 only (C)}\n";
print $latex_fh "\\vspace{-3mm}\n";
output_category_report($library_dir."/logs/parse_C.log", "C");

print $latex_fh "\\subsection*{Fragments where neither read contain the junction adaptor (D)}\n";
print $latex_fh "\\vspace{-3mm}\n";
output_category_report($library_dir."/logs/parse_D.log", "D");

write_gc_content_section();

print $latex_fh "\\clearpage\n";

write_shortest_pair_section();
write_clipped_read_length_section();
write_duplication_section();
write_ambiguous_bases_section();
write_external_adaptor_section();
write_notes_section();
output_latex_footer();

close($parseable_fh);
close($latex_fh);

log_and_screen("Trying to build PDF\n");

make_pdf();

exit(0);

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

# ----------------------------------------------------------------------
# Function: output_latex_header
# Purpose:  Ouput preamble to LaTeX file
# ----------------------------------------------------------------------
sub output_latex_header
{
    print $latex_fh "\\documentclass[a4paper,11pt,oneside]{article}\n";
    print $latex_fh "\\usepackage{graphicx}\n";
    print $latex_fh "\\usepackage{url}\n";
    print $latex_fh "\\usepackage{subcaption}\n";
    print $latex_fh "\\usepackage{rotating}\n";
    print $latex_fh "\\usepackage{color}\n";
    print $latex_fh "\\usepackage[compact]{titlesec}\n";
    print $latex_fh "\\usepackage[portrait,top=1cm, bottom=2cm, left=1cm, right=1cm]{geometry}\n";
    print $latex_fh "\\begin{document}\n";
    print $latex_fh "\\renewcommand*{\\familydefault}{\\sfdefault}\n";
}

# ----------------------------------------------------------------------
# Function: output_latex_footer
# Purpose:  End LaTeX document
# ----------------------------------------------------------------------
sub output_latex_footer
{
    print $latex_fh "\\end{document}\n";
}

# ----------------------------------------------------------------------
# Function: make_pdf
# Purpose:  Run pdflatex
# ----------------------------------------------------------------------
sub make_pdf
{
    my $pdf_file=$latex_file;
    $pdf_file =~ s/\.tex/\.pdf/;

    log_and_screen "Writing PDF $pdf_file\n";
    system("pdflatex --file-line-error --shell-escape --output-directory ".$output_dir." --interaction=batchmode ".$latex_file);

    log_and_screen "\nNextClip analysis complete.\n";
}

# ----------------------------------------------------------------------
# Function: output_category_report
# Purpose:  Output LaTeX for each category (A, B, C etc.)
# ----------------------------------------------------------------------
sub output_category_report
{
    my $sumfile = $_[0];
    my $suffix = $_[1];

    log_and_screen "ERROR: Can't find $sumfile - something must have gone wrong with an earlier process.\n" unless (-e $sumfile);

    undef %att;
    
    open(SUMFILE, $sumfile) or die "ERROR: Can't open $sumfile\n";
    while(<SUMFILE>) {
        chomp(my $line = $_);
        if ($line =~ /MP limit:\t(\d+)/) {
            $att{"mp_limit"} = $1;
        } elsif ($line =~ /PE limit:\t(\d+)/) {
            $att{"pe_limit"} = $1;
        } elsif ($line =~ /Tandem limit:\t(\d+)/) {
            $att{"tandem_limit"} = $1;
        } elsif ($line =~ /MAPQ threshold:\t(\d+)/) {
            $att{"mapq_threshold"} = $1;
        } elsif ($line =~ /Reference minimum contig size:\t(\d+)/) {
            $att{"ref_min"} = $1;
        } elsif ($line =~ /Number of alignments:\t(\d+)/) {
            $att{"num_alignments"} = $1;
        } elsif ($line =~ /Number of mate pair:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/) {
            $att{"n_mp"} = $1;
            $att{"pc_mp"} = $2;
            $att{"mp_r1_av"} = $3;
            $att{"mp_r2_av"} = $4;
        } elsif ($line =~ /Number of mate pair in limit:\t(\d+)\t(\S+)/) {
            $att{"n_mp_in_range"} = $1;
            $att{"pc_mp_in_range"} = $2;
        } elsif ($line =~ /Number of mate pair out of limit:\t(\d+)\t(\S+)/) {
            $att{"n_mp_out_of_range"} = $1;
            $att{"pc_mp_out_of_range"} = $2;
        } elsif ($line =~ /Number of pair end:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/) {
            $att{"n_pe"} = $1;
            $att{"pc_pe"} = $2;
            $att{"pe_r1_av"} = $3;
            $att{"pe_r2_av"} = $4;
        } elsif ($line =~ /Number of pair end in limit:\t(\d+)\t(\S+)/) {
            $att{"n_pe_in_range"} = $1;
            $att{"pc_pe_in_range"} = $2;
        } elsif ($line =~ /Number of pair end out of limit:\t(\d+)\t(\S+)/) {
            $att{"n_pe_out_of_range"} = $1;
            $att{"pc_pe_out_of_range"} = $2;
        } elsif ($line =~ /Number of tandem:\t(\d+)\t(\S+)\t(\S+)\t(\S+)/) {
            $att{"n_tandem"} = $1;
            $att{"pc_tandem"} = $2;
            $att{"tandem_r1_av"} = $3;
            $att{"tandem_r2_av"} = $4;
        } elsif ($line =~ /Number of tandem in limit:\t(\d+)\t(\S+)/) {
            $att{"n_tandem_in_range"} = $1;
            $att{"pc_tandem_in_range"} = $2;
        } elsif ($line =~ /Number of tandem out of limit:\t(\d+)\t(\S+)/) {
            $att{"n_tandem_out_of_range"} = $1;
            $att{"pc_tandem_out_of_range"} = $2;
        } elsif ($line =~ /Number with both unmapped:\t(\d+)\t(\S+)/) {
            $att{"n_unmapped"} = $1;
            $att{"pc_unmapped"} = $2;
        } elsif ($line =~ /Number with one read unmapped:\t(\d+)\t(\S+)/) {
            $att{"n_singlemapped"} = $1;
            $att{"pc_singlemapped"} = $2;
        } elsif ($line =~ /Number that map differently:\t(\d+)\t(\S+)/) {
            $att{"n_differentmapped"} = $1;
            $att{"pc_differentmapped"} = $2;
        } elsif ($line =~ /Number with poor mapping:\t(\d+)\t(\S+)/) {
            $att{"n_poormapping"} = $1;
            $att{"pc_poormapping"} = $2;
        } elsif ($line =~ /Total number with good mapping:\t(\d+)\t(\S+)/) {
            $att{"good_maps"} = $1." (".$2."\\\%)";
            $att{"n_good_maps"} = $1;
            $att{"pc_good_maps"} = $2;
        } elsif ($line =~ /Total good maps in limit:\t(\d+)\t(\S+)/) {
            $att{"good_maps_in_range"} = $1." (".$2."\\\%)";
            $att{"n_good_maps_in_range"} = $1;
            $att{"pc_good_maps_in_range"} = $2;
        } elsif ($line =~ /Total good maps out of limit:\t(\d+)\t(\S+)/) {
            $att{"good_maps_out_of_range"} = $1." (".$2."\\\%)";
            $att{"n_good_maps_out_of_range"} = $1;
            $att{"pc_good_maps_out_of_range"} = $2;
        } elsif ($line =~ /Total number with bad mapping:\t(\d+)\t(\S+)/) {
            $att{"bad_maps"} = $1." (".$2."\\\%)";
            $att{"n_bad_maps"} = $1;
            $att{"pc_bad_maps"} = $2;
        }
    }
    close(SUMFILE);
   
    if ((not defined $att{"n_mp"}) || (not defined $att{"n_pe"}) || (not defined $att{"n_tandem"}) || (not defined $att{"bad_maps"})) {
        log_and_screen("\nERROR: Missing fields in ".$sumfile." - did BWA run correctly?\n");
        die;
    }
 
    if ($suffix eq "A") {
        print $parseable_fh "MPLimit:".$att{"mp_limit"}."\n";
        print $parseable_fh "PELimit:".$att{"pe_limit"}."\n";
        print $parseable_fh "TandemLimit:".$att{"tandem_limit"}."\n";
        print $parseable_fh "MapQThreshold:".$att{"mapq_threshold"}."\n";
        print $parseable_fh "RefMin:".$att{"ref_min"}."\n";
    }
    
    print $parseable_fh "Cat".$suffix."NumAlignments:".$att{"num_alignments"}."\n";
    print $parseable_fh "Cat".$suffix."NumMatePair:".$att{"n_mp"}."\n";
    print $parseable_fh "Cat".$suffix."PcMatePair:".$att{"pc_mp"}."\n";
    print $parseable_fh "Cat".$suffix."MatePairR1Av:".$att{"mp_r1_av"}."\n";
    print $parseable_fh "Cat".$suffix."MatePairR2Av:".$att{"mp_r2_av"}."\n";
    print $parseable_fh "Cat".$suffix."NumMatePairInRange:".$att{"n_mp_in_range"}."\n";
    print $parseable_fh "Cat".$suffix."PcMatePairInRange:".$att{"pc_mp_in_range"}."\n";
    print $parseable_fh "Cat".$suffix."NumMatePairOutOfRange:".$att{"n_mp_out_of_range"}."\n";
    print $parseable_fh "Cat".$suffix."PcMatePairOutOfRange:".$att{"pc_mp_out_of_range"}."\n";
    print $parseable_fh "Cat".$suffix."NumPairEnd:".$att{"n_pe"}."\n";
    print $parseable_fh "Cat".$suffix."PcPairEnd:".$att{"pc_pe"}."\n";
    print $parseable_fh "Cat".$suffix."PairEndR1Av:".$att{"pe_r1_av"}."\n";
    print $parseable_fh "Cat".$suffix."PairEndR2Av:".$att{"pe_r2_av"}."\n";
    print $parseable_fh "Cat".$suffix."NumPairEndInRange:".$att{"n_pe_in_range"}."\n";
    print $parseable_fh "Cat".$suffix."PcPairEndInRange:".$att{"pc_pe_in_range"}."\n";
    print $parseable_fh "Cat".$suffix."NumPairEndOutOfRange:".$att{"n_pe_out_of_range"}."\n";
    print $parseable_fh "Cat".$suffix."PcPairEndOutOfRange:".$att{"pc_pe_out_of_range"}."\n";
    print $parseable_fh "Cat".$suffix."NumTandem:".$att{"n_tandem"}."\n";
    print $parseable_fh "Cat".$suffix."PcTandem:".$att{"pc_tandem"}."\n";
    print $parseable_fh "Cat".$suffix."TandemR1Av:".$att{"tandem_r1_av"}."\n";
    print $parseable_fh "Cat".$suffix."TandemR2Av:".$att{"tandem_r2_av"}."\n";
    print $parseable_fh "Cat".$suffix."NumTandemInRange:".$att{"n_tandem_in_range"}."\n";
    print $parseable_fh "Cat".$suffix."PcTandemInRange:".$att{"pc_tandem_in_range"}."\n";
    print $parseable_fh "Cat".$suffix."NumTandemOutOfRange:".$att{"n_tandem_out_of_range"}."\n";
    print $parseable_fh "Cat".$suffix."PcTandemOutOfRange:".$att{"pc_tandem_out_of_range"}."\n";
    print $parseable_fh "Cat".$suffix."NumUnmapped:".$att{"n_unmapped"}."\n";
    print $parseable_fh "Cat".$suffix."PcUnmapped:".$att{"pc_unmapped"}."\n";
    print $parseable_fh "Cat".$suffix."NumSingleMapped:".$att{"n_singlemapped"}."\n";
    print $parseable_fh "Cat".$suffix."PcSingleMapped:".$att{"pc_singlemapped"}."\n";
    print $parseable_fh "Cat".$suffix."NumDifferentMapped:".$att{"n_differentmapped"}."\n";
    print $parseable_fh "Cat".$suffix."PcDifferentMapped:".$att{"pc_differentmapped"}."\n";
    print $parseable_fh "Cat".$suffix."NumPoorlyMapped:".$att{"n_poormapping"}."\n";
    print $parseable_fh "Cat".$suffix."PcPoorlyMapped:".$att{"pc_poormapping"}."\n";
    print $parseable_fh "Cat".$suffix."NumGoodMaps:".$att{"n_good_maps"}."\n";
    print $parseable_fh "Cat".$suffix."PcGoodMaps:".$att{"pc_good_maps"}."\n";
    print $parseable_fh "Cat".$suffix."NumGoodMapsInRange:".$att{"n_good_maps_in_range"}."\n";
    print $parseable_fh "Cat".$suffix."PcGoodMapsInRange:".$att{"pc_good_maps_in_range"}."\n";
    print $parseable_fh "Cat".$suffix."NumGoodMapsOutOfRange:".$att{"n_good_maps_out_of_range"}."\n";
    print $parseable_fh "Cat".$suffix."PcGoodMapsOutOfRange:".$att{"pc_good_maps_out_of_range"}."\n";
    print $parseable_fh "Cat".$suffix."NumBadMaps:".$att{"n_bad_maps"}."\n";
    print $parseable_fh "Cat".$suffix."PcBadMaps:".$att{"pc_bad_maps"}."\n";

    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{l c c c c c}\n";
    print $latex_fh "& {\\bf Total} & {\\bf In range} & {\\bf Out of range} & {\\bf R1 mean} & {\\bf R2 mean} \\\\\n";
    print $latex_fh "{\\bf Pairs producing good mappings} & ".$att{"good_maps"}." & ".$att{"good_maps_in_range"}." & ".$att{"good_maps_out_of_range"}." & & \\\\\n";
    print $latex_fh "\\hspace{5mm}In Mate Pair orientation & ".$att{"n_mp"}." (".$att{"pc_mp"}."\\\%) & ".$att{"n_mp_in_range"}." (".$att{"pc_mp_in_range"}."\\\%) & ".$att{"n_mp_out_of_range"}." (".$att{"pc_mp_out_of_range"}."\\\%) & ".$att{"mp_r1_av"}." & ".$att{"mp_r2_av"}." \\\\\n";
    print $latex_fh "\\hspace{5mm}In Pair End orientation & ".$att{"n_pe"}." (".$att{"pc_pe"}."\\\%) & ".$att{"n_pe_in_range"}." (".$att{"pc_pe_in_range"}."\\\%) & ".$att{"n_pe_out_of_range"}." (".$att{"pc_pe_out_of_range"}."\\\%) & ".$att{"pe_r1_av"}." & ".$att{"pe_r2_av"}." \\\\\n";
    print $latex_fh "\\hspace{5mm}In Tandem orientation & ".$att{"n_tandem"}." (".$att{"pc_tandem"}."\\\%) & ".$att{"n_tandem_in_range"}." (".$att{"pc_tandem_in_range"}."\\\%) & ".$att{"n_tandem_out_of_range"}." (".$att{"pc_tandem_out_of_range"}."\\\%) & ".$att{"tandem_r1_av"}." & ".$att{"tandem_r2_av"}." \\\\\n";
    print $latex_fh "{\\bf Pairs producing bad mappings} & ".$att{"bad_maps"}." & & & \\\\\n";
    print $latex_fh "\\hspace{5mm}Both reads unmapped & ".$att{"n_unmapped"}." (".$att{"pc_unmapped"}."\\\%) & & & \\\\\n";
    print $latex_fh "\\hspace{5mm}One read unmapped & ".$att{"n_singlemapped"}." (".$att{"pc_singlemapped"}."\\\%) & & & \\\\\n";
    print $latex_fh "\\hspace{5mm}Reads map to different ID & ".$att{"n_differentmapped"}." (".$att{"pc_differentmapped"}."\\\%) & & & \\\\\n";
    print $latex_fh "\\hspace{5mm}Reads with MAPQ {\\textless} ".$att{"mapq_threshold"}." & ".$att{"n_poormapping"}." (".$att{"pc_poormapping"}."\\\%) & & & \\\\\n";
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";
    print $latex_fh "\\vspace{-5mm}\n";

    print $latex_fh "\\begin{figure}[h!]\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\begin{subfigure}{.3\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_mp_".$suffix.".pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Mate pair}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\begin{subfigure}{.3\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_pe_".$suffix.".pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Pair end}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\begin{subfigure}{.3\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_tandem_".$suffix.".pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Tandem}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\end{figure}\n";
}

# ----------------------------------------------------------------------
# Function: get_main_stats
# Purpose:  Get main NextClip stats
# ----------------------------------------------------------------------
sub get_main_stats
{
    log_and_screen "ERROR: Can't find $nextclip_file - did NextClip fail to run?\n" unless (-e $nextclip_file);
    
    open(NCFILE, $nextclip_file) or die "ERROR: Can't open $nextclip_file\n";
    while(<NCFILE>) {
        chomp(my $line = $_);
        if ($line =~ /Minimum read size: (\d+)/) {
            $attributes{"minimum_read_size"} = $1;
        } elsif ($line =~ /Number of read pairs: (\d+)/) {
            $attributes{"number_of_pairs"} = $1;
        } elsif ($line =~ /R1 Num reads with adaptor: (\d+)\t(\S+)/) {
            $attributes{"r1_with_adaptor"} = $1." (".$2."\\\%)";
            $attributes{"n_r1_with_adaptor"} = $1;
            $attributes{"pc_r1_with_adaptor"} = $2;
        } elsif ($line =~ /R1 Num with external also: (\d+)\t(\S+)/) {
            $attributes{"r1_with_adaptor_and_external"} = $1." (".$2."\\\%)";
            $attributes{"n_r1_with_adaptor_and_external"} = $1;
            $attributes{"pc_r1_with_adaptor_and_external"} = $2;
        } elsif ($line =~ /R1 long adaptor reads: (\d+)\t(\S+)/) {
            $attributes{"r1_with_adaptor_long"} = $1." (".$2."\\\%)";
            $attributes{"n_r1_with_adaptor_long"} = $1;
            $attributes{"pc_r1_with_adaptor_long"} = $2;
        } elsif ($line =~ /R1 reads too short: (\d+)\t(\S+)/) {
            $attributes{"r1_with_adaptor_short"} = $1." (".$2."\\\%)";
            $attributes{"n_r1_with_adaptor_short"} = $1;
            $attributes{"pc_r1_with_adaptor_short"} = $2;
        } elsif ($line =~ /R1 Num reads no adaptor: (\d+)\t(\S+)/) {
            $attributes{"r1_without_adaptor"} = $1." (".$2."\\\%)";
            $attributes{"n_r1_without_adaptor"} = $1;
            $attributes{"pc_r1_without_adaptor"} = $2;
        } elsif ($line =~ /R1 no adaptor but external: (\d+)\t(\S+)/) {
            $attributes{"r1_without_adaptor_but_external"} = $1." (".$2."\\\%)";
            $attributes{"n_r1_without_adaptor_but_external"} = $1;
            $attributes{"pc_r1_without_adaptor_but_external"} = $2;
        } elsif ($line =~ /R2 Num reads with adaptor: (\d+)\t(\S+)/) {
            $attributes{"r2_with_adaptor"} = $1." (".$2."\\\%)";
            $attributes{"n_r2_with_adaptor"} = $1;
            $attributes{"pc_r2_with_adaptor"} = $2;
        } elsif ($line =~ /R2 Num with external also: (\d+)\t(\S+)/) {
            $attributes{"r2_with_adaptor_and_external"} = $1." (".$2."\\\%)";
            $attributes{"n_r2_with_adaptor_and_external"} = $1;
            $attributes{"pc_r2_with_adaptor_and_external"} = $2;
        } elsif ($line =~ /R2 long adaptor reads: (\d+)\t(\S+)/) {
            $attributes{"r2_with_adaptor_long"} = $1." (".$2."\\\%)";
            $attributes{"n_r2_with_adaptor_long"} = $1;
            $attributes{"pc_r2_with_adaptor_long"} = $2;
        } elsif ($line =~ /R2 reads too short: (\d+)\t(\S+)/) {
            $attributes{"r2_with_adaptor_short"} = $1." (".$2."\\\%)";
            $attributes{"n_r2_with_adaptor_short"} = $1;
            $attributes{"pc_r2_with_adaptor_short"} = $2;
        } elsif ($line =~ /R2 Num reads no adaptor: (\d+)\t(\S+)/) {
            $attributes{"r2_without_adaptor"} = $1." (".$2."\\\%)";
            $attributes{"n_r2_without_adaptor"} = $1;
            $attributes{"pc_r2_without_adaptor"} = $2;
        } elsif ($line =~ /R2 no adaptor but external: (\d+)\t(\S+)/) {
            $attributes{"r2_without_adaptor_but_external"} = $1." (".$2."\\\%)";
            $attributes{"n_r2_without_adaptor_but_external"} = $1;
            $attributes{"pc_r2_without_adaptor_but_external"} = $2;
        } elsif ($line =~ /Total pairs in category A: (\d+)\t(\S+)/) {
            $attributes{"cat_a_total"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_a_total"} = $1;
            $attributes{"pc_cat_a_total"} = $2;
        } elsif ($line =~ /A pairs long enough: (\d+)\t(\S+)/) {
            $attributes{"cat_a_long"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_a_long"} = $1;
            $attributes{"pc_cat_a_long"} = $2;
        } elsif ($line =~ /A pairs too short: (\d+)\t(\S+)/) {
            $attributes{"cat_a_short"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_a_short"} = $1;
            $attributes{"pc_cat_a_short"} = $2;
        } elsif ($line =~ /A external clip in 1 or both: (\d+)\t(\S+)/) {
            $attributes{"cat_a_external"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_a_external"} = $1;
            $attributes{"pc_cat_a_external"} = $2;
        } elsif ($line =~ /A bases before clipping: (\d+)/) {
            $attributes{"cat_a_bases_before_clipping"} = $1;
        } elsif ($line =~ /A total bases written: (\d+)/) {
            $attributes{"cat_a_bases_written"} = $1;
        } elsif ($line =~ /Total pairs in category B: (\d+)\t(\S+)/) {
            $attributes{"cat_b_total"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_b_total"} = $1;
            $attributes{"pc_cat_b_total"} = $2;
        } elsif ($line =~ /B pairs long enough: (\d+)\t(\S+)/) {
            $attributes{"cat_b_long"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_b_long"} = $1;
            $attributes{"pc_cat_b_long"} = $1;
        } elsif ($line =~ /B pairs too short: (\d+)\t(\S+)/) {
            $attributes{"cat_b_short"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_b_short"} = $1;
            $attributes{"pc_cat_b_short"} = $2;
        } elsif ($line =~ /B external clip in 1 or both: (\d+)\t(\S+)/) {
            $attributes{"cat_b_external"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_b_external"} = $1;
            $attributes{"pc_cat_b_external"} = $2;
        } elsif ($line =~ /B bases before clipping: (\d+)/) {
            $attributes{"cat_b_bases_before_clipping"} = $1;
        } elsif ($line =~ /B total bases written: (\d+)/) {
            $attributes{"cat_b_bases_written"} = $1;
        } elsif ($line =~ /Total pairs in category C: (\d+)\t(\S+)/) {
            $attributes{"cat_c_total"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_c_total"} = $1;
            $attributes{"pc_cat_c_total"} = $2;
        } elsif ($line =~ /C pairs long enough: (\d+)\t(\S+)/) {
            $attributes{"cat_c_long"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_c_long"} = $1;
            $attributes{"pc_cat_c_long"} = $2;
        } elsif ($line =~ /C pairs too short: (\d+)\t(\S+)/) {
            $attributes{"cat_c_short"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_c_short"} = $1;
            $attributes{"pc_cat_c_short"} = $2;
        } elsif ($line =~ /C external clip in 1 or both: (\d+)\t(\S+)/) {
            $attributes{"cat_c_external"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_c_external"} = $1;
            $attributes{"pc_cat_c_external"} = $2;
        } elsif ($line =~ /C bases before clipping: (\d+)/) {
            $attributes{"cat_c_bases_before_clipping"} = $1;
        } elsif ($line =~ /C total bases written: (\d+)/) {
            $attributes{"cat_c_bases_written"} = $1;
        } elsif ($line =~ /Total pairs in category D: (\d+)\t(\S+)/) {
            $attributes{"cat_d_total"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_d_total"} = $1;
            $attributes{"pc_cat_d_total"} = $2;
        } elsif ($line =~ /D pairs long enough: (\d+)\t(\S+)/) {
            $attributes{"cat_d_long"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_d_long"} = $1;
            $attributes{"pc_cat_d_long"} = $2;
        } elsif ($line =~ /D pairs too short: (\d+)\t(\S+)/) {
            $attributes{"cat_d_short"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_d_short"} = $1;
            $attributes{"pc_cat_d_short"} = $2;
        } elsif ($line =~ /D external clip in 1 or both: (\d+)\t(\S+)/) {
            $attributes{"cat_d_external"} = $1." (".$2."\\\%)";
            $attributes{"n_cat_d_external"} = $1;
            $attributes{"pc_cat_d_external"} = $2;
        } elsif ($line =~ /D bases before clipping: (\d+)/) {
            $attributes{"cat_d_bases_before_clipping"} = $1;
        } elsif ($line =~ /D total bases written: (\d+)/) {
            $attributes{"cat_d_bases_written"} = $1;
        } elsif ($line =~ /All categories too short: (\d+)\t(\S+)/) {
            $attributes{"all_too_short"} = $1." (".$2."\\\%)";
            $attributes{"n_all_too_short"} = $1;
            $attributes{"pc_all_too_short"} = $2;
        } elsif ($line =~ /All long enough: (\d+)\t(\S+)/) {
            $attributes{"all_long_enough"} = $1." (".$2."\\\%)";
            $attributes{"n_all_long_enough"} = $1;
            $attributes{"pc_all_long_enough"} = $2;
        } elsif ($line =~ /Total usable pairs: (\d+)\t(\S+)/) {
            $attributes{"total_usable"} = $1." (".$2."\\\%)";
            $attributes{"n_total_usable"} = $1;
            $attributes{"pc_total_usable"} = $2;
        } elsif ($line =~ /Number of duplicate pairs: (\d+)\t(\S+)/) {
            $attributes{"number_of_duplicates"} = $1." (".$2."\\\%)";
            $attributes{"n_number_of_duplicates"} = $1;
            $attributes{"pc_number_of_duplicates"} = $2;
        } elsif ($line =~ /Number of pairs containing N: (\d+)\t(\S+)/) {
            $attributes{"number_of_pairs_containing_n"} = $1." (".$2."\\\%)";
            $attributes{"n_number_of_pairs_containing_n"} = $1;
            $attributes{"pc_number_of_pairs_containing_n"} = $2;
        } elsif ($line =~ /Overall GC content: (\S+)/) {
            $attributes{"gc_content"} = $1;
        }
    }

    if (not defined $attributes{"gc_content"}) {
        log_and_screen("\nERROR: Can't read nextclip log - did something go wrong?\n");
        die;
    }
    
    $attributes{"bases_writen_abcd"} = $attributes{"cat_a_bases_written"} + $attributes{"cat_b_bases_written"} + $attributes{"cat_c_bases_written"} + $attributes{"cat_d_bases_written"};
    $attributes{"bases_writen_abc"} = $attributes{"cat_a_bases_written"} + $attributes{"cat_b_bases_written"} + $attributes{"cat_c_bases_written"};
    $attributes{"reference_coverage"} = $attributes{"bases_writen_abc"} / $attributes{"reference_size"};
    
    close(NCFILE);
}

# ----------------------------------------------------------------------
# Function: write_overall_section
# Purpose:  Output overall stats
# ----------------------------------------------------------------------
sub write_overall_section
{
    print $parseable_fh "MinReadSize:".$attributes{"minimum_read_size"}."\n";
    print $parseable_fh "R1NumWithJunctionAdaptor:".$attributes{"n_r1_with_adaptor"}."\n";
    print $parseable_fh "R1PcWithJunctionAdaptor:".$attributes{"pc_r1_with_adaptor"}."\n";
    print $parseable_fh "R1NumWithJunctionAdaptorAndExternalAdaptor:".$attributes{"n_r1_with_adaptor_and_external"}."\n";
    print $parseable_fh "R1PcWithJunctionAdaptorAndExternalAdaptor:".$attributes{"pc_r1_with_adaptor_and_external"}."\n";
    print $parseable_fh "R1NumWithJunctionAdaptorLongEnough:".$attributes{"n_r1_with_adaptor_long"}."\n";
    print $parseable_fh "R1PcWithJunctionAdaptorLongEnough:".$attributes{"pc_r1_with_adaptor_long"}."\n";
    print $parseable_fh "R1NumWithJunctionAdaptorTooShort:".$attributes{"n_r1_with_adaptor_short"}."\n";
    print $parseable_fh "R1PcWithJunctionAdaptorTooShort:".$attributes{"pc_r1_with_adaptor_short"}."\n";
    print $parseable_fh "R1NumWithoutJunctionAdaptor:".$attributes{"n_r1_without_adaptor"}."\n";
    print $parseable_fh "R1PcWithoutJunctionAdaptor:".$attributes{"pc_r1_without_adaptor"}."\n";
    print $parseable_fh "R1NumWithoutJunctionAdaptorButWithExternalAdaptor:".$attributes{"n_r1_without_adaptor_but_external"}."\n";
    print $parseable_fh "R1PcWithoutJunctionAdaptorButWithExternalAdaptor:".$attributes{"pc_r1_without_adaptor_but_external"}."\n";
    print $parseable_fh "R2NumWithJunctionAdaptor:".$attributes{"n_r2_with_adaptor"}."\n";
    print $parseable_fh "R2PcWithJunctionAdaptor:".$attributes{"pc_r2_with_adaptor"}."\n";
    print $parseable_fh "R2NumWithJunctionAdaptorAndExternalAdaptor:".$attributes{"n_r2_with_adaptor_and_external"}."\n";
    print $parseable_fh "R2PcWithJunctionAdaptorAndExternalAdaptor:".$attributes{"pc_r2_with_adaptor_and_external"}."\n";
    print $parseable_fh "R2NumWithJunctionAdaptorLongEnough:".$attributes{"n_r2_with_adaptor_long"}."\n";
    print $parseable_fh "R2PcWithJunctionAdaptorLongEnough:".$attributes{"pc_r2_with_adaptor_long"}."\n";
    print $parseable_fh "R2NumWithJunctionAdaptorTooShort:".$attributes{"n_r2_with_adaptor_short"}."\n";
    print $parseable_fh "R2PcWithJunctionAdaptorTooShort:".$attributes{"pc_r2_with_adaptor_short"}."\n";
    print $parseable_fh "R2NumWithoutJunctionAdaptor:".$attributes{"n_r2_without_adaptor"}."\n";
    print $parseable_fh "R2PcWithoutJunctionAdaptor:".$attributes{"pc_r2_without_adaptor"}."\n";
    print $parseable_fh "R2NumWithoutJunctionAdaptorButWithExternalAdaptor:".$attributes{"n_r2_without_adaptor_but_external"}."\n";
    print $parseable_fh "R2PcWithoutJunctionAdaptorButWithExternalAdaptor:".$attributes{"pc_r2_without_adaptor_but_external"}."\n";
    print $parseable_fh "CatATotalNum:".$attributes{"n_cat_a_total"}."\n";
    print $parseable_fh "CatATotalPc:".$attributes{"pc_cat_a_total"}."\n";
    print $parseable_fh "CatANumLongEnough:".$attributes{"n_cat_a_long"}."\n";
    print $parseable_fh "CatAPcLongEnough:".$attributes{"pc_cat_a_long"}."\n";
    print $parseable_fh "CatANumTooShort:".$attributes{"n_cat_a_short"}."\n";
    print $parseable_fh "CatAPcTooShort:".$attributes{"pc_cat_a_short"}."\n";
    print $parseable_fh "CatANumTrimmedForExternal:".$attributes{"n_cat_a_external"}."\n";
    print $parseable_fh "CatAPcTrimmedForExternal:".$attributes{"pc_cat_a_external"}."\n";
    print $parseable_fh "CatBTotalNum:".$attributes{"n_cat_b_total"}."\n";
    print $parseable_fh "CatBTotalPc:".$attributes{"pc_cat_b_total"}."\n";
    print $parseable_fh "CatBNumLongEnough:".$attributes{"n_cat_b_long"}."\n";
    print $parseable_fh "CatBPcLongEnough:".$attributes{"pc_cat_b_long"}."\n";
    print $parseable_fh "CatBNumTooShort:".$attributes{"n_cat_b_short"}."\n";
    print $parseable_fh "CatBPcTooShort:".$attributes{"pc_cat_b_short"}."\n";
    print $parseable_fh "CatBNumTrimmedForExternal:".$attributes{"n_cat_b_external"}."\n";
    print $parseable_fh "CatBPcTrimmedForExternal:".$attributes{"pc_cat_b_external"}."\n";
    print $parseable_fh "CatCTotalNum:".$attributes{"n_cat_c_total"}."\n";
    print $parseable_fh "CatCTotalPc:".$attributes{"pc_cat_c_total"}."\n";
    print $parseable_fh "CatCNumLongEnough:".$attributes{"n_cat_c_long"}."\n";
    print $parseable_fh "CatCPcLongEnough:".$attributes{"pc_cat_c_long"}."\n";
    print $parseable_fh "CatCNumTooShort:".$attributes{"n_cat_c_short"}."\n";
    print $parseable_fh "CatCPcTooShort:".$attributes{"pc_cat_c_short"}."\n";
    print $parseable_fh "CatCNumTrimmedForExternal:".$attributes{"n_cat_c_external"}."\n";
    print $parseable_fh "CatCPcTrimmedForExternal:".$attributes{"pc_cat_c_external"}."\n";
    print $parseable_fh "CatDTotalNum:".$attributes{"n_cat_d_total"}."\n";
    print $parseable_fh "CatDTotalPc:".$attributes{"pc_cat_d_total"}."\n";
    print $parseable_fh "CatDNumLongEnough:".$attributes{"n_cat_d_long"}."\n";
    print $parseable_fh "CatDPcLongEnough:".$attributes{"pc_cat_d_long"}."\n";
    print $parseable_fh "CatDNumTooShort:".$attributes{"n_cat_d_short"}."\n";
    print $parseable_fh "CatDPcTooShort:".$attributes{"pc_cat_d_short"}."\n";
    print $parseable_fh "CatDNumTrimmedForExternal:".$attributes{"n_cat_d_external"}."\n";
    print $parseable_fh "CatDPcTrimmedForExternal:".$attributes{"pc_cat_d_external"}."\n";
    print $parseable_fh "AllNumTooShort:".$attributes{"n_all_too_short"}."\n";
    print $parseable_fh "AllPcTooShort:".$attributes{"pc_all_too_short"}."\n";
    print $parseable_fh "AllNumLongEnough:".$attributes{"n_all_long_enough"}."\n";
    print $parseable_fh "AllPcLongEnough".$attributes{"pc_all_long_enough"}."\n";
    print $parseable_fh "NumTotalUsable:".$attributes{"n_total_usable"}."\n";
    print $parseable_fh "PcTotalUsable:".$attributes{"pc_total_usable"}."\n";
    print $parseable_fh "NumDuplicates:".$attributes{"n_number_of_duplicates"}."\n";
    print $parseable_fh "PcDuplicates:".$attributes{"pc_number_of_duplicates"}."\n";
    print $parseable_fh "NumPairsWithN:".$attributes{"n_number_of_pairs_containing_n"}."\n";
    print $parseable_fh "PcPairsWithN:".$attributes{"pc_number_of_pairs_containing_n"}."\n";
    print $parseable_fh "GCContent:".$attributes{"gc_content"}."\n";
    
    print $latex_fh "\\subsection*{Overall}\n";
    print $latex_fh "\\vspace{-3mm}\n";
    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{c c c}\n";
    print $latex_fh "{\\bf R1} & & {\\bf R2} \\\\\n";
    print $latex_fh $attributes{"number_of_pairs"}." & Number of reads & ".$attributes{"number_of_pairs"}." \\\\\n";
    print $latex_fh $attributes{"r1_with_adaptor"}." & With junction adaptor & ".$attributes{"r2_with_adaptor"}." \\\\\n";
    print $latex_fh $attributes{"r1_with_adaptor_long"}." & of which long enough (\$\\ge\$ ".$attributes{"minimum_read_size"}.") & ".$attributes{"r2_with_adaptor_long"}." \\\\\n";
    print $latex_fh $attributes{"r1_with_adaptor_short"}." & and too short (\$<\$ ".$attributes{"minimum_read_size"}.") & ".$attributes{"r2_with_adaptor_short"}." \\\\\n";
    print $latex_fh $attributes{"r1_without_adaptor"}." & Without junction adaptor & ".$attributes{"r2_without_adaptor"}." \\\\\n";
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";
    print $latex_fh "\\vspace{-5mm}\n";
    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{l c c c c}\n";
    print $latex_fh "{\\bf Category} & {\\bf Number of pairs} & {\\bf Too short (\$<\$ ".$attributes{"minimum_read_size"}.")} & {\\bf Long enough (\$\\ge\$ ".$attributes{"minimum_read_size"}.")} & {\\bf Bases written} \\\\\n";
    print $latex_fh "Adaptor in R1 and R2 (A) & ".$attributes{"cat_a_total"}." & ".$attributes{"cat_a_short"}." & ".$attributes{"cat_a_long"}." & ".$attributes{"cat_a_bases_written"}." \\\\\n";
    print $latex_fh "Adaptor in R2 only (B) & ".$attributes{"cat_b_total"}." & ".$attributes{"cat_b_short"}." & ".$attributes{"cat_b_long"}." & ".$attributes{"cat_b_bases_written"}." \\\\\n";
    print $latex_fh "Adaptor in R1 only (C) & ".$attributes{"cat_c_total"}." & ".$attributes{"cat_c_short"}." & ".$attributes{"cat_c_long"}." & ".$attributes{"cat_c_bases_written"}." \\\\\n";
    print $latex_fh "Adaptor in neither (D) & ".$attributes{"cat_d_total"}." & ".$attributes{"cat_d_short"}." & ".$attributes{"cat_d_long"}." & ".$attributes{"cat_d_bases_written"}." \\\\\n";
    print $latex_fh "All categories & ".$attributes{"number_of_pairs"}." (100\\\%) & ".$attributes{"all_too_short"}." & ".$attributes{"all_long_enough"}." & ".$attributes{"bases_writen_abcd"}." \\\\\n";
    print $latex_fh "Total usable (A,B,C) & & & ".$attributes{"total_usable"}." & ".$attributes{"bases_writen_abc"};
    printf $latex_fh " (%.1f x) \\\\\n", $attributes{"reference_coverage"};
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";
}

# ----------------------------------------------------------------------
# Function: write_gc_content_section
# Purpose:  Output GC stats
# ----------------------------------------------------------------------
sub write_gc_content_section
{
    print $latex_fh "\\subsection*{GC content}\n";
 
    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{l c}\n";
    print $latex_fh "Overall GC content & ".$attributes{"gc_content"}."\\\% \\\\\n";
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";

    print $latex_fh "\\vspace{-3mm}\n";
    print $latex_fh "\\begin{figure}[h!]\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\begin{subfigure}{.33\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_R1_gc.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize R1}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\begin{subfigure}{.33\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_R2_gc.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize R2}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\end{figure}\n";
}

# ----------------------------------------------------------------------
# Function: write_shortest_pair_section
# Purpose:  Output shortest pair stats
# ----------------------------------------------------------------------
sub write_shortest_pair_section
{
    print $latex_fh "\\subsection*{Shortest pair length}\n";
    print $latex_fh "\\vspace{-3mm}\n";
    print $latex_fh "\\begin{figure}[h!]\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\begin{subfigure}{.22\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_cumulative_pairs_A.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Category A pairs}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\begin{subfigure}{.22\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_cumulative_pairs_B.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Category B pairs}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\begin{subfigure}{.22\\textwidth}\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_cumulative_pairs_C.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Category C pairs}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\end{figure}\n";
}

# ----------------------------------------------------------------------
# Function: write_clipped_read_length_section
# Purpose:  Output clipped read length stats
# ----------------------------------------------------------------------
sub write_clipped_read_length_section
{
    print $latex_fh "\\subsection*{Clipped read lengths}\n";
    print $latex_fh "\\vspace{-3mm}\n";
    print $latex_fh "\\begin{figure}[h!]\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\begin{subfigure}{.22\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_lengths_A_R1.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Category A R1 lengths}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\begin{subfigure}{.22\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_lengths_A_R2.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Category A R2 lengths}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\begin{subfigure}{.22\\textwidth}\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_lengths_B_R2.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Category B R2 lengths}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\begin{subfigure}{.22\\textwidth}\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.9\\linewidth]{".$library_dir."/graphs/".$library_name."_lengths_C_R1.pdf}\n";
    print $latex_fh "\\caption*{\\footnotesize Category C R1 lengths}\n";
    print $latex_fh "\\end{subfigure}\n";
    print $latex_fh "\\end{figure}\n";
}

# ----------------------------------------------------------------------
# Function: write_duplication_section
# Purpose:  Output PCR duplication stats
# ----------------------------------------------------------------------
sub write_duplication_section
{
    print $latex_fh "\\subsection*{Duplication}\n";
    print $latex_fh "\\vspace{-3mm}\n";
#    
    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{l c}\n";
    print $latex_fh "Number of read pairs in library & ".$attributes{"number_of_pairs"}." \\\\\n";
    print $latex_fh "PCR duplicates & ".$attributes{"number_of_duplicates"}." \\\\\n";
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";
#    
    print $latex_fh "\\vspace{-5mm}\n";
    print $latex_fh "\\begin{figure}[h!]\n";
    print $latex_fh "\\centering\n";
    print $latex_fh "\\includegraphics[width=.35\\linewidth]{".$library_dir."/graphs/".$library_name."_duplicates.pdf}\n";
    print $latex_fh "\\end{figure}\n";
}

# ----------------------------------------------------------------------
# Function: write_external_adaptor_section
# Purpose:  Output external adaptor stats
# ----------------------------------------------------------------------
sub write_external_adaptor_section
{
    print $latex_fh "\\subsection*{External adaptor}\n";
    print $latex_fh "\\vspace{-3mm}\n";
    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{c c c}\n";
    print $latex_fh "{\\bf R1} & & {\\bf R2} \\\\\n";
    print $latex_fh $attributes{"number_of_pairs"}." & Number of reads & ".$attributes{"number_of_pairs"}." \\\\\n";
    print $latex_fh $attributes{"r1_with_adaptor_and_external"}." & With junction adaptor and external adaptor & ".$attributes{"r2_with_adaptor_and_external"}." \\\\\n";
    print $latex_fh $attributes{"r1_without_adaptor_but_external"}." & Without junction adaptor, but with external adaptor & ".$attributes{"r2_without_adaptor_but_external"}." \\\\\n";
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";
    print $latex_fh "\\vspace{-5mm}\n";
    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{l c}\n";
    print $latex_fh "Category A pairs also trimmed for external adaptor & ".$attributes{"cat_a_external"}." \\\\\n";
    print $latex_fh "Category B pairs also trimmed for external adaptor & ".$attributes{"cat_b_external"}." \\\\\n";
    print $latex_fh "Category C pairs also trimmed for external adaptor & ".$attributes{"cat_c_external"}." \\\\\n";
    print $latex_fh "Category D pairs also trimmed for external adaptor & ".$attributes{"cat_d_external"}." \\\\\n";
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";
}

# ----------------------------------------------------------------------
# Function: write_ambiguous_bases_section
# Purpose:  Output ambiguous base stats
# ----------------------------------------------------------------------
sub write_ambiguous_bases_section
{
    print $latex_fh "\\subsection*{Ambiguous bases}\n";
    print $latex_fh "\\vspace{-3mm}\n";
    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{l c}\n";
    print $latex_fh "Number of pairs containing Ns & ".$attributes{"number_of_pairs_containing_n"}." \\\\\n";
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";
}

# ----------------------------------------------------------------------
# Function: write_notes_section
# Purpose:  Output notes
# ----------------------------------------------------------------------
sub write_notes_section
{
    print $latex_fh "\\subsection*{Notes}\n";
    print $latex_fh "\\vspace{-3mm}\n";
    print $latex_fh "\\begin{table}[h!]\n";
    print $latex_fh "{\\footnotesize\n";
    print $latex_fh "\\fontsize{9pt}{11pt}\\selectfont\n";
    print $latex_fh "\\begin{tabular}{l c}\n";
    print $latex_fh "Minimum contig size for alignment to reference & ".$att{"ref_min"}." \\\\\n";
    print $latex_fh "Maximum allowed MP insert & ".$att{"mp_limit"}." \\\\\n";
    print $latex_fh "Maximum allowed PE insert & ".$att{"pe_limit"}." \\\\\n";
    print $latex_fh "Maximum allowed tandem insert & ".$att{"tandem_limit"}." \\\\\n";
    print $latex_fh "Reference size & ".$attributes{"reference_size"}." bp \\\\\n";
    print $latex_fh "\\end{tabular}\n";
    print $latex_fh "}\n";
    print $latex_fh "\\end{table}\n";
}

# ----------------------------------------------------------------------
# Function: load_reference_index
# Purpose:  Find size of reference
# ----------------------------------------------------------------------
sub load_reference_index
{
    my $reference_index = $reference.".nextclip";
    my $total_size = 0;
    
    log_and_screen("Opening index file $reference_index\n");
    log_and_screen("ERROR: nextclip_sam_parse.pl: Can't find $reference_index - have you run nextclip_index_reference.pl as per the manual?\n") unless (-e $reference_index);
    
    open(INDEXFILE, $reference_index) or die "Can't open reference index $reference_index\n";
    while(<INDEXFILE>) {
        chomp(my $line = $_);
        my @arr = split(/\t/, $line);
        $total_size += $arr[1];
    }
    close(INDEXFILE);

    $attributes{"reference_size"} = $total_size;
}

