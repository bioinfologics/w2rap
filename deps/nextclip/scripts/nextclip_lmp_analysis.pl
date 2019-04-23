#!/usr/bin/perl -w

# Script:  nextclip_lmp_analysis.pl
# Purpose: NextClip Long Mate Pair analysis pipeline
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

use warnings;
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;

# ============= CAN BE HARDCODED PATHS IF NECESSARY FOR YOUR ENVIRONMENT =============
my $script_dir;
my $nextclip_tool="nextclip";
# ====================================================================================

# Command line options
my $configure_file;
my $bwa_threads=1;
my $min_map_q=10;
my $scheduler="NONE";
my $help;

# Things in the configure file
my $library_name;
my $organism;
my $output_dir;
my $read_one;
my $read_two;
my $reference;
my $minimum_contig_alignment_size;
my $number_of_pairs;
my $min_length=25;
my $trim_ends=0;
my $nextclip_options;
my $start_stage=0;

# We'll work these out later
my $logdir;
my $readsdir;
my $bwadir;
my $analysisdir;
my $read_one_clipped;
my $read_two_clipped;
my $read_one_noadapt;
my $read_two_noadapt;
my $read_length;
my $machine;
my $log_filename;
my %job_id_table;
my $queue;

# Get command line options
&GetOptions(
'config:s'     => \$configure_file,
'bwathreads:i' => \$bwa_threads,
'minmapq:i'    => \$min_map_q,
'nextclip:s'   => \$nextclip_tool,
'scriptdir:s'  => \$script_dir,
'scheduler:s'  => \$scheduler,
'queue:s'      => \$queue,
'stage:i'      => \$start_stage,
'help|h'       => \$help
);

print "\nNextClip LMP Analysis v1.3\n\n";

# Help message
if ($help) {
    print "Pipeline for alignment and analysis of LMP libraries with NextClip\n\n";
    print "Syntax: nextclip_lmp_analysis.pl -config <file> [options]\n\n";
    print "Options\n";
    print "    -config <file> - Specify name of configuration file\n";
    print "    -scheduler <name> - Specify job scheduler to use (default NONE)\n";
    print "    -bwathreads <int> - Specify number of threads to use for BWA (default 1)\n";
    print "    -minmapq <int> - Specify minimum mapping quality (default 10)\n";
    print "    -queue <name> - Optionally specify scheduler queue name (default queue used if ommitted)\n";
    print "    -stage <int> - Start at pipeline stage X, where X is:\n";
    print "                      1 = clipping\n";
    print "                      2 = alignment\n";
    print "                      3 = parsing alignment\n";
    print "                      4 = graph plotting\n";
    print "                      5 = report generation\n";
    print "\n";
    exit(0);
}

if (defined $scheduler) {
    $scheduler = uc($scheduler);
}

# Check we know everything
die "Error: You must specify a config file\n" if not defined $configure_file;
die "Error: Scheduler must be either LSF, PBS or none\n" if (($scheduler ne "NONE") && ($scheduler ne "LSF") && ($scheduler ne "PBS"));
die "Error: BWA threads must be between 1 and 32\n" if (($bwa_threads <1) || ($bwa_threads > 32));
die "Error: minmapq must be between 0 and 255\n" if (($min_map_q < 0) || ($min_map_q > 255));

# Check script dir
if (! defined $script_dir) {
    $script_dir=dirname($0);
}

# Make sure there is a / at the end of script_dir
if ($script_dir !~ /\/$/) {
    $script_dir = $script_dir."/";
}

# Main stages of pipeline
print "Scripts dir is $script_dir\n";

read_config_file($configure_file);
check_software();
find_read_length_and_machine($read_one);

if ($start_stage <= 1) {
    run_clipping();
}

if ($start_stage <= 2) {
    run_alignment();
}

if ($start_stage <= 3) {
    parse_alignment();
}

if ($start_stage <= 4) {
    plot_graphs();
}

if ($start_stage <= 5) {
    make_report();
}

log_print("\n");
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

# Check software versions
sub check_software
{
    my @nextclip_output = readpipe($nextclip_tool." -h 2>&1");
    my @bwa_output = readpipe("bwa 2>&1");
    my $r_output = readpipe("Rscript --version 2>&1");
    my @tex_output = readpipe("pdflatex --version 2>&1");
    my $nextclip_found = 0;
    my $bwa_found = 0;
    my $r_found = 0;
    my $tex_found = 0;
    my $line;

    die "Error: can't find scripts directory (path specified $script_dir)!" unless (-d $script_dir);

    foreach $line (@nextclip_output) {
        if ($line =~ /^NextClip/) {
            $nextclip_found = 1;
            chomp($line);
            log_and_screen $line, "\n";
        }
    }

    die "Error: can't find NextClip tool!" unless ($nextclip_found == 1);
    
    foreach $line (@bwa_output) {
        if ($line =~ /^Version: (\S+)/) {
            my $version = $1;
            $bwa_found = 1;
            if ($version =~ /(\d+).(\d+).(\d+)-(\S+)$/) {
                my $old = 0;
                log_and_screen "BWA version ",$version,"\n";
                if ((($1 == 0) && ($2 < 6)) ||
                    (($1 == 0) && ($2 == 6) && ($3 < 1))) {
                        print"\nWARNING: BWA version appears to be too old (needs to be 0.6.1 or greater).\n";
                        print"         Will try continuing, but may not work.\n\n"
                }
            }
        }
    }

    if (defined $r_output) {
        if ($r_output =~ /version (\d+).(\d+).(\d+)/) {
            log_and_screen $r_output;
            $r_found = 1;
            if (($1 < 2) ||
                (($1 == 2) && ($2 < 12)) ||
                (($1 == 2) && ($2 == 12) && ($3 < 2))) {
                    print"\nWARNING: R version appears to be too old (needs to be 2.12.2 or greater).\n";
                    print"         Will try continuing, but may not work.\n\n"
                }
        }
    }
    
    foreach $line (@tex_output) {
        if ($line =~ /^pdfTeX/) {
            $tex_found = 1;
            log_and_screen $line, "\n";
        }
    }
    

    if ($nextclip_found ==  0) {
        print"\nWARNING: Could not find NextClip. Will try continuing, but probably won't work!\n\n";
    }
    
    if ($bwa_found ==  0) {
        print"\nWARNING: Could not find BWA. Will try continuing, but probably won't work!\n\n";
    }
    
    if ($r_found ==  0) {
        print"\nWARNING: Could not find Rscript. Will try continuing, but probably won't work!\n\n";
    }
    
    if ($tex_found ==  0) {
        print"\nWARNING: Could not find pdflatex. Will try continuing, but probably won't work!\n\n";
    }
}

# Read config file, check files exist, make dirs as appropriate
sub read_config_file
{
    my $filename = $_[0];
    open(CONFIGFILE, $filename) or die "Can't open $filename\n";
    while (<CONFIGFILE>) {
        chomp(my $line = $_);
        my @arr=split(':', $line);
        if (defined $arr[0]) {
            if ($arr[0] eq "library_name") {
                $library_name = $arr[1];
            } elsif ($arr[0] eq "organism") {
                $organism = $arr[1];
            } elsif ($arr[0] eq "output_dir") {
                $output_dir = $arr[1];
            } elsif ($arr[0] eq "read_one") {
                $read_one = $arr[1];
            } elsif ($arr[0] eq "read_two") {
                $read_two = $arr[1];
            } elsif ($arr[0] eq "reference") {
                $reference = $arr[1];
            } elsif ($arr[0] eq "minimum_contig_alignment_size") {
                $minimum_contig_alignment_size = $arr[1];
            } elsif ($arr[0] eq "number_of_pairs") {
                $number_of_pairs = $arr[1];
            } elsif ($arr[0] eq "min_length") {
                $min_length = $arr[1];
            } elsif ($arr[0] eq "trim_ends") {
                $trim_ends = $arr[1];
            } elsif ($arr[0] eq "nextclip_options") {
                $nextclip_options = $arr[1];
            }
        }
    }
    close(CONFIGFILE);
    
    die "Error: Config file contained no library_name\n" if not defined $library_name;
    die "Error: Config file contained no organism\n" if not defined $organism;
    die "Error: Config file contained no output_dir\n" if not defined $output_dir;
    die "Error: Config file contained no read_one\n" if not defined $read_one;
    die "Error: Config file contained no read_two\n" if not defined $read_two;
    die "Error: Config file contained no reference\n" if not defined $reference;
    die "Error: Config file contained no minimum_contig_alignment_size\n" if not defined $minimum_contig_alignment_size;

    die "Error: read 1 and read 2 are the same!\n" if ($read_one eq $read_two);
    
    die "Error: read 1 file does not exist\n" unless (-e $read_one);
    die "Error: read 2 file does not exist\n" unless (-e $read_two);
    
    unless ( -d $output_dir) {
        print "Creating directory $output_dir\n";
        mkdir($output_dir);
    }
    
    die "Error: couldn't create directory $output_dir" unless (-d $output_dir);

    my @dirs = qw(bwa reads logs graphs analysis latex);
    foreach my $dir (@dirs) {
        my $newdir = $output_dir."/".$dir;
        unless (-d $newdir) {
            print "Creating directory $newdir\n";
            mkdir ($newdir);
            die "Error: couldn't create directory $newdir" unless (-d $newdir);
        }
    }
    
    $logdir=$output_dir."/logs";
    $readsdir=$output_dir."/reads";
    $bwadir=$output_dir."/bwa";
    $analysisdir=$output_dir."/analysis";
    $read_one_clipped = $readsdir."/R1_clipped.fastq";
    $read_two_clipped = $readsdir."/R2_clipped.fastq";
    $read_one_noadapt = $readsdir."/R1_noadapt.fastq";
    $read_two_noadapt = $readsdir."/R2_noadapt.fastq";
    $log_filename = $logdir."/nextclip_analysis.log";

    # Initialise log
    open(my $log_fh, ">".$log_filename);
    log_print "NextClip LMP analysis\n\n";
    close($log_fh);
    
    # Output variables
    print "\n";
    log_and_screen "    Library name: $library_name\n";
    log_and_screen "        Organism: $organism\n";
    log_and_screen "      Output dir: $output_dir\n";
    log_and_screen "          Read 1: $read_one\n";
    log_and_screen "          Read 2: $read_two\n";
    log_and_screen "       Reference: $reference\n";
    log_and_screen "  Read 1 clipped: $read_one_clipped\n";
    log_and_screen "  Read 2 clipped: $read_two_clipped\n";
    log_and_screen "  Read 1 noadapt: $read_one_noadapt\n";
    log_and_screen "  Read 2 noadapt: $read_two_noadapt\n";
    log_and_screen "        Log file: $log_filename\n";
    log_and_screen "      Min length: $min_length\n";
    log_and_screen "       Trim ends: $trim_ends\n";
    if (defined $nextclip_options) {
        log_and_screen "Nextclip options: $nextclip_options\n";
    }
    log_and_screen "\n";
}

# Find read length and machine name
sub find_read_length_and_machine
{
    my $filename = $_[0];
    open(READFILE, $filename) or die "Can't open $filename\n";
    chomp(my $header = <READFILE>);
    chomp(my $read = <READFILE>);
    close(READFILE);
    
    $read_length = length($read);
    my @arr = split(':', $header);
    $machine = $arr[0];
    log_and_screen "     Read length: $read_length\n";
    log_and_screen "    Machine name: $machine\n";

    if (not defined $number_of_pairs) {
        if ($start_stage <= 1) {
            print "Counting number of input pairs...\n\n";
            my $result = readpipe("wc -l ".$read_one); 
            my @arr = split(' ', $result);
            my $lines = $arr[0];
            
            if ($lines > 0) {
                $number_of_pairs = $lines / 4;
            } else {
                log_and_screen "Something went wrong with line count - defaulting to 100000000\n";
                $number_of_pairs = 100000000;
            }
        } else {
            $number_of_pairs = 0;
        }
    }

    log_and_screen " Number of pairs: $number_of_pairs\n\n";
}

# Handle PBS
#
# job_id_table handles the mapping between PBS IDs and the job IDs used internally in this program.
sub submit_pbs
{
    my $command = $_[0];
    my $this_id = $_[1];
    my $log = $_[2];
    my $previous_id = $_[3];
    my $threads = $_[4];
    my $memory = $_[5];
    my $cwd = getcwd();
    my $system_command = "echo \"cd ".$cwd." ; ".$command."\" | qsub -V -l ncpus=".$threads.",mem=".$memory."mb,walltime=4800:00:00 -j oe -o ".$log; #." -N ".$this_id;

    if (defined $queue) {
        $system_command = $system_command." -q ".$queue;
    }
    
    if ($previous_id ne "") {
        if ($previous_id =~ /\*$/) {
            my $count=0;
            $previous_id =~ s/\*$//;
            foreach my $id ( keys %job_id_table )
            {
                if ($id =~ /^$previous_id/) {
                    if ($count == 0) {
                        $system_command = $system_command." -W depend=afterany"
                    }
                    $system_command = $system_command.":".$job_id_table{$id};
                    $count++;
                }
            }
        } else {
            die "Can't find previous dependency $previous_id for PBS.\n" if (not defined $job_id_table{$previous_id});
            $system_command = $system_command." -W depend=afterany:".$job_id_table{$previous_id};
        }
    }
    
    log_print "PBS command: $system_command\n\n";
    my $newid = readpipe($system_command);
    chomp ($newid);
    $job_id_table{$this_id} = $newid;
}

# Handle LSF
sub submit_lsf
{
    my $command = $_[0];
    my $this_id = $_[1];
    my $log = $_[2];
    my $previous_id = $_[3];
    my $threads = $_[4];
    my $memory = $_[5];
    my $system_command;

    $system_command="bsub -oo ".$log." -J ".$this_id;
    
    if (defined $queue) {
        $system_command = $system_command." -q ".$queue;
    }
    
    if ($previous_id ne "") {
        $system_command = $system_command." -w \"ended(".$previous_id.")\"";
    }
    if ($threads > 1) {
        $system_command = $system_command." -R \"rusage[mem=".$memory."] span[hosts=1]\" -n ".$threads;
    } else {
        $system_command = $system_command." -R \"rusage[mem=".$memory."]\"";
    }
    $system_command=$system_command." \"".$command."\"";
    log_print "LSF command: $system_command\n\n";
    system($system_command);
}

# Submit a job
# To add a different job scheduler:
# 1. Add an if statement below
# 2. Code a submit_job_XYZ function, where XYZ is the scheduler.
# 3. Amend the die statement near the start of the program that checks for valid scheduler values
sub submit_job
{
    my $command = $_[0];
    my $this_id = $_[1];
    my $log = $_[2];
    my $previous_id = $_[3];
    my $threads = $_[4];
    my $memory = $_[5];
    
    if ($scheduler eq "LSF") {
        submit_lsf($command, $this_id, $log, $previous_id, $threads, $memory);
    } elsif ($scheduler eq "PBS") {
        submit_pbs($command, $this_id, $log, $previous_id, $threads, $memory);
    } elsif ($scheduler eq "NONE") {
        my $system_command;

        # Minor bodge to ensure nextclip output ends up in a log when not running with scheduler
        if (($this_id =~ /_sam_/) ||
            ($this_id =~ /_sai_/) ||
            ($this_id =~ /_latex$/)) {
            $system_command = $command;
        } elsif ($this_id =~/clip$/) {
            $system_command = $command." 2>&1 | tee ".$log;
        } else {
            $system_command = $command." > ".$log;
        }

        log_and_screen "Executing $system_command\n\n";
        system($system_command);
    } else {
        die "Unknown job scheduler $scheduler\n";
    }
    
}

# Run NextClip
sub run_clipping
{
    my $num_pairs = $number_of_pairs * 2;
    my $command=$nextclip_tool." --input_one ".$read_one." --input_two ".$read_two." --output_prefix ".$readsdir."/".$library_name." --log ".$logdir."/nextclip_alignment.log --min_length ".$min_length." --number_of_reads ".$num_pairs." --trim_ends ".$trim_ends;
    my $logfile=$logdir."/nextclip.log";
    my $previous="";
    my $job_id=$library_name."clip";
    my $memory=4000;

    if (defined $nextclip_options) {
        $command = $command." ".$nextclip_options;
    }

    # Find memory requirements
    my @output = readpipe($nextclip_tool." --number_of_reads ".$num_pairs." --memory_requirements");
    foreach my $line (@output) {
        if ($line =~ /Memory required: (\d+) MB/) {
            $memory = $1;
            log_and_screen "Queried required memory for NextClip, received: $memory\n";
        }
    }
    
    log_and_screen "Memory required for NextClip: $memory\n";
    submit_job($command, $job_id, $logfile, $previous, 1, $memory);
}

# Run BWA
sub run_alignment
{
    my @types = qw(A B C D);
    foreach my $type (@types) {
        for (my $read=1; $read<=2; $read++) {
            my $fastqfile=$readsdir."/".$library_name."_".$type."_R".$read.".fastq";
            my $saifile=$bwadir."/".$library_name."_".$type."_R".$read.".sai";
            my $samfile=$bwadir."/".$library_name."_".$type."_R".$read.".sam";
    
            if ($read_length > 101) {
                my $command="bwa bwasw -M -t ".$bwa_threads." ".$reference." ".$fastqfile." > ".$samfile;
                my $logfile=$logdir."/sam_".$type."_R".$read.".log";
                my $previous=$library_name."clip";
                my $job_id=$library_name."_sam_".$type."_R".$read;
                submit_job($command, $job_id, $logfile, $previous, $bwa_threads, 4000);
            } else {
                my $command="bwa aln -t ".$bwa_threads." ".$reference." ".$fastqfile." > ".$saifile;
                my $logfile=$logdir."/sai_".$type."_R".$read.".log";
                my $previous=$library_name."clip";
                my $job_id=$library_name."_sai_".$type."_R".$read;
                submit_job($command, $job_id, $logfile, $previous, $bwa_threads, 4000);
                $command="bwa samse ".$reference." ".$saifile." ".$fastqfile." > ".$samfile;
                $logfile=$logdir."/sam_".$type."_R".$read.".log";
                $previous=$job_id;
                $job_id=$library_name."_sam_".$type."_R".$read;
                submit_job($command, $job_id, $logfile, $previous, $bwa_threads, 4000);
            }
        }
    }
}

# Parse BWA output
sub parse_alignment
{
    my @types = qw(A B C D);
    foreach my $type (@types) {
        my $sam_one=$bwadir."/".$library_name."_".$type."_R1.sam";
        my $sam_two=$bwadir."/".$library_name."_".$type."_R2.sam";
        my $out_prefix=$analysisdir."/".$library_name."_".$type;
        my $command="perl ".$script_dir."nextclip_sam_parse.pl -one ".$sam_one." -two ".$sam_two." -out ".$out_prefix." -reference ".$reference." -refminsize ".$minimum_contig_alignment_size." -minmapq ".$min_map_q." -log ".$log_filename;
        my $logfile=$logdir."/parse_".$type.".log";
        my $previous=$library_name."_sam_".$type."*";
        my $job_id=$library_name."_parse_".$type;
        submit_job($command, $job_id, $logfile, $previous, 1, 4000);
    }
}

# Plot graphs
sub plot_graphs
{
    my $command=$script_dir."nextclip_plot_graphs.sh ".$output_dir." ".$library_name." ".$log_filename." \"".$script_dir."\"";
    my $logfile=$logdir."/plot_graphs.log";
    my $previous=$library_name."_parse*";
    my $job_id=$library_name."_plot_graphs";
    submit_job($command, $job_id, $logfile, $previous, 1, 4000);
}

# Make LaTeX report
sub make_report
{
    my $command=$script_dir."nextclip_make_report.pl -libdir ".$output_dir." -libname ".$library_name." -reference '".$reference."' -organism '".$organism."'"." -log ".$log_filename." -nextclip ".$logdir."/nextclip.log";
    my $logfile=$logdir."/makereport.log";
    my $previous=$library_name."_plot_graphs";
    my $job_id=$library_name."_latex";
    submit_job($command, $job_id, $logfile, $previous, 1, 4000);
}
