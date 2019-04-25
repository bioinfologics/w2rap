/*----------------------------------------------------------------------*
 * File:    nextclip.c                                                  *
 * Purpose: Analysis and read preparation for Nextera LMP               *
 * Author:  Richard Leggett                                             *
 *          The Genome Analysis Centre (TGAC), Norwich, UK              *
 *          richard.leggett@tgac.ac.uk    								*
 *----------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "binary_kmer.h"
#include "element.h"
#include "hash_table.h"

/* Experiment that I decided against using
 * #define USE_MULTIPLE_HASHES
 * #define NUMBER_OF_HASHES 3
 * Who knows, may use in the future!
 */

/*----------------------------------------------------------------------*
 * Constants
 *----------------------------------------------------------------------*/
#define NEXTCLIP_VERSION "1.3.2"
#define MAX_PATH_LENGTH 1024
#define NUMBER_OF_CATEGORIES 5
#define SEPARATE_KMER_SIZE 11
#define TOTAL_KMER_SIZE (4 * SEPARATE_KMER_SIZE)
#define MAX_DUPLICATES 1000
#define FIRST_KMER_OFFSET 20

// The kmer-based PCR duplication assessment won't work so well with very small reads, so have set this limit
#define MINIMUM_INPUT_READ_SIZE 64

/*----------------------------------------------------------------------*
 * Structures
 *----------------------------------------------------------------------*/
typedef struct {
    char read_header[1024];
    char read[MAX_READ_LENGTH];
    char quality_header[1024];
    char qualities[MAX_READ_LENGTH];
    int read_size;
    int trim_at_base;
    boolean trimmed_for_junction_adaptor;
    boolean trimmed_for_external_adaptor;
} FastQRead;

typedef struct {
    int read_size;
    int score;
    int total_matches;
    double total_identity;
    int total_alignment_length;
    int position;
    int matches[2];
    int mismatches[2];
    int alignment_length[2];
    double identity[2];
    int read_start;
    int read_end;
    int query_start;
    int query_end;
    int accepted;
} JunctionAdaptorAlignment;

typedef struct {
    char adaptor[100];
    int read_size;
    int score;
    int position;
    int matches;
    int mismatches;
    int alignment_length;
    double identity;
    int read_start;
    int read_end;
    int query_start;
    int query_end;
    int accepted;
} GenericAdaptorAlignment;

typedef struct {
    int read_length;
    FILE* input_fp[2];
    FILE* output_fp[NUMBER_OF_CATEGORIES][2];
    FILE* log_fp;
    FILE* duplicates_fp;
    char input_filenames[2][MAX_PATH_LENGTH];
    char output_prefix[MAX_PATH_LENGTH];
    char output_filenames[NUMBER_OF_CATEGORIES][MAX_PATH_LENGTH];
    char log_filename[MAX_PATH_LENGTH];
    char duplicates_log_filename[MAX_PATH_LENGTH];
    int num_read_pairs;
    int count_adaptor_found[2];
    int count_adaptor_and_external_found[2];
    int count_no_adaptor[2];
    int count_external_only_found[2];
    int count_long_enough[2];
    int count_too_short[2];
    double percent_adaptor_found[2];
    double percent_no_adaptor[2];
    double percent_long_enough[2];
    double percent_too_short[2];
    double percent_adaptor_and_external_found[2];
    double percent_external_only_found[2];
    int count_by_category[NUMBER_OF_CATEGORIES];
    int count_by_category_long_enough[NUMBER_OF_CATEGORIES];
    int count_by_category_too_short[NUMBER_OF_CATEGORIES];
    int count_by_category_relaxed_hit[NUMBER_OF_CATEGORIES];
    int count_by_category_external_clipped[NUMBER_OF_CATEGORIES];
    double percent_by_category[NUMBER_OF_CATEGORIES];
    double percent_by_category_long_enough[NUMBER_OF_CATEGORIES];
    double percent_by_category_too_short[NUMBER_OF_CATEGORIES];
    double percent_by_category_relaxed_hit[NUMBER_OF_CATEGORIES];
    double percent_by_category_external_clipped[NUMBER_OF_CATEGORIES];
    int total_too_short;
    double percent_total_too_short;
    int total_long_enough;
    int total_usable;
    double percent_usable;
    double percent_total_long_enough;
    int read_length_counts[NUMBER_OF_CATEGORIES][2][MAX_READ_LENGTH];
    int read_pair_length_counts[NUMBER_OF_CATEGORIES][MAX_READ_LENGTH];
    int n_duplicates;
    int n_invalid_for_duplicate;
    double percent_duplicates;
    int duplicates_not_written;
    double percent_duplicates_not_written;
    long int at_bases;
    long int gc_bases;
    int gc_content[2][101];
    double percent_gc;
    long int pairs_containing_n;
    double percent_pairs_containing_n;
    long int bases_before_clipping[NUMBER_OF_CATEGORIES];
    long int bases_written[NUMBER_OF_CATEGORIES];
} MPStats;

/*----------------------------------------------------------------------*
 * Globals
 *----------------------------------------------------------------------*/
char single_junction_adaptor[128] = "CTGTCTCTTATACACATCT";
char duplicate_junction_adaptor[128] = "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG";
char* external_adaptors[2] = {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"};
int minimum_read_size = 25;
int strict_double_match = 34;
int strict_single_match = 18;
int relaxed_double_match = 32;
int relaxed_single_match = 17;
int discarded = 0;
int use_category_e = 0;
int remove_duplicates = 0;
int num_categories = NUMBER_OF_CATEGORIES;
int trim_ends = 19;
int approximate_reads = 20000000;
int output_memory_requirements = false;
int duplicate_only_mode = false;

/*
 * Single hash option algorithm
 *
 * We take a kmer from the start of each read and one from the middle, eg.
 *
 *   R1: TGACTGACTGACTGACTGACTGACTGACTG     R2: TGACTGACTGACTGACTGACTGACTGACTG
 * kmer: AAA           BBB                      aaa           bbb
 *
 * Then we make a hash AAABBBaaabbb.
 *
 *
 * Multiple hash option algorithm
 *
 * We take a kmer from each read, from a third of the way and from two thirds of the way through, eg.
 *
 *   R1: TGACTGACTGACTGACTGACTGACTGACTG     R2: TGACTGACTGACTGACTGACTGACTGACTG
 * kmer: AAA       BBB       CCC                aaa       bbb       ccc
 *
 * Then we make three kmer hashes:
 * kmer_hashes[0] = AAABBBaaabbb
 * kmer_hashes[1] = AAACCCaaaccc
 * kmer_hashes[2] = BBBCCCbbbccc
 *
 * Then, for something to be declared a duplicate, we ask for a match to any of the three hash tables.
 * This allows one of the three kmers to not match - eg. because of miscalled base.
 */

#ifdef USE_MULTIPLE_HASHES
HashTable* kmer_hashes[NUMBER_OF_HASHES];
int kmer_offsets[NUMBER_OF_HASHES][2];
#else
HashTable* duplicate_hash = NULL;
#endif

/*----------------------------------------------------------------------*
 * Function:   initialise_stats
 * Purpose:    Initialise MPStats structure
 * Parameters: stats -> an MPStats structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void initialise_stats(MPStats* stats)
{
    int i, j;
    
    for (i=0; i<2; i++) {
        stats->input_filenames[i][0] = 0;
        stats->count_adaptor_found[i] = 0;
        stats->count_no_adaptor[i] = 0;
        stats->count_too_short[i] = 0;
        stats->count_long_enough[i] = 0;
        stats->input_fp[i] = 0;
        stats->count_adaptor_and_external_found[i] = 0;
        stats->count_external_only_found[i] = 0;
    }

    for (i=0; i<NUMBER_OF_CATEGORIES; i++) {
        stats->output_filenames[i][0] = 0;
        stats->output_fp[i][0] = 0;
        stats->output_fp[i][1] = 0;
        stats->count_by_category[i] = 0;
        stats->count_by_category_long_enough[i] = 0;
        stats->count_by_category_too_short[i] = 0;
        stats->count_by_category_relaxed_hit[i] = 0;
        stats->count_by_category_external_clipped[i] = 0;
        stats->bases_written[i] = 0;
        stats->bases_before_clipping[i] = 0;

        for (j=0; j<MAX_READ_LENGTH; j++) {
            stats->read_length_counts[i][0][j] = 0;
            stats->read_length_counts[i][1][j] = 0;
            stats->read_pair_length_counts[i][j] = 0;
        }
    }

    for (i=0; i<=100; i++) {
        stats->gc_content[0][i] = 0;
        stats->gc_content[1][i] = 0;
    }

    stats->read_length = 0;
    stats->log_filename[0] = 0;
    stats->duplicates_log_filename[0] = 0;
    stats->output_prefix[0] = 0;
    stats->log_fp = 0;
    stats->duplicates_fp = 0;
    stats->num_read_pairs = 0;
    stats->n_duplicates = 0;
    stats->n_invalid_for_duplicate = 0;
    stats->duplicates_not_written = 0;
    stats->total_usable = 0;
    stats->gc_bases = 0;
    stats->at_bases = 0;
    stats->pairs_containing_n = 0;
}

/*----------------------------------------------------------------------*
 * Function:   initialise_junction_adaptor_alignment
 * Purpose:    Initialise structure
 * Parameters: result -> a JunctionAdaptorAlignment structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void initialise_junction_adaptor_alignment(JunctionAdaptorAlignment* result)
{
    result->score = -1;
    result->position = -1;
    result->total_matches = 0;
    result->total_alignment_length = 0;
    result->total_identity = 0;
    result->mismatches[0] = 0;
    result->matches[0] = 0;
    result->mismatches[1] = 0;
    result->matches[1] = 0;
    result->read_start= 0;
    result->read_end = 0;
    result->query_start = 0;
    result->query_end = 0;
    result->accepted = 0;
}

/*----------------------------------------------------------------------*
 * Function:   initialise_generic_adaptor_alignment
 * Purpose:    Initialise structure
 * Parameters: result -> a GenericAdaptorAlignment structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void initialise_generic_adaptor_alignment(GenericAdaptorAlignment* result)
{
    result->score = 0;
    result->position = 0;
    result->alignment_length = 0;
    result->identity = 0;
    result->mismatches = 0;
    result->matches = 0;
    result->read_start= 0;
    result->read_end = 0;
    result->query_start = 0;
    result->query_end = 0;
    result->accepted = 0;
}

/*----------------------------------------------------------------------*
 * Function:   reverse_compliment
 * Purpose:    Make reverse compliment of sequence.
 * Parameters: fwd -> sequence to produce reverse compliment of,
 *             rev -> resulting reverse compliment.
 * Returns:    Pointer to reverse compliment
 *----------------------------------------------------------------------*/
char* reverse_compliment(char* fwd, char* rev)
{
    int i;
    
    for (i=0; i<strlen(fwd); i++) {
        switch(fwd[strlen(fwd) - i - 1]) {
            case 'A': rev[i] = 'T'; break;
            case 'C': rev[i] = 'G'; break;
            case 'G': rev[i] = 'C'; break;
            case 'T': rev[i] = 'A'; break;
            default: printf("ERROR: Bad nucleotide in %s.\n", fwd); exit(2); break; 
        }
    }
    
    rev[strlen(fwd)] = 0;
    
    return rev;
}

/*----------------------------------------------------------------------*
 * Function:   size_hashtable
 * Purpose:    Calculate memory required for hashtable
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void size_hashtable(int* n_return, int* b_return)
{
    double utilisation = 0.8;
    double required_entries = approximate_reads / utilisation;
    int b = 100;
    double c = log(required_entries/b)/log(2);
    int n = ceil(c);
    int entries = pow(2.0, (double)n)*b;
    long int memory_bytes = entries*sizeof(Element);
    long int memory_meg = memory_bytes/(1024*1024);

#ifdef USE_MULTIPLE_HASHES
    memory_meg = memory_meg * NUMBER_OF_HASHES;
#endif
    
    memory_meg += 256;
    
    printf("                n: %d\n", n);
    printf("                b: %d\n", b);
    printf("          Entries: %d\n", entries);
    printf("       Entry size: %ld\n", sizeof(Element));
    
#ifdef USE_MULTIPLE_HASHES
    printf(" Number of tables: %d\n", NUMBER_OF_HASHES);
#endif

    printf("  Memory required: %ld MB\n\n", memory_meg);
    
    *n_return = n;
    *b_return = b;
}

/*----------------------------------------------------------------------*
 * Function:   usage
 * Purpose:    Report program usage.
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void usage(void)
{
    printf("Clip and analyse Illumina Nextera Long Mate Pair reads\n" \
           "\nSyntax: nextclip [-i r1.fastq] [-j r2.fastq] [-o prefix] [options]\n" \
           "\nOptions:\n" \
           "    [-d | --remove_duplicates] Remove PCR duplicates\n"
           "    [-e | --use_category_e] Use category E\n"
           "    [-h | --help] This help screen\n" \
           "    [-i | --input_one] Input FASTQ R1 file\n" \
           "    [-j | --input_two] Input FASTQ R2 file\n" \
           "    [-l | --log] Log filename\n" \
           "    [-m | --min_length] Minimum usable read length (default 25)\n" \
           "    [-n | --number_of_reads] Approximate number of reads (default 20,000,000)\n" \
           "    [-o | --output_prefix] Prefix for output files\n" \
           "    [-p | --only_duplicates] Only remove duplicates, don't trim\n" \
           "    [-q | --duplicates_log] PCR duplicates log filename\n" \
           "    [-r | --memory_requirements] Output memory requirements for specified number of reads\n" \
           "    [-t | --trim_ends] Trim ends of non-matching reads by amount (default 19)\n" \
           "    [-x | --strict_match] Strict alignment matches (default '34,18')\n" \
           "    [-y | --relaxed_match] Relaxed alignment matches (default '32,17')\n" \
           "\nComments/suggestions to richard.leggett@earlham.ac.uk\n" \
           "\n");
}

/*----------------------------------------------------------------------*
 * Function:   chomp
 * Purpose:    Remove hidden characters from end of line
 * Parameters: str -> String to change
 * Returns:    None
 *----------------------------------------------------------------------*/
void chomp(char* str)
{
    int i = strlen(str) - 1;
    
    while ((i > 0) && (str[i] < ' ')) {
        str[i--] = 0;
    }    
}

/*----------------------------------------------------------------------*
 * Function:   parse_csv_params
 * Purpose:    Parse comma separated arguments - eg. X,Y
 * Parameters: s -> string to parse
 *             a -> argument one
 *             b -> argument two
 * Returns:    None
 *----------------------------------------------------------------------*/
void parse_csv_params(char* s, int* a, int* b)
{
    char* comma = strchr(s, ',');
    *comma = 0;
    *a = atoi(s);
    *b = atoi(comma+1);    
}

/*----------------------------------------------------------------------*
 * Function:   parse_command_line
 * Purpose:    Parse command line options
 * Parameters: argc = number of arguments
 *             argv -> array of arguments
 * Returns:    None
 *----------------------------------------------------------------------*/
void parse_command_line(int argc, char* argv[], MPStats* stats)
{
    static struct option long_options[] = {
        {"remove_duplicates", no_argument, NULL, 'd'},
        {"use_category_e", no_argument, NULL, 'e'},
        {"help", no_argument, NULL, 'h'},
        {"input_one", required_argument, NULL, 'i'},
        {"input_two", required_argument, NULL, 'j'},
        {"log", required_argument, NULL, 'l'},
        {"min_length", required_argument, NULL, 'm'},
        {"number_of_reads", required_argument, NULL, 'n'},
        {"output_prefix", required_argument, NULL, 'o'},
        {"only_duplicates", no_argument, NULL, 'p'},
        {"duplicates_log", required_argument, NULL, 'q'},
        {"memory_requirements", no_argument, NULL, 'r'},
        {"adaptor_sequence", required_argument, NULL, 's'},
        {"trim_ends", required_argument, NULL, 't'},
        {"strict_match", required_argument, NULL, 'x'},
        {"relaxed_match", required_argument, NULL, 'y'},
        {0, 0, 0, 0}
    };
    int opt;
    int longopt_index;
    int i;
    
    if (argc == 1) {
        usage();
        exit(0);
    }
    
    while ((opt = getopt_long(argc, argv, "dehi:j:l:m:n:o:pq:rs:t:x:y:z:", long_options, &longopt_index)) > 0)
    {
        switch(opt) {
            case 'd':
                remove_duplicates=1;
                break;
            case 'e':
                use_category_e = 1;
                break;
            case 'h':
                usage();
                exit(0);
                break;
            case 'i':
                if (optarg==NULL) {
                    printf("Error: [-i | --input_one] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->input_filenames[0], optarg);
                break;
            case 'j':
                if (optarg==NULL) {
                    printf("Error: [-j | --input_two] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->input_filenames[1], optarg);
                break;
            case 'l':
                if (optarg==NULL) {
                    printf("Error: [-l | --log] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->log_filename, optarg);
                break;                
            case 'm':
                if (optarg==NULL) {
                    printf("Error: [-m | --minium_length] option requires an argument.\n");
                    exit(1);
                }
                minimum_read_size = atoi(optarg);
                break;
            case 'n':
                if (optarg==NULL) {
                    printf("Error: [-n | --number_of_reads] option requires an argument.\n");
                    exit(1);
                }
                approximate_reads = atoi(optarg);
                break;
            case 'o':
               if (optarg==NULL) {
                    printf("Error: [-o | --output_prefix] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->output_prefix, optarg);
                break;
            case 'p':
                duplicate_only_mode = true;
                remove_duplicates=1;
                break;
            case 'q':
                if (optarg==NULL) {
                    printf("Error: [-q | --duplicates_log] option requires an argument.\n");
                    exit(1);
                }
                strcpy(stats->duplicates_log_filename, optarg);
                break;
            case 'r':
                output_memory_requirements = true;
                break;
            case 's':
                if (optarg==NULL) {
                    printf("Error: [-s | --adaptor_sequence] option requires an argument.\n");
                    exit(1);
                }
                strcpy(single_junction_adaptor, optarg);
                break;
            case 't':
                if (optarg==NULL) {
                    printf("Error: [-t | --trim_ends] option requires an argument.\n");
                    exit(1);
                }
                trim_ends = atoi(optarg);
                break;                
            case 'x':
                if (optarg==NULL) {
                    printf("Error: [-x | --strict_match] option requires an argument.\n");
                    exit(1);
                }
                if (strchr(optarg, ',') == 0) {
                    printf("Error: [-x | --strict_match] option is of the format 'a,b'.\n");
                    exit(1);                    
                }
                parse_csv_params(optarg, &strict_double_match, &strict_single_match);
                break;
            case 'y':
                if (optarg==NULL) {
                    printf("Error: [-x | --relaxed_match] option requires an argument.\n");
                    exit(1);
                }
                if (strchr(optarg, ',') == 0) {
                    printf("Error: [-x | --relaxed_match] option is of the format 'a,b'.\n");
                    exit(1);                    
                }
                parse_csv_params(optarg, &relaxed_double_match, &relaxed_single_match);
                break;
            default:
                printf("Error: Unknown option %c\n", opt);
                exit(1);
                break;
        }
    }
    
    if (output_memory_requirements == true) {
        int n, b;
        size_hashtable(&n, &b);
        exit(0);
    }

    if (use_category_e == 1) {
        num_categories = 5;
    } else {
        num_categories = 4;
    }

    if (stats->output_prefix[0] != 0) {
        for (i=0; i<num_categories; i++) {
            sprintf(stats->output_filenames[i], "%s_%c", stats->output_prefix, i+'A');
        }
    }
    
    for (i=0; i<2; i++) {
        if (stats->input_filenames[i][0] == 0) {
            printf("Error: you must specify two input filenames\n");
            exit(2);
        }
    }
    
    for (i=0; i<num_categories; i++) {
        if (stats->output_filenames[i][0] == 0) {
            printf("Error: you must specify four output filenames\n");
            exit(2);
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:   strict_check
 * Purpose:    Perform the strict decision making on an alignment result
 * Parameters: result -> alignment result
 * Returns:    None
 *----------------------------------------------------------------------*/
void strict_check(JunctionAdaptorAlignment* result)
{
    if (result->total_matches >= strict_double_match) {
        result->accepted = 1;
        return;
    } 
   
    if ((result->matches[0] >= strict_single_match) || (result->matches[1] >= strict_single_match)) {
        result->accepted = 1;
        return;
    }
    
    result->accepted = 0;
}


/*----------------------------------------------------------------------*
 * Function:   relaxed_check
 * Purpose:    Perform the relaxed decision making on an alignment result
 * Parameters: result -> alignment result
 * Returns:    None
 *----------------------------------------------------------------------*/
void relaxed_check(JunctionAdaptorAlignment* result)
{
    if (result->total_matches >= relaxed_double_match) {
        result->accepted = 1;
        return;
    } 
    
    if ((result->matches[0] >= relaxed_single_match) || (result->matches[1] >= relaxed_single_match)) {
        result->accepted = 1;
        return;
    }    
}

/*----------------------------------------------------------------------*
 * Function:   find_sequence_in_read
 * Purpose:    Find a sequence within a read, eg. external adaptor
 * Parameters: read -> read to find sequence in
 *             sequence -> sequence to look for
 *             result -> alignment result
 * Returns:    None
 *----------------------------------------------------------------------*/
void find_sequence_in_read(FastQRead* read, char* sequence, GenericAdaptorAlignment* result)
{
    int x, p;
    int seq_length = strlen(sequence);
    
    // Initialise a result structure to store information
    initialise_generic_adaptor_alignment(result);
    result->read_size = read->read_size;
    
    // Start searching for the sequence... x is the position in the read where we start to compare the sequence
    for (x=-seq_length+5; x<read->read_size-5; x++) {
        int score = 0;
        int read_start = -1;
        int read_end = -1;
        int query_start = -1;
        int query_end = -1;
        int matches = 0;
        int mismatches = 0;
        
        // Go through each base of sequence, count matches and store start and end of match
        for (p=0; p<seq_length; p++) {
            if (((x+p) >= 0) && ((x+p) < read->read_size)) {
                if (sequence[p] == read->read[x+p]) {
                    matches++;
                    if (read_start == -1) {
                        read_start = x+p;
                        query_start = p;
                    } else {
                        read_end = x+p;
                        query_end = p;
                    }
                } else {
                    mismatches++;
                }
            }
        }
        
        // Score is simply matches for now
        score = matches;
        
        // Is this the best score yet?
        if (score > result->score) {
            result->score = score;
            result->matches = matches;
            result->mismatches = mismatches;
            result->position = x;
            result->read_start = read_start;
            result->read_end = read_end;
            result->query_start = query_start;
            result->query_end = query_end;
            result->alignment_length = 1 + (query_end - query_start);
            result->identity = 100.0 * result->matches / result->alignment_length;
        }
    }
    
    strcpy(result->adaptor, sequence);
    if ((result->alignment_length > 20) && (result->identity > 90)) {
        result->accepted = 1;
    } else {
        result->accepted = 0;
    }
}


/*----------------------------------------------------------------------*
 * Function:   find_junction_adaptors
 * Purpose:    Find junction adaptors in a read
 * Parameters: read -> read to find adaptors in
 *             result -> alignment result
 * Returns:    None
 *----------------------------------------------------------------------*/
void find_junction_adaptors(FastQRead* read, JunctionAdaptorAlignment* result)
{
    int x, p;
    int adaptor_length = strlen(duplicate_junction_adaptor);

    // Initialise a result structure to store information
    initialise_junction_adaptor_alignment(result);
    result->read_size = read->read_size;
    
    // Start searching for the transposon... x is the position in the read where we start to compare the transposon
    for (x=-adaptor_length+5; x<read->read_size-5; x++) {
        int matches[2] = {0, 0};
        int mismatches[2] = {0, 0};
        int score = 0;
        int read_start = -1;
        int read_end = -1;
        int query_start = -1;
        int query_end = -1;
        int better_result = 0;
 
        // Go through each base of transposon, count matches and store start and end of match
        // For interest, we store the 19nt sequence and it's reverse as a part 1 and part 2!
        for (p=0; p<adaptor_length; p++) {
            if (((x+p) >= 0) && ((x+p) < read->read_size)) {
                if (duplicate_junction_adaptor[p] == read->read[x+p]) {
                    matches[p < 19 ? 0:1]++;
                    if (read_start == -1) {
                        read_start = x+p;
                        query_start = p;
                    } else {
                        read_end = x+p;
                        query_end = p;
                    }
                } else {
                    mismatches[p < 19 ? 0:1]++;
                }
            }
        }
        
        // Score is simply matches for part 1 and 2
        score = matches[0] + matches[1];
        
        // Decide if it's a better match
        
        // If we've found a double match...
        if (score >= strict_double_match) {
            // Then see if it's a better double match than we've already got...
            if (score > result->score) {
                better_result = 1;
            }
        // Failing a double match, have we found a single match?
        } else if ((matches[0] >= strict_single_match) || (matches[1] >= strict_single_match)) {
            // Only consider if we haven't already found a double match
            if (result->score < strict_double_match) {
                int best_current_match = matches[0] > matches[1] ? matches[0]:matches[1];
                int best_result_match = result->matches[0] > result->matches[1] ? result->matches[0]:result->matches[1];
                
                // It's a better match if there are more bases of a single adaptor
                // OR, the same number of bases of a single adaptor, but a higher score.
                
                if (best_current_match > best_result_match) {
                    better_result = 1;
                } else if ((best_current_match == best_result_match) &&
                           (score > result->score)) {
                    better_result = 1;
                }
            }
        } else {
            // No matches... well, at least store best match, if not already got a good double or single match
            if ((result->score < strict_double_match) &&
                (result->matches[0] < strict_single_match) &&
                (result->matches[1] < strict_single_match)) {
                    if (score > result->score) {
                        better_result = 1;
                    }
            }
        }
    
        // Is this the best score yet?
        if (better_result == 1) {
            result->score = score;
            result->matches[0] = matches[0];
            result->mismatches[0] = mismatches[0];
            result->matches[1] = matches[1];
            result->mismatches[1] = mismatches[1];
            result->position = x;
            result->read_start = read_start;
            result->read_end = read_end;
            result->query_start = query_start;
            result->query_end = query_end;
            result->alignment_length[0] = query_start < 19 ? 19 - query_start : 0;
            result->alignment_length[1] = query_end >= 19 ? query_end-18 : 0;
            result->total_matches = matches[0] + matches[1];
            result->total_alignment_length = 1 + (query_end - query_start);
        }
    }
    
    // If we found a match, do some calculation
    if (result->score > 0) {
        result->total_identity = 100.0 * result->total_matches / result->total_alignment_length;
        result->identity[0] = 100.0 * result->matches[0] / result->alignment_length[0];
        result->identity[1] = 100.0 * result->matches[1] / result->alignment_length[1];
        strict_check(result);
    }
}

/*----------------------------------------------------------------------*
 * Function:   log_output_alignment
 * Purpose:    Output an alignment to the log
 * Parameters: stats -> MPStats structure
 *             read -> read adaptors were aligned to
 *             result -> junction adaptor alignment result
 *             external_adaptor_result -> external adaptor alignment result
 * Returns:    None
 *----------------------------------------------------------------------*/
void log_output_alignment(MPStats* stats, FastQRead* read, JunctionAdaptorAlignment* result, GenericAdaptorAlignment* external_adaptor_result)
{
    char s1[MAX_READ_LENGTH];
    char s2[MAX_READ_LENGTH];
    int i;
    int duplicate_junction_adaptor_length = strlen(duplicate_junction_adaptor);
    int external_adaptor_length =  strlen(external_adaptor_result->adaptor);

    fprintf(stats->log_fp, "\n---------- Read ID: %s ----------\n", read->read_header);
    fprintf(stats->log_fp, "%s\n", read->read);

    if ((result->score < 0) && (external_adaptor_result->score < 17)) {
        fprintf(stats->log_fp, "No alignment\n");
    } else {
        
        for (i=0; i<read->read_size; i++) {
            s1[i] = ' ';
            s2[i] = ' ';
        }
        
        if (result->score > 0) {
            for (i=0; i<read->read_size; i++) {
                int p = i - result->position;
                
                if ((p >= 0) && (p < duplicate_junction_adaptor_length)) {
                    if (duplicate_junction_adaptor[p] == read->read[i]) {
                        s1[i] = '|';
                    }
                    s2[i] = duplicate_junction_adaptor[p];
                }
            }
        }
        
        if (external_adaptor_result->score >= 17) {
            for (i=0; i<read->read_size; i++) {
                int p = i - external_adaptor_result->position;
                
                if ((p >= 0) && (p < external_adaptor_length)) {
                    if (external_adaptor_result->adaptor[p] == read->read[i]) {
                        s1[i] = '^';
                    }
                    s2[i] = external_adaptor_result->adaptor[p];
                }
            }
        }
   
        s1[read->read_size] = 0;
        s2[read->read_size] = 0;
        
        fprintf(stats->log_fp, "%s\n", s1);
        fprintf(stats->log_fp, "%s\n", s2);
    }
    
    fprintf(stats->log_fp, "Junction adaptor: Read %d-%d Adaptor %d-%d Score %d Id %.1f ", result->read_start, result->read_end, result->query_start, result->query_end, result->score, result->total_identity);
    fprintf(stats->log_fp, "(Matches %d,%d Mismatches %d,%d Lengths %d,%d Identity %.1f,%.1f)\n",
           result->matches[0], result->matches[1],
           result->mismatches[0], result->mismatches[1],
           result->alignment_length[0], result->alignment_length[1],
           result->identity[0], result->identity[1]);
    fprintf(stats->log_fp, "                  JUNCTION ADAPTOR %s\n", result->accepted == 1 ? "GOOD ALIGNMENT":"BAD ALIGNMENT");
    
    fprintf(stats->log_fp, "External adaptor: Read %d-%d Adaptor %d-%d Score %d Id %.1f\n", external_adaptor_result->read_start, external_adaptor_result->read_end, external_adaptor_result->query_start, external_adaptor_result->query_end, external_adaptor_result->score, external_adaptor_result->identity);
    fprintf(stats->log_fp, "                  EXTERNAL ADAPTOR %s\n", external_adaptor_result->accepted == 1 ? "GOOD ALIGNMENT":"BAD ALIGNMENT");
}

/*----------------------------------------------------------------------*
 * Function:   get_read
 * Purpose:    Read from a file into a FastQRead structure
 * Parameters: fp -> file to read from
 *             read -> structure to read into
 * Returns:    None
 *----------------------------------------------------------------------*/
int get_read(FILE* fp, FastQRead* read)
{
    int got_read = 1;
    
    if (!fgets(read->read_header, 1024, fp)) {
        got_read = 0;
    }

    if (!fgets(read->read, MAX_READ_LENGTH, fp)) {
        got_read = 0;
    }
    
    if (!fgets(read->quality_header, 1024, fp)) {
        got_read = 0;
    }
    
    if (!fgets(read->qualities, MAX_READ_LENGTH, fp)) {
        got_read = 0;
    }
    
    if (strlen(read->read) < MINIMUM_INPUT_READ_SIZE) {
        printf("Warning: read shorter than minimum read size (%d) - ignoring\n", MINIMUM_INPUT_READ_SIZE);
        got_read = 0;
    }

    if (got_read == 1) {
        chomp(read->read_header);
        chomp(read->read);
        chomp(read->quality_header);
        chomp(read->qualities);
        read->read_size=strlen(read->read);
        read->trim_at_base = read->read_size;
        read->trimmed_for_external_adaptor = false;
        read->trimmed_for_junction_adaptor = false;
    }
    
    return got_read;
}

/*----------------------------------------------------------------------*
 * Function:   write_read
 * Purpose:    Write read structure to file
 * Parameters: read -> FastQRead to write
 *             fp -> file to write to
 * Returns:    None
 *----------------------------------------------------------------------*/
void write_read(FastQRead* read, FILE *fp)
{
    fprintf(fp, "%s\n", read->read_header);
    fprintf(fp, "%s\n", read->read);
    fprintf(fp, "%s\n", read->quality_header);
    fprintf(fp, "%s\n", read->qualities);
}

/*----------------------------------------------------------------------*
 * Function:   decide_category
 * Purpose:    Decide what category a read pair is
 * Parameters: read_one -> FastQRead structure for read 1
 *             result_one -> alignment result for read 1
 *             read_two -> FastQRead structure for read 2
 *             result_two -> alignment result for read 22
 *             stats -> MPStats structure
 * Returns:    0 for category A, 1 for category B etc.
 *----------------------------------------------------------------------*/
int decide_category(MPStats* stats, FastQRead* read_one, JunctionAdaptorAlignment* result_one, FastQRead* read_two, JunctionAdaptorAlignment* result_two)
{
    int category = -1;
    
    if ((result_one->accepted == 1) && (result_two->accepted == 1)) {
        // Category A
        category = 0;
    } else if ((result_one->accepted == 0) && (result_two->accepted == 1)) {
        // Category B
        category = 1;
        // If using category E, check if we can get a relaxed hit on R1
        if (use_category_e == 1) {
            relaxed_check(result_one);
            if (result_one->accepted == 1) {
                if (result_one->read_start < read_one->trim_at_base) {
                    read_one->trim_at_base = result_one->read_start;
                }
                stats->count_by_category_relaxed_hit[category]++;
                category=4;
            }
        }
    } else if ((result_one->accepted == 1) && (result_two->accepted == 0)) {
        // Category C
        category = 2;
        // If using category E, check if we can get a relaxed hit on R2
        if (use_category_e == 1) {
            relaxed_check(result_two);
            if (result_two->accepted == 1) {
                if (result_two->read_start < read_two->trim_at_base) {
                    read_two->trim_at_base = result_two->read_start;
                }
                stats->count_by_category_relaxed_hit[category]++;
                category=4;
            }
        }
    } else {
        // Category D
        category = 3;
    }

    if ((read_one->trimmed_for_external_adaptor) || (read_two->trimmed_for_external_adaptor)) {
        stats->count_by_category_external_clipped[category]++;
    }
    
    return category;
}

/*----------------------------------------------------------------------*
 * Function:   trim_and_write_pair
 * Purpose:    Trim read and write to output file
 * Parameters: read_one -> FastQRead structure for read 1
 *             read_two -> FastQRead structure for read 2
 *             stats -> MPStats structure
 *             category = category of pair (0=A, 1=B etc.)
 * Returns:    None
 *----------------------------------------------------------------------*/
void trim_and_write_pair(MPStats* stats, int category, FastQRead* read_one, FastQRead* read_two)
{
    int l_one;
    int l_two;
    int smallest;
    
    // Count bases
    stats->bases_before_clipping[category] += strlen(read_one->read);
    stats->bases_before_clipping[category] += strlen(read_two->read);
    
    // Trim reads
    if (duplicate_only_mode == false) {
        read_one->read[read_one->trim_at_base] = 0;
        read_one->qualities[read_one->trim_at_base] = 0;
        read_two->read[read_two->trim_at_base] = 0;
        read_two->qualities[read_two->trim_at_base] = 0;
    }
    
    // Keep count
    stats->count_by_category[category]++;
    
    if (stats->log_fp != 0) {
        fprintf(stats->log_fp, "\nCategory %c\n\n", 'A' + category);
    }
    
    l_one = strlen(read_one->read);
    l_two = strlen(read_two->read);
    smallest = l_one < l_two ? l_one:l_two;
    
    stats->read_length_counts[category][0][l_one]++;
    stats->read_length_counts[category][1][l_two]++;
    
    stats->read_pair_length_counts[category][smallest]++;
    
    if ((strlen(read_one->read) < (minimum_read_size)) ||
        (strlen(read_two->read) < (minimum_read_size))) {
        stats->count_by_category_too_short[category]++;
        if (stats->log_fp != 0) {
            fprintf(stats->log_fp, "Too short\n");
        }
    } else {
        stats->count_by_category_long_enough[category]++;
        // Write reads
        write_read(read_one, stats->output_fp[category][0]);
        write_read(read_two, stats->output_fp[category][1]);
        stats->bases_written[category] += strlen(read_one->read);
        stats->bases_written[category] += strlen(read_two->read);
    }

}

/*----------------------------------------------------------------------*
 * Function:   check_read_ids
 * Purpose:    Check read 1 and read 2 have same ID
 * Parameters: read_one -> FastQRead structure for read 1
 *             read_two -> FastQRead structure for read 2
 *             stats -> MPStats structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void check_read_ids(MPStats* stats, FastQRead* read_one, FastQRead* read_two)
{
    int p=0;
    
    while ((p < strlen(read_one->read_header)) &&
           (p < strlen(read_two->read_header))) {
        if (read_one->read_header[p] != read_two->read_header[p]) {
            printf("Error: headers don't match up: %s and %s\n", read_one->read_header, read_two->read_header);
            exit(2);
        }
        
        if (read_one->read_header[p] <= ' ') {
            break;
        }
        else if (read_one->read_header[p] == '/') {
            break;
        } else {
            p++;
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:   valid_bases
 * Purpose:    Check if string is only A, C, G, T
 * Parameters: string -> string to check
 * Returns:    true if only A, C, T, G
 *----------------------------------------------------------------------*/
boolean valid_bases(char *string)
{
    int i;
    
    for (i=0; i<strlen(string); i++) {
        if ((string[i] != 'A') && (string[i] != 'C') && (string[i] != 'G') && (string[i] != 'T')) {
            return false;
        }
    }
    
    return true;
}

/*----------------------------------------------------------------------*
 * Function:   check_valid_bases_and_gc_content
 * Purpose:    Check string is only A, C, T and G, plus calculate GC
 * Parameters: string -> string to check
 *             gc -> GC return value
 * Returns:    true if only A, C, T, G
 *----------------------------------------------------------------------*/
boolean check_valid_bases_and_gc_content(char *string, int* gc)
{
    int i; 

    *gc = 0;

    for (i=0; i<strlen(string); i++) {
        if ((string[i] == 'G') || (string[i] == 'C')) {
            *gc = (*gc) + 1;
        } else if  ((string[i] != 'A') && (string[i] != 'T')) {
            return false;
        }
    }

    return true;
}

#ifdef USE_MULTIPLE_HASHES
/*----------------------------------------------------------------------*
 * Function:   check_pcr_duplicates
 * Purpose:    Calculate kmer signature and check for duplicates
 * Parameters: read_one -> FastQRead structure for read 1
 *             read_two -> FastQRead structure for read 2
 *             stats -> MPStats structure
 * Returns:    true if duplicate, false otherwise
 *----------------------------------------------------------------------*/
boolean check_pcr_duplicates(FastQRead* read_one, FastQRead* read_two, MPStats* stats)
{
    char kmer_string[TOTAL_KMER_SIZE+1];
    BinaryKmer kmer;
    BinaryKmer tmp_kmer;
    Element* e;
    boolean is_duplicate = 0;
    boolean found = false;
    int i;

    if ((!valid_bases(read_one->read)) || (!valid_bases(read_two->read))) {
        stats->pairs_containing_n++;
        return false;
    }
    
    if (kmer_offsets[0][0] == -1) {
        kmer_offsets[0][0] = 0;
        kmer_offsets[0][1] = read_one->read_size / 3;
        kmer_offsets[1][0] = 0;
        kmer_offsets[1][1] = 2 * read_one->read_size / 3;
        kmer_offsets[2][0] = read_one->read_size / 3;
        kmer_offsets[2][1] = 2 * read_one->read_size / 3;
    }
    
    for (i=0; i<NUMBER_OF_HASHES; i++) {
        strncpy(kmer_string,                        read_one->read + kmer_offsets[i][0], SEPARATE_KMER_SIZE);
        strncpy(kmer_string+(1*SEPARATE_KMER_SIZE), read_one->read + kmer_offsets[i][1], SEPARATE_KMER_SIZE);
        strncpy(kmer_string+(2*SEPARATE_KMER_SIZE), read_two->read + kmer_offsets[i][0], SEPARATE_KMER_SIZE);
        strncpy(kmer_string+(3*SEPARATE_KMER_SIZE), read_two->read + kmer_offsets[i][1], SEPARATE_KMER_SIZE);
        kmer_string[TOTAL_KMER_SIZE]=0;

        seq_to_binary_kmer(kmer_string, TOTAL_KMER_SIZE, &kmer);
        Key key = element_get_key(&kmer, TOTAL_KMER_SIZE, &tmp_kmer);
        e = hash_table_find_or_insert(key, &found, kmer_hashes[i]);
        if (e == NULL) {
            printf("Error: Hash table not big enough! Try specifying a larger number of reads.")
            exit(101);
        }
        
        if (found) {
            is_duplicate = true;
            e->count++;
        } else {
            e->flags = 1;
            e->count = 1;
        }
    }
    
    if (is_duplicate) {
        stats->n_duplicates++;
        if (stats->duplicates_fp) {
            fprintf(stats->duplicates_fp, "Match: %s\n", kmer_string);
            fprintf(stats->duplicates_fp, "   R1: %s\n", read_one->read);
            fprintf(stats->duplicates_fp, "   R2: %s\n\n", read_two->read);
        }
    }
    
    return is_duplicate;
}
#else
/*----------------------------------------------------------------------*
 * Function:   check_pcr_duplicates
 * Purpose:    Calculate kmer signature and check for duplicates
 * Parameters: read_one -> FastQRead structure for read 1
 *             read_two -> FastQRead structure for read 2
 *             stats -> MPStats structure
 * Returns:    true if duplicate, false otherwise
 *----------------------------------------------------------------------*/
boolean check_pcr_duplicates(FastQRead* read_one, FastQRead* read_two, MPStats* stats)
{
    char kmer_string[TOTAL_KMER_SIZE+1];
    BinaryKmer kmer;
    BinaryKmer tmp_kmer;
    Element* e;
    boolean is_duplicate = false;
    boolean found = false;
    int gc_one = 0;
    int gc_two = 0;
    
    if ((!check_valid_bases_and_gc_content(read_one->read, &gc_one)) || (!check_valid_bases_and_gc_content(read_two->read, &gc_two))) {
        stats->pairs_containing_n++;
        stats->n_invalid_for_duplicate++;
        return false;
    }

    stats->gc_bases+=gc_one;
    stats->gc_bases+=gc_two;
    stats->at_bases+=(read_one->read_size - gc_one);
    stats->at_bases+=(read_two->read_size - gc_two);
    
    if (gc_one > 0) {
        gc_one=(gc_one*100)/read_one->read_size;
    }
    
    if (gc_two > 0) {
        gc_two=(gc_two*100)/read_two->read_size;
    }
    
    stats->gc_content[0][gc_one]++;
    stats->gc_content[1][gc_two]++;

    strncpy(kmer_string, read_one->read + FIRST_KMER_OFFSET, SEPARATE_KMER_SIZE);
    strncpy(kmer_string+(1*SEPARATE_KMER_SIZE), (read_one->read) + (read_one->read_size / 2), SEPARATE_KMER_SIZE);
    strncpy(kmer_string+(2*SEPARATE_KMER_SIZE), read_two->read, SEPARATE_KMER_SIZE);
    strncpy(kmer_string+(3*SEPARATE_KMER_SIZE), (read_two->read) + (read_two->read_size / 2), SEPARATE_KMER_SIZE);
    //strncpy(kmer_string, read_one->read, TOTAL_KMER_SIZE);
    kmer_string[TOTAL_KMER_SIZE]=0;
    
    seq_to_binary_kmer(kmer_string, TOTAL_KMER_SIZE, &kmer);
    Key key = element_get_key(&kmer, TOTAL_KMER_SIZE, &tmp_kmer);
    e = hash_table_find_or_insert(key, &found, duplicate_hash);
    if (e == NULL) {
        printf("Error: Hash table not big enough! Try specifying a larger number of reads.");
        exit(101);
    }
    
    if (found) {
        e->count++;
        is_duplicate = true;
        if (stats->duplicates_fp) {
            fprintf(stats->duplicates_fp, "Match: %s\n", kmer_string);
            fprintf(stats->duplicates_fp, "   R1: %s\n", read_one->read);
            fprintf(stats->duplicates_fp, "   R2: %s\n\n", read_two->read);
        }
    } else {
        e->flags = 1;
        e->count = 1;
    }

    return is_duplicate;
}
#endif

/*----------------------------------------------------------------------*
 * Function:   process_files
 * Purpose:    Main function to process FASTQ files
 * Parameters: stats -> MPStats structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void process_files(MPStats* stats)
{
    FastQRead reads[2];
    JunctionAdaptorAlignment junction_adaptor_alignments[2];
    GenericAdaptorAlignment external_adaptor_alignments[2];
    int i, j;
    int is_duplicate;
    
    if (stats->log_filename[0] != 0) {
        stats->log_fp = fopen(stats->log_filename, "w");
        if (!stats->log_fp) {
            printf("Error: Can't open %s\n", stats->log_filename);
            exit(102);
        }
    }
    
    if (stats->duplicates_log_filename[0] != 0) {
        stats->duplicates_fp = fopen(stats->duplicates_log_filename, "w");
        if (!stats->duplicates_fp) {
            printf("Error: Can't open %s\n", stats->duplicates_log_filename);
            exit(102);
        }
    }
    
    
    // Open input files
    for (i=0; i<2; i++) {
        printf("Opening input filename %s\n", stats->input_filenames[i]);
        stats->input_fp[i] = fopen(stats->input_filenames[i], "r");
        if (!stats->input_fp[i]) {
            printf("Error: can't open file %s\n", stats->input_filenames[i]);
            exit(2);
        }
    }

    // Open output files
    for (i=0; i<num_categories; i++) {
        for (j=0; j<2; j++) {
            char filename[MAX_PATH_LENGTH];
            sprintf(filename, "%s_R%d.fastq", stats->output_filenames[i], j+1);
            printf("Opening output file %s\n", filename);
            stats->output_fp[i][j] = fopen(filename, "w");
            if (!stats->output_fp[i][j]) {
                printf("Error: can't open file %s\n", filename);
                exit(2);
            }
        }
    }
    
    // Read each entry in FASTQ files
    while (!feof(stats->input_fp[0])) {
        // Go next pair of reads
        int n_reads = 0;
        int category = -1;
        
        for (i=0; i<2; i++) {
            if (get_read(stats->input_fp[i], &reads[i]) == 1) {
                n_reads++;
            }

            if (stats->read_length == 0) {
                stats->read_length = strlen(reads[i].read);
            }
        }
        
        // Process pair
        if (n_reads == 2) {
            // Check read IDs match up
            check_read_ids(stats, &reads[0], &reads[1]);

            // Count pairs
            stats->num_read_pairs++;
            
            // Handle PCR duplicates
            is_duplicate = check_pcr_duplicates(&reads[0], &reads[1], stats);
            
            if ((remove_duplicates == 0) ||
                ((remove_duplicates == 1) && (is_duplicate == 0))) {

                if (stats->log_fp != 0) {
                    fprintf(stats->log_fp, "==================== New read pair ====================\n");
                }
                
                if (duplicate_only_mode == true) {
                    category = 3; // D
                } else {
                    for (i=0; i<2; i++) {
                        // Find junction adaptor
                        find_junction_adaptors(&reads[i], &junction_adaptor_alignments[i]);

                        // Look for external adaptor
                        find_sequence_in_read(&reads[i], external_adaptors[i], &external_adaptor_alignments[i]);

                        // Display log information
                        if (stats->log_fp != 0) {
                            log_output_alignment(stats, &reads[i], &junction_adaptor_alignments[i], &external_adaptor_alignments[i]);
                        }
                                        
                        // If junction adaptor found...
                        if (junction_adaptor_alignments[i].accepted == 1) {
                            // Count
                            stats->count_adaptor_found[i]++;
                            
                            // Trim
                            reads[i].trim_at_base = junction_adaptor_alignments[i].read_start;
                            reads[i].trimmed_for_junction_adaptor = true;
                        } else {
                            if (trim_ends > 0) {
                                reads[i].trim_at_base = reads[i].read_size - trim_ends;
                            }
                        }

                        // If external adaptor found...?
                        if (external_adaptor_alignments[i].accepted == 1) {
                            if (external_adaptor_alignments[i].read_start < reads[i].trim_at_base) {
                                reads[i].trim_at_base = external_adaptor_alignments[i].read_start;
                                reads[i].trimmed_for_external_adaptor = true;
                                if (reads[i].trimmed_for_junction_adaptor) {
                                    if (stats->log_fp != 0) {
                                        fprintf(stats->log_fp, "                  EXTERNAL ADAPTOR BEFORE JUNCTION ADAPTOR\n");
                                    }
                                }
                            }
                            
                            if (junction_adaptor_alignments[i].accepted == 1) {
                                stats->count_adaptor_and_external_found[i]++;
                            } else {
                                stats->count_external_only_found[i]++;
                            }
                        
                        }
                        
                        if (junction_adaptor_alignments[i].accepted == 1) {
                            if (reads[i].trim_at_base < (minimum_read_size)) {
                                stats->count_too_short[i]++;
                            } else {
                                stats->count_long_enough[i]++;
                            }
                        } else {
                            stats->count_no_adaptor[i]++;
                        }
                        
                    }
                    
                    // Decide category (A, B, C, D, E)
                    category = decide_category(stats, &reads[0], &junction_adaptor_alignments[0], &reads[1], &junction_adaptor_alignments[1]);
                }
                
                // Trim and write reads
                trim_and_write_pair(stats, category, &reads[0], &reads[1]);
            } else {
                stats->duplicates_not_written++;
            }
        } else if (n_reads == 1) {
            printf("Warning: Only managed to get one read - pair ignored\n");
        }
    }
    
    // Close files
    for (i=0; i<2; i++) {        
        fclose(stats->input_fp[i]);
    }
    
    for (i=0; i<num_categories; i++) {
        for (j=0; j<2; j++) {
            fclose(stats->output_fp[i][j]);
        }
    }
    
    if (stats->log_fp != 0) {
        fprintf(stats->log_fp, "\nDONE\n");
        fclose(stats->log_fp);
    }
    
    if (stats->duplicates_fp != 0) {
        fclose(stats->duplicates_fp);
    }
}

/*----------------------------------------------------------------------*
 * Function:   process_adaptor
 * Purpose:    Make double junction adaptor from adaptor and reverse
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void process_adaptor(void)
{
    char reverse[1024];

    reverse_compliment(single_junction_adaptor, reverse);
    strcpy(duplicate_junction_adaptor, single_junction_adaptor);
    strcat(duplicate_junction_adaptor, reverse);

    printf("Adaptor: %s\n\n", duplicate_junction_adaptor);
}

/*----------------------------------------------------------------------*
 * Function:   calculate_stats
 * Purpose:    Calculate percentages for stats
 * Parameters: stats -> MPStats structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void calculate_stats(MPStats* stats)
{
    int i;
    
    if (stats->num_read_pairs < 1) {
        printf("Error: Number of read pairs < 1!\n");
        exit(2);
    }
    
    if (stats->pairs_containing_n > 0) {
        stats->percent_pairs_containing_n = (100.0 * (double)stats->pairs_containing_n) / (double)stats->num_read_pairs;
    } else {
        stats->percent_pairs_containing_n = 0;
    }
    
    for (i=0; i<2; i++) {
        if (stats->count_adaptor_found[i] > 0) {
            stats->percent_adaptor_found[i] = (100.0 * (double)stats->count_adaptor_found[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_adaptor_found[i] = 0;
        }
        
        if (stats->count_no_adaptor[i] > 0) {
            stats->percent_no_adaptor[i] = (100.0 * (double)stats->count_no_adaptor[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_no_adaptor[i] = 0;
        }
        
        if (stats->count_long_enough[i] > 0) {
            stats->percent_long_enough[i] = (100.0 * (double)stats->count_long_enough[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_long_enough[i] = 0;
        }
        
        if (stats->count_too_short[i] > 0) {            
            stats->percent_too_short[i] = (100.0 * (double)stats->count_too_short[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_too_short[i] = 0;
        }
        
        if (stats->count_adaptor_and_external_found[i] > 0) {
            stats->percent_adaptor_and_external_found[i] = (100.0 * (double)stats->count_adaptor_and_external_found[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_adaptor_and_external_found[i] = 0;
        }

        if (stats->count_external_only_found[i] > 0) {
            stats->percent_external_only_found[i] = (100.0 * (double)stats->count_external_only_found[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_external_only_found[i] = 0;
        }
    }
    
    stats->total_too_short = 0;
    stats->total_long_enough = 0;   
 
    for (i=0; i<num_categories; i++) {
        if (stats->count_by_category[i] > 0) {
            stats->percent_by_category[i] = (100.0 * (double)stats->count_by_category[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_by_category[i] = 0;
        }
        
        if (stats->count_by_category_long_enough[i] > 0) {
            stats->percent_by_category_long_enough[i] = (100.0 * (double)stats->count_by_category_long_enough[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_by_category_long_enough[i] = 0;
        }
        
        if (stats->count_by_category_too_short[i] > 0) {
            stats->percent_by_category_too_short[i] = (100.0 * (double)stats->count_by_category_too_short[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_by_category_too_short[i] = 0;
        }
        
        if (stats->count_by_category_relaxed_hit[i] > 0) {
            stats->percent_by_category_relaxed_hit[i] = (100.0 * (double)stats->count_by_category_relaxed_hit[i]) / (double)stats->num_read_pairs;            
        } else {
            stats->percent_by_category_relaxed_hit[i] = 0;
        }
        
        if (stats->count_by_category_external_clipped[i] > 0) {
            stats->percent_by_category_external_clipped[i] = (100.0 * (double)stats->count_by_category_external_clipped[i]) / (double)stats->num_read_pairs;
        } else {
            stats->percent_by_category_external_clipped[i] = 0;
        }
        
        stats->total_too_short += stats->count_by_category_too_short[i];
        stats->total_long_enough += stats->count_by_category_long_enough[i];
    }
    
    // Total usable is categories A, B, C, E which are long enough
    stats->total_usable = stats->count_by_category_long_enough[0] + stats->count_by_category_long_enough[1] + stats->count_by_category_long_enough[2] + stats->count_by_category_long_enough[4];
    if (stats->total_usable > 0) {
        stats->percent_usable = (100.0 * (double)stats->total_usable) / (double)stats->num_read_pairs;
    } else {
        stats->percent_usable = 0;
    }
    
    if (stats->total_too_short > 0) {
        stats->percent_total_too_short = (100.0 * (double)stats->total_too_short) / (double)stats->num_read_pairs;
    } else {
        stats->percent_total_too_short = 0;
    }

    if (stats->total_long_enough > 0) {
        stats->percent_total_long_enough = (100.0 * (double)stats->total_long_enough) / (double)stats->num_read_pairs;
    } else {
        stats->percent_total_long_enough = 0;
    }
    
    if (stats->duplicates_not_written > 0) {
        stats->percent_duplicates_not_written = (100.0 * (double)stats->duplicates_not_written) / (double)stats->num_read_pairs;
    } else {
        stats->percent_duplicates_not_written = 0;
    }
    
    // GC content
    printf("GC bases: %ld  AT bases: %ld\n", stats->gc_bases, stats->at_bases);
    if (stats->gc_bases == 0) {
        stats->percent_gc = 0;
    } else if (stats->at_bases == 0) {
        stats->percent_gc = 100;
    } else {
        stats->percent_gc = (100.0 * stats->gc_bases) / (stats->gc_bases + stats->at_bases);
    }
}

/*----------------------------------------------------------------------*
 * Function:   report_stats
 * Purpose:    Output report of stats
 * Parameters: stats -> MPStats structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void report_stats(MPStats* stats)
{
    int i;
    
    printf("\nSUMMARY\n\n");
    printf("     Strict match parameters: %d, %d\n", strict_double_match, strict_single_match);

    if (use_category_e == 1) {
        printf("    Relaxed match parameters: %d, %d\n", relaxed_double_match, relaxed_single_match);
    }
    
    printf("           Minimum read size: %d\n", minimum_read_size);
    printf("                   Trim ends: %d\n", trim_ends);
    printf("\n");
    printf("        Number of read pairs: %d\n", stats->num_read_pairs);
    printf("   Number of duplicate pairs: %d\t%.2f %%\n", stats->n_duplicates, stats->percent_duplicates);
    printf("Number of pairs containing N: %ld\t%.2f %%\n", stats->pairs_containing_n, stats->percent_pairs_containing_n);
    
    for (i=0; i<2; i++) {
        printf("\n");
        printf("   R%d Num reads with adaptor: %d\t%.2f %%\n", i+1, stats->count_adaptor_found[i], stats->percent_adaptor_found[i]);
        printf("   R%d Num with external also: %d\t%.2f %%\n", i+1, stats->count_adaptor_and_external_found[i], stats->percent_adaptor_and_external_found[i]);
        printf("       R%d long adaptor reads: %d\t%.2f %%\n", i+1, stats->count_long_enough[i], stats->percent_long_enough[i]);
        printf("          R%d reads too short: %d\t%.2f %%\n", i+1, stats->count_too_short[i], stats->percent_too_short[i]);
        printf("     R%d Num reads no adaptor: %d\t%.2f %%\n", i+1, stats->count_no_adaptor[i], stats->percent_no_adaptor[i]);
        printf("  R%d no adaptor but external: %d\t%.2f %%\n", i+1, stats->count_external_only_found[i], stats->percent_external_only_found[i]);
    }
        
    for (i=0; i<num_categories; i++) {
        printf("\n");
        printf("   Total pairs in category %c: %d\t%.2f %%\n", 'A'+i, stats->count_by_category[i], stats->percent_by_category[i]);
        printf("         %c pairs long enough: %d\t%.2f %%\n", 'A'+i, stats->count_by_category_long_enough[i], stats->percent_by_category_long_enough[i]);
        printf("           %c pairs too short: %d\t%.2f %%\n", 'A'+i, stats->count_by_category_too_short[i], stats->percent_by_category_too_short[i]);
        printf("%c external clip in 1 or both: %d\t%.2f %%\n", 'A'+i, stats->count_by_category_external_clipped[i], stats->percent_by_category_external_clipped[i]);
        printf("     %c bases before clipping: %ld\n", 'A'+i, stats->bases_before_clipping[i]);
        printf("       %c total bases written: %ld\n", 'A'+i, stats->bases_written[i]);
    }
 
    printf("\n");
    printf("          Total usable pairs: %d\t%.2f %%\n", stats->total_usable, stats->percent_usable);
    printf("             All long enough: %d\t%.2f %%\n", stats->total_long_enough, stats->percent_total_long_enough);
    printf("    All categories too short: %d\t%.2f %%\n", stats->total_too_short, stats->percent_total_too_short);
    printf("      Duplicates not written: %d\t%.2f %%\n", stats->duplicates_not_written, stats->percent_duplicates_not_written);

    if (use_category_e == 1) {
        printf("\n");
        printf("         Category B became E: %d\t%.2f %%\n", stats->count_by_category_relaxed_hit[1], stats->percent_by_category_relaxed_hit[1]);
        printf("         Category C became E: %d\t%.2f %%\n", stats->count_by_category_relaxed_hit[2], stats->percent_by_category_relaxed_hit[2]);
    }

    printf("          Overall GC content: %.2f %%\n", stats->percent_gc);
    
    printf("\n");
}

/*----------------------------------------------------------------------*
 * Function:   output_histograms
 * Purpose:    Output historgrams for read lengths, GC etc.
 * Parameters: stats -> MPStats structure
 * Returns:    None
 *----------------------------------------------------------------------*/
void output_histograms(MPStats* stats)
{
    char filename[MAX_PATH_LENGTH];
    FILE* fp;
    int c, r, i;
    
    for (c=0; c<num_categories; c++) {
        for (r=0; r<2; r++) {            
            sprintf(filename, "%s_%c_R%d_hist.txt", stats->output_prefix, 'A'+c, r+1);
            fp = fopen(filename, "w");
            if (fp) {
                for (i=1; i<=stats->read_length; i++) {
                    fprintf(fp, "%d\t%d\n", i, stats->read_length_counts[c][r][i]);
                }
                fclose(fp);
            } else {
                printf("Error: can't open histogram file %s\n", filename);
            }
        }
        
        sprintf(filename, "%s_%c_pair_hist.txt", stats->output_prefix, 'A'+c);
        fp = fopen(filename, "w");
        if (fp) {
            long int cumulative[stats->read_length+1];
            
            // Make cumulative totals
            for (i=stats->read_length; i>=1; i--) {
                if (i==stats->read_length) {
                    cumulative[i]=0;
                } else {
                    cumulative[i]=cumulative[i+1];
                }
                cumulative[i]+=stats->read_pair_length_counts[c][i];
            }
            
            for (i=1; i<=stats->read_length; i++) {
                fprintf(fp, "%d\t%d\t%ld\n", i, stats->read_pair_length_counts[c][i], cumulative[i]);
            }
            fclose(fp);
        } else {
            printf("Error: can't open histogram file %s\n", filename);
        }
    }

    for (r=0; r<2; r++) {
        sprintf(filename, "%s_R%d_gc.txt", stats->output_prefix, r+1);
        fp = fopen(filename, "w");
        if (fp) {
            for (i=0; i<=99; i+=2) {
                fprintf(fp, "%d\t%d\n", i, stats->gc_content[r][i]+stats->gc_content[r][i+1]);
            }
            fclose(fp);
        } else {
            printf("Error: can't open GC histogram file %s\n", filename);
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:   calculate_pcr_duplicate_stats
 * Purpose:    Calculate PCR duplication statistics
 * Parameters: stats -> MPStats structure
 * Returns:    None
 *----------------------------------------------------------------------*/
#ifdef USE_MULTIPLE_HASHES
void calculate_pcr_duplicate_stats(MPStats* stats)
{
    int n_duplicates;
    double pc_duplicates;
    int t;

    void store_duplicates(Element * node) {
        if (node->count > 1) {
            n_duplicates+=(node->count-1);
        } 
    }

    printf("\n");
    for (t=0; t<3; t++) {
        printf("Counting duplicates for table %d\n", t);
        hash_table_print_stats(kmer_hashes[t]);
        n_duplicates=0;
        hash_table_traverse(&store_duplicates, kmer_hashes[t]);
        pc_duplicates=(double)((100.0*n_duplicates)/(double)stats->num_read_pairs);
        printf("\nDuplicates: %d (%.2f %%)\n", n_duplicates, pc_duplicates);
    }
 
    stats->percent_duplicates=(double)((100.0*stats->n_duplicates)/(double)stats->num_read_pairs);
}
#else
void calculate_pcr_duplicate_stats(MPStats* stats)
{
    FILE* fp;
    char filename[MAX_PATH_LENGTH];
    int duplicate_counts[MAX_DUPLICATES];
    int largest_count = 0;
    int i;
    int extra_duplicates = 0;
    
    if (stats->num_read_pairs <= 1) {
        printf("Error: Number of read pairs is too low!\n");
        exit(1);
    }
    
    void store_duplicates(Element * node) {
        int count = node->count;

        if (count > largest_count) {
            largest_count = count;
        }
        
        if (count >= MAX_DUPLICATES) {
            count = MAX_DUPLICATES - 1;
            printf("Warning: count (%d) exceeds maximum - treated as %d\n", count, MAX_DUPLICATES-1);
            extra_duplicates += (count + 1 - MAX_DUPLICATES);
        }
        
        duplicate_counts[count]++;
        
        if (count > 1) {
            stats->n_duplicates+=(count-1);
        }
    }
    
    printf("\n");
    hash_table_print_stats(duplicate_hash);
    printf("\nCounting duplicates...\n");
    
    for (i=0; i<MAX_DUPLICATES; i++) {
        duplicate_counts[i] = 0;
    }
    
    stats->n_duplicates=0;
    
    hash_table_traverse(&store_duplicates, duplicate_hash);
    
    if (stats->n_duplicates > 0 ) {
        stats->percent_duplicates = (double)((100.0*stats->n_duplicates)/(double)stats->num_read_pairs);
    } else {
        stats->percent_duplicates = 0;
    }

    sprintf(filename, "%s_duplicates.txt", stats->output_prefix);

    if (stats->num_read_pairs == extra_duplicates) {
        printf("Error: number of read pairs is equal to extra duplicates!\n");
        exit(1);
    }
    
    fp = fopen(filename, "w");
    if (fp) {
        fprintf(fp, "n\tCount\tPercent\n");
        for (i=1; ((i<MAX_DUPLICATES) && (i<=largest_count)); i++) {
            int count = (i*duplicate_counts[i]);
            double percent = 0;
            
            if (count > 0) {
                percent = (100.0*count)/(double)(stats->num_read_pairs - extra_duplicates);
            }
            
            if (i == 1) {
                count += stats->n_invalid_for_duplicate;
            }
            
            fprintf(fp, "%d\t%d\t%.2f\n", i, count, percent);
        }
        fclose(fp);
    } else {
        printf("Error: can't open duplicates file %s\n", filename);
    }
    
    if (extra_duplicates > 0) {
        printf("Note: %d pairs represented by duplicates of more than %d\n", extra_duplicates, MAX_DUPLICATES);
    }
}
#endif 

/*----------------------------------------------------------------------*
 * Function:   create_hash_table
 * Purpose:    Create hash table for storing PCR duplicates
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void create_hash_table(void)
{
    int b = 0;
    int n = 0;

    size_hashtable(&n, &b);

    printf("Creating hash tables for duplicate storage...\n");

#ifdef USE_MULTIPLE_HASHES
    int i;

    for (i=0; i<NUMBER_OF_HASHES; i++) {
        kmer_hashes[i] = hash_table_new(n, b, 25, TOTAL_KMER_SIZE);
        hash_table_print_stats(kmer_hashes[i]);
        kmer_offsets[i][0] = -1;
        kmer_offsets[i][1] = -1;
    }
#else
    duplicate_hash = hash_table_new(n, b, 25, TOTAL_KMER_SIZE);
    
    if (duplicate_hash == NULL) {
        printf("Error: No memory for hash table\n");
        exit(101);
    }
    
    hash_table_print_stats(duplicate_hash);
#endif
}


/*----------------------------------------------------------------------*
 * Function:   main
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    MPStats stats;
    time_t start, end;
    double seconds;
    
    time(&start);
    
    printf("\nNextClip v%s\n\n", NEXTCLIP_VERSION);
    
    initialise_stats(&stats);
    parse_command_line(argc, argv, &stats);
    create_hash_table();

    process_adaptor();
    process_files(&stats);
    calculate_stats(&stats);
    calculate_pcr_duplicate_stats(&stats);
    report_stats(&stats);
    output_histograms(&stats);
    
    time(&end);
    seconds = difftime(end, start);

    printf("\nDone. Completed in %.0f seconds.\n", seconds);
    
    return 0;
}
