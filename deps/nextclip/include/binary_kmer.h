/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo 
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team: 
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */

#ifndef BINARY_KMER_H_
#define BINARY_KMER_H_
#include <global.h>
#include <stdio.h>
#include <nucleotide.h>
#include <seq.h>

/*** Basic constants to allow compile-time configuration of what kmer-sizes we support ****/
#define ALLOW_KMERS_UP_TO_31 1
#define ALLOW_KMERS_UP_TO_63 2
#define ALLOW_KMERS_UP_TO_95 3
#define ALLOW_KMERS_UP_TO_127 4
#define ALLOW_KMERS_UP_TO_160 5
#define ALLOW_KMERS_UP_TO_192 6
#define ALLOW_KMERS_UP_TO_223 7
#define ALLOW_KMERS_UP_TO_255 8

#define LAST_NUCLEOTIDE_MASK 3ULL // 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0011
#define SINGLE_NUCLEOTIDE_SHIFT 2
#define BINVERSION 5

#ifdef INCLUDE_QUALITY_SCORES
#define TYPE_UNKNOWN  0
#define TYPE_SANGER   1
#define TYPE_SOLEXA   2
#define TYPE_ILLUMINA 3
#endif

// change this if you want to support different kmer-ranges
#ifndef NUMBER_OF_BITFIELDS_IN_BINARY_KMER
#define NUMBER_OF_BITFIELDS_IN_BINARY_KMER  ALLOW_KMERS_UP_TO_31 //Moved to the makefile
#endif

//platform specific
typedef uint64_t bitfield_of_64bits;

typedef bitfield_of_64bits BinaryKmer[NUMBER_OF_BITFIELDS_IN_BINARY_KMER]; //think of this as the number of long longs we encode the kmer in

typedef struct {
	int nkmers;
	BinaryKmer * kmer;
#ifdef INCLUDE_QUALITY_SCORES
	QualityString * quality_strings;
#endif
} KmerSlidingWindow;

/**
* A set of KmerArrays. 
* Strunct that holds the binary representation of a sequence. 
* To iterate and extract the binary kmers on a multithreaded
* environment, each thread needs to have it's own sliding window
* (which makes sense, since you don't want more than one thread 
* working in each sequence). 
*
*/

typedef struct {
	short kmer_size;
	int nwindows;
	int max_nwindows;
    int max_kmers;
	int i;
	int j;
	BinaryKmer  current;
	KmerSlidingWindow * window;
} KmerSlidingWindowSet;

// basic BinaryKmer operations
void binary_kmer_initialise_to_zero(BinaryKmer* bkmer);
void binary_kmer_assignment_operator(BinaryKmer left, BinaryKmer right);
boolean binary_kmer_comparison_operator(const BinaryKmer const left,
		const BinaryKmer const right);
boolean binary_kmer_less_than(const BinaryKmer const left,
		const BinaryKmer const right, short kmer_size);
void binary_kmer_right_shift(BinaryKmer* kmer, int num_bits_to_shift);
void binary_kmer_left_shift(BinaryKmer* kmer, int num_bits_to_shift,
		short kmer_size);
void binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(
		BinaryKmer* bkmer, Nucleotide n, short kmer_size);
boolean binary_kmer_modify_base(BinaryKmer* bkmer, Nucleotide n,
		short kmer_size, short possition);
Nucleotide binary_kmer_get_base_at_possition(BinaryKmer * bkmer,
		short possition);

Nucleotide char_to_binary_nucleotide(char c);

char binary_nucleotide_to_char(Nucleotide n);

char * nucleotides_to_string(Nucleotide * nucleotides, int length,
		char * string);

Nucleotide reverse_binary_nucleotide(Nucleotide n);



char * binary_kmer_to_seq(BinaryKmer* kmer, short kmer_size, char * seq);

BinaryKmer* seq_to_binary_kmer(char * seq, short kmer_size,
		BinaryKmer* prealloced_kmer);

char * seq_reverse_complement(char * in, int length, char * out);

BinaryKmer* binary_kmer_reverse_complement(BinaryKmer* kmer, short kmer_size,
		BinaryKmer* prealloc_reverse_kmer);

void binary_kmer_to_stream(BinaryKmer * bkmer, short kmer_size, FILE * f);

Nucleotide binary_kmer_get_last_nucleotide(BinaryKmer* kmer);

Nucleotide binary_kmer_get_first_nucleotide(BinaryKmer* kmer, short kmer_size);

//BinaryKmer binary_kmer_add_nucleotide_shift(BinaryKmer kmer,Nucleotide nucleotide, short kmer_size);

char reverse_char_nucleotide(char c);

void binary_kmer_right_shift_and_insert_new_base_at_left_end(BinaryKmer * kmer, Nucleotide n, short kmer_size);


void nucleotide_iterator(void(*f)(Nucleotide));

/**
* Block of sliding windows functions
* 
*/


void binary_kmer_free_kmers(KmerSlidingWindow * *);


/**
 *Block of functions for the sliding windows set. Ideally, the functions bellow are the 
 * "public" API.  
 */
KmerSlidingWindowSet * binary_kmer_sliding_window_set_new( int max_windows, int max_kmers);

KmerSlidingWindowSet * binary_kmer_sliding_window_set_new_from_read_length(short kmer_size,int max_read_length) ;


//get overlapping kmers from sequence
int get_sliding_windows_from_sequence(char * seq,  char * qualities, int length, char quality_cut_off, short kmer_size, KmerSlidingWindowSet * windows, 
				      int max_windows, int max_kmers, boolean break_homopolymers, int homopolymer_cutoff);

int get_sliding_windows_from_sequence_with_quality(char * seq,  char * qualities, int length, char quality_cut_off, short kmer_size, KmerSlidingWindowSet * windows, 
                                                   int max_windows, int max_kmers);

int get_sliding_windows(Sequence * seq, char quality_cutoff, KmerSlidingWindowSet * windows);

//use preallocated sliding window, and get all the kmers from the passed-in sequence. Any kmer that would have contained an N is returned as NULL
int get_single_kmer_sliding_window_from_sequence(char * seq, int length,
		short kmer_size, KmerSlidingWindow* kmer_window);
		
void binary_kmer_alloc_kmers_set(KmerSlidingWindowSet * windows,
		int max_windows, int max_kmers);
		
void binary_kmer_free_kmers_set(KmerSlidingWindowSet * *);

boolean binary_kmer_sliding_window_set_get_next(KmerSlidingWindowSet * ksws);

void binary_kmers_sliding_window_reset_iterator(KmerSlidingWindowSet * ksws);

KmerSlidingWindowSet * new_kmer_sliding_window_set(short kmer_size,int max_read_length);

#endif /* BINARY_KMER_H_ */
