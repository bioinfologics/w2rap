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
 
#ifndef STDSEQ_H_
#define STDSEQ_H_

#include <stdio.h>
#include <global.h>
//#include <binary_kmer.h>



typedef struct
{
    char * name;
    int  start,end, length;
    char * seq;  // sequence 
    char * qual; // qualities 
    //TODO validate that the sequence never gets out of bounds...
    int max_length;
    int max_name_length;
#ifdef ALIGN
    int count;
#endif
    boolean upper_case;
    boolean check_quality_values;
    char qual_offset;
    int lowest_expected_value;
    int highest_expected_value;
    int highest_raw_read_value;
    void * header;
   
} Sequence;

typedef struct{
    void (* header_parser)(Sequence * seq); 
    char * (* get_index)(Sequence * seq); 
}header_function;


typedef struct{
	Sequence ** sequences;
	int capacity; 
	int total;
} SequenceArray;

typedef char dinucleotide;

typedef struct{
  long long dinculeotides[25][MAX_READ_LENGTH];	
  long long nucleotides[5][MAX_READ_LENGTH];
  long long qual_dinculeotides[25][MAX_READ_LENGTH];	
  long long qual_nucleotides[5][MAX_READ_LENGTH];
  long long total_bases;
  int  max_length;
  int reads;
}SequenceStats;

void sequence_set_quality_parameters(Sequence *seq, int quality_offset);

void clean_stats(SequenceStats * total);

void sequence_stats(SequenceStats * total, Sequence * sequence );

double base_content(Nucleotide n, SequenceStats * total);

boolean base_is_valid(char base, char code);

void  print_stats(FILE * f, SequenceStats * total );

//returns length of sequence read

//note that fastq format dosn't support "partial" reading of a long entry -- ask Zam or Mario
//that means we can only read full fastq entries
int read_sequence_from_fastq(FILE * fp, Sequence * seq, int max_read_length);

//this routine can read long sequences (eg full chromosomes) , this is implemented by reading the sequence in chunks
int read_sequence_from_fasta(FILE * fp, Sequence * seq, int max_chunk_length, boolean new_entry, boolean * full_entry, int offset);

//Reads a sequnece file and a qualities file
int read_sequence_from_fasta_and_qual(FILE * fp, FILE * fq, Sequence * seq, int max_read_length);

void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length, char offset);

void free_sequence(Sequence ** );

void shift_last_kmer_to_start_of_sequence(Sequence * sequence, int length, short kmer_size);

void sequence_clean(Sequence * seq);

void sequence_copy(Sequence * to, Sequence * from);

void sequence_set_name(char * name, Sequence * seq);

void sequence_append_name(char * name, Sequence * seq);

void sequence_append(char * s, char * q, Sequence * seq);

void sequence_to_upper_case(Sequence * seq);

void sequence_add_base(char s, char q, Sequence * seq);

void sequence_reverse(Sequence * seq);

void sequence_complement(Sequence * seq);

void sequence_reverse_complement(Sequence * seq);

char sequence_get_base(int pos, Sequence * seq);

char sequence_get_qual(int pos, Sequence * seq);

int sequence_get_length(Sequence * seq);

void sequence_mask(int from, int to, Sequence * seq);

void sequence_remove_low_quality(Sequence * seq, char threshold);

Sequence * sequence_new(int max_read_length, int max_name_length, char offset);

void sequence_print_fasta(FILE * f, Sequence * seq);

void sequence_print_fasta_subseq(FILE * f,int start, int end, Sequence * seq);

void sequence_print_fastq(FILE * f,  Sequence * seq);

void sequence_iterator(void(*f)(char, int), Sequence * seq);

int sequence_count_gaps(Sequence * seq, int max);

int sequence_count_homopolymer(boolean forward, int pos, Sequence * seq);

void sequence_remove_base(int pos, Sequence * seq);

void sequence_remove_base_up_to_limit(int pos,int limit ,Sequence * seq);

void sequence_remove_missing_last_bases(Sequence * seq);

void sequence_insert_base_up_to_limit(char base, int pos,int limit ,Sequence * seq);

void sequence_insert_base_up_to_limit(char base, int pos,int limit ,Sequence * seq);

void sequence_trim(int l, Sequence * seq);

int sequence_differences_with_mask(Sequence * seq, Sequence * ref);

int sequence_compare_with_ambiguity(Sequence * seq1, Sequence * seq2);

void sequence_merge_removing_ambiguity(Sequence * seq1, Sequence * seq2);

//Functions to deal with the array of sequences.... 
SequenceArray * sequence_array_new(int capacity);

void sequence_array_destroy(SequenceArray ** sa);

Sequence * sequence_array_get_sequence(int pos, SequenceArray * sa);

void sequence_array_add(int max_read_length, int max_name_length, SequenceArray * sa);

void sequence_array_clean(SequenceArray *sa);


//Sequences for the alignment algorithms
int sequence_prev_anchor_base(int pos, Sequence * seq);

int sequence_next_anchor_base(int pos, Sequence * seq);

int sequence_next_hompoplymer(int pos, int max_righ, Sequence *seq);

int sequence_prev_hompoplymer(int pos, int max_left, Sequence *seq);


#endif /* STDSEQ_H_ */
