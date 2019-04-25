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
/*
 element.c -- implements the nodes of the dBruijn graph
 */

#include <stdlib.h>
#include <global.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <string.h>

#include <global.h>
#include <flags.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <element.h>

long long int visited_count = 0;

void element_assign(Element * e1, Element * e2)
{
	binary_kmer_assignment_operator((*e1).kmer, (*e2).kmer);
	e1->count = e2->count;
}

boolean element_is_key(Key key, Element e, short kmer_size)
{
	if (key == NULL) {
		printf
		    ("Do not call element_is_key wth a NULL pointer. Exiting\n");
		exit(1);
	}

    //char temp_left[1024];
    //char temp_right[1024];
    //printf("Comparing %s with %s\n", binary_kmer_to_seq(key, 31, temp_left), binary_kmer_to_seq(&(e.kmer), 31, temp_right));
    
	return binary_kmer_comparison_operator(*key, e.kmer);
}

#ifdef INCLUDE_QUALITY_SCORES
void element_preallocate_quality_strings(Element * e, int c, int n)
{
	QualityStringArray *qa = &e->quality_string_arrays[c];
    
    if (n > 0) {
        qa->quality_strings = calloc(n, sizeof(QualityString));
        if (qa->quality_strings == NULL) {
            printf("Can't allocate memory for quality strings.\n");
            exit(1);
        }
    }
    
	qa->number_of_strings = 0;
	qa->limit = n;
}

void element_add_quality_string(Element * e, short c, char *q)
{
	QualityStringArray *qa = &e->quality_string_arrays[c];

    // Sanity checking
    assert(e != NULL);
    assert(q != NULL);
    assert(c >= 0);
    assert(qa->number_of_strings >= 0);    
    
	if (qa->number_of_strings == qa->limit) {
		qa->limit += 1;
		if (qa->limit == 1) {
			qa->quality_strings = malloc(qa->limit * sizeof(QualityString));
		} else {
			qa->quality_strings = realloc(qa->quality_strings, qa->limit * sizeof(QualityString));
		}
		if (qa->quality_strings == NULL) {
			printf("Can't allocate memory for quality strings.\n");
			exit(1);
		}
	}
    
	qa->quality_strings[qa->number_of_strings].quality = malloc(strlen(q) + 1);
	if (qa->quality_strings[qa->number_of_strings].quality == NULL) {
		printf("Can't allocate memory for quality strings.\n");
		exit(1);
	}
	strcpy(qa->quality_strings[qa->number_of_strings].quality, q);
	qa->number_of_strings++;
}
#endif

/* No longer required?
boolean element_smaller(Element e1, Element e2) {
	//return e1.edges < e2.edges;
	return db_node_get_edges(e1) < db_node_get_edges(e2);
}
*/

//TODO - make API safer - this gets contents of hash table, not  copy
BinaryKmer *element_get_kmer(Element * e)
{
	return &(e->kmer);
}

Key element_get_key(BinaryKmer * kmer, short kmer_size, Key preallocated_key)
{

	BinaryKmer local_rev_kmer;
	binary_kmer_initialise_to_zero(&local_rev_kmer);

	binary_kmer_reverse_complement(kmer, kmer_size, &local_rev_kmer);

	if (binary_kmer_less_than(local_rev_kmer, *kmer, kmer_size)) {
		binary_kmer_assignment_operator(*((BinaryKmer *)
						  preallocated_key),
						local_rev_kmer);
	} else {
		binary_kmer_assignment_operator(*((BinaryKmer *)
						  preallocated_key), *kmer);
	}

	return preallocated_key;

}

void element_initialise(Element * e, BinaryKmer * kmer, short kmer_size)
{

	BinaryKmer tmp_kmer;

	binary_kmer_initialise_to_zero(&tmp_kmer);
	binary_kmer_assignment_operator(e->kmer, *(element_get_key(kmer, kmer_size, &tmp_kmer)));

	e->count=0;
}

Orientation db_node_get_orientation(BinaryKmer * k, Element * e, short kmer_size)
{
	if (binary_kmer_comparison_operator(e->kmer, *k) == true) {
		return forward;
	}

	BinaryKmer tmp_kmer;

	if (binary_kmer_comparison_operator(e->kmer,
					    *(binary_kmer_reverse_complement
					      (k, kmer_size,
					       &tmp_kmer))) == true) {
		return reverse;
	}

	printf
	    ("programming error - you have called  db_node_get_orientation with a kmer that is neither equal to the kmer in this node, nor its rev comp\n");
	char tmpseq1[kmer_size];
	char tmpseq2[kmer_size]; 
	printf("Arg 1 Kmer is %s and Arg 2 node kmer is %s\n",
	       binary_kmer_to_seq(k, kmer_size, tmpseq1),
	       binary_kmer_to_seq(&(e->kmer), kmer_size, tmpseq2));
	exit(1);

}

char binary_nucleotide_to_edge(Nucleotide base){
    return 1 << base;
}

Orientation opposite_orientation(Orientation o)
{
   // assert(o==forward || o == reverse); //TODO: temporary assert, to expensive to have it all the time. 
    
	return o ^ 1;

}

boolean element_check_for_flag_ALL_OFF(Element * node)
{
	return node->flags == ALL_OFF;
}

