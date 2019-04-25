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
 hash_table.h

 open hash table implementation - ie every bucket has a predifined size
 overloads results in rehashing
 all the routines as prefixed with hash_table
 */

#ifndef HASH_H_
#define HASH_H_

#define MAGIC_TEXT "BINARY_HASH"
#define HASH_VERSION 1
#ifdef ENABLE_READ_PAIR
struct read_pair_descriptor_array;
#endif
typedef struct {
    long long number_buckets;
    long long unique_kmers;
    long long pruned_kmers;
    int bucket_size;
    Element * table;
    long long * collisions;
    int * next_element; //keeps index of the next free element in bucket
    int max_rehash_tries;
    int number_of_threads;
    short kmer_size;
    short max_double_y_complexity;
    short max_coverage_for_branches;
    boolean calculated;
    int number_of_reads[NUMBER_OF_COLOURS];
} HashTable;

HashTable * hash_table_new(int number_bits, int bucket_size,
		int max_rehash_tries, short kmer_size);

void hash_table_free(HashTable * * hash_table);

//if the key is present applies f otherwise adds a new element for kmer
boolean hash_table_apply_or_insert(Key key, void(*f)(Element*), HashTable *);

//applies f to every element of the table
void hash_table_traverse(void(*f)(Element *), HashTable *);

//applies f to every element of the table, we can pass arguments to the option. The array is designed for
//the multithreaded version, only the first element in the array of arguments is used. The signature
//compatibility is designed to facilitate the adoption of multithreading. 
void hash_table_traverse_with_args( void (*f)(Element *, void * args), void ** args, HashTable * hash_table);

//if the element is not in table create an element with key and adds it
Element * hash_table_find_or_insert(Key key, boolean * found,
		HashTable * hash_table);
Element * hash_table_insert(Key key, HashTable * hash_table);

void hash_table_print_stats(HashTable * db_hash);

float hash_table_percentage_occupied(HashTable * hash_table);

long long hash_table_get_unique_kmers(HashTable *);

void hash_table_dump_memory(char * filename, HashTable * db_graph);

//return entry for kmer
Element * hash_table_find(Key key, HashTable * hash_table);

//returns the index of the element in the hash.
long long hash_table_array_index_of_element(Element *, HashTable *);

//Sets the possible number of threads. It validates that is a multiple of the size
//of the hash table.
void hash_table_set_number_of_threads(int threads, HashTable * hash_table);


#ifdef THREADS
//Iterator that splits the hash table on the number of threads given to the
//table. Default is 1. 
void hash_table_threaded_traverse_with_vars( void (*f)(Element *, void * args), void ** args, HashTable * hash_table);

//Iterator that splits the hash table on the number of threads given to the
//table. Default is 1. 
void hash_table_threaded_traverse( void (*f)(Element *), HashTable * hash_table);
#endif

HashTable *  hash_table_read_dumped_memory(char * filename );


void hash_table_set_number_of_reads(long long read_count, short colour, HashTable * hash);
void hash_table_add_number_of_reads(long long read_count, short colour,  HashTable * hash);
long long hash_table_get_number_of_reads( short colour, HashTable * hash);

#endif /* HASH_H_ */
