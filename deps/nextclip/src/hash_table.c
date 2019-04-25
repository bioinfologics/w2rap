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
 hash_table.c -- implementation
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <locale.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <global.h>
#include <flags.h>
#include <string.h>
#include <binary_kmer.h>
#include <element.h>
#include <hash_table.h>
#include <hash_value.h>
#include <logger.h>


HashTable * hash_table_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size){ 
	assert(kmer_size > 0);
    assert(kmer_size < NUMBER_OF_BITFIELDS_IN_BINARY_KMER * 32 );
	//HashTable *hash_table = malloc(sizeof(HashTable));
	HashTable *hash_table = calloc(1,sizeof(HashTable));
	if (hash_table == NULL) {
		fprintf(stderr,"ERROR: could not allocate hash table of size %qd\n", (long long) bucket_size * (1 << number_bits));
        exit(1);
		//return NULL;
	}
	
	hash_table->collisions = calloc(max_rehash_tries, sizeof(long long));
	if (hash_table->collisions == NULL) {
		fprintf(stderr,"ERROR: could not allocate memory\n");
		exit(1);
        //return NULL;
		//exit(1);
	}
	
	hash_table->unique_kmers = 0;
	hash_table->max_rehash_tries = max_rehash_tries;
	hash_table->number_buckets = (long long) 1 << number_bits;
	hash_table->bucket_size   = bucket_size;

	//calloc is vital - we want to make sure initialised to zero
	hash_table->table = calloc(hash_table->number_buckets * hash_table->bucket_size, sizeof(Element));
	
	if (hash_table->table == NULL) {
		fprintf(stderr,"ERROR: could not allocate hash table of size %qd\n",hash_table->number_buckets * hash_table->bucket_size);
		//return NULL;
		exit(1);
	}
	
	hash_table->next_element = calloc(hash_table->number_buckets, sizeof(int));
	if (hash_table->table == NULL) {
		fprintf(stderr,"ERROR: could not allocate array of pointers for next available element in buckets [%qd]\n",hash_table->number_buckets);
		//return NULL;
		exit(1);
	}
	hash_table->kmer_size      = kmer_size;
	hash_table->number_of_threads = 1;
#ifdef ENABLE_MARK_PAIR
    hash_table->supernode_link = NULL;
#endif
	return hash_table;
}

void hash_table_free(HashTable ** hash_table)
{ 
	free((*hash_table)->table);
	free((*hash_table)->next_element);
	free((*hash_table)->collisions);
	free(*hash_table);
	*hash_table = NULL;
}


// Lookup for key in bucket defined by the hash value. 
// If key is in bucket, returns true and the position of the key/element in current_pos.
// If key is not in bucket, and bucket is not full, returns the next available position in current_pos (and overflow is returned as false)
// If key is not in bucket, and bucket is full, returns overflow=true
boolean hash_table_find_in_bucket(Key key, long long * current_pos, boolean * overflow, HashTable * hash_table, int rehash){
	
	
	//add the rehash to the final bitfield in the BinaryKmer
	BinaryKmer bkmer_with_rehash_added;
	binary_kmer_initialise_to_zero(&bkmer_with_rehash_added);
	binary_kmer_assignment_operator(bkmer_with_rehash_added, *key);
	bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] =   bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]+ (bitfield_of_64bits) rehash;
	
	int hashval = hash_value(&bkmer_with_rehash_added,hash_table->number_buckets);
	
	
	boolean found = false;
	int i=0;                     //position in bucket
	*overflow    = false;
	*current_pos   = (long long) hashval * hash_table->bucket_size;   //position in hash table
	
    //printf("Start of bucket: %lld\n", *current_pos);
    //printf("Flag: %d\n", element_check_for_flag_ALL_OFF(&hash_table->table[*current_pos]));
    
	while( (i<hash_table->bucket_size) &&   // still within the bucket
		  (!element_check_for_flag_ALL_OFF(&hash_table->table[*current_pos]) )  && // not yet reached an empty space
		  (!found)
		  )
    {
        //printf("Loopy\n");
		
		//sanity check -- to avoid out of boundary access
		if (*current_pos >= hash_table->number_buckets * hash_table->bucket_size || *current_pos<0)
		{
			printf("out of bounds problem found in hash table_find_with_position\n");
			exit(1);
		}
		
		//element found
		
		
		if (element_is_key(key,hash_table->table[*current_pos], hash_table->kmer_size))
		{
			found = true;
		}
		else
		{
			(*current_pos)++;
			i++;
		}
		
    }
	
	
	if (i == hash_table->bucket_size)
    {
		*overflow = true;
    }
	
    //if (found == true) {
        //printf("Found element at %lld\n", *current_pos);
    //} else {
    //    printf("Overflowed\n");
    //}
    
	assert(!found || !(*overflow));
	return found;
}


//currently not used, and must add a test
boolean hash_table_apply_or_insert(Key key, void (*f)(Element *), HashTable * hash_table){
	if (hash_table == NULL) {
		puts("NULL table!");
		exit(1);
	}
	
	long long current_pos;
	Element element;
	boolean overflow;
	int rehash=0;
	boolean found;
	do
    {
		found = hash_table_find_in_bucket(key,&current_pos,&overflow, hash_table,rehash);
		
		if (!found)
		{
			if (!overflow)
			{
				//sanity check
				if (!element_check_for_flag_ALL_OFF(&hash_table->table[current_pos])){
					printf("Out of bounds - trying to insert new node beyond end of bucket\n");
					exit(1);
				}
				
				element_initialise(&element,key, hash_table->kmer_size);
				element_assign( &(hash_table->table[current_pos]),  &element); 
				hash_table->unique_kmers++;
			}
			else//overflow
			{
				rehash++;
				if (rehash>hash_table->max_rehash_tries)
				{
					fprintf(stderr,"too much rehashing!! Rehash=%d\n", rehash);
					exit(1);
				}
			}
		}
		else
		{
			f(&hash_table->table[current_pos]);
		}
		
		
    }while(overflow);
	
	return found;
	
}


void hash_table_traverse(void (*f)(Element *),HashTable * hash_table){
	long long i;
	
	printf("\n");
	
	long long one_percent =  (hash_table->number_buckets * hash_table->bucket_size) / 100 ;
	int percent = 0;
	
	log_progress_bar(0);
	for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
		
		if (!element_check_for_flag_ALL_OFF(&hash_table->table[i])){
			f(&hash_table->table[i]);
		}
		
		if(one_percent > 0){
			
			if(i % one_percent == 0){
				percent = ((double)i / (double)(hash_table->number_buckets * hash_table->bucket_size)) *100;
				log_progress_bar(percent);
			} 
		}
	}
	log_progress_bar(100);
	printf("\n");
	
}

void hash_table_traverse_with_args(void (*f)(Element *, void *),void ** args, HashTable * hash_table){
	long long i;
	
	printf("\n");
	
	long long one_percent =  (hash_table->number_buckets * hash_table->bucket_size) / 100 ;
	int percent = 0;
	
	log_progress_bar(0);
	for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
		
		if (!element_check_for_flag_ALL_OFF(&hash_table->table[i])){
			f(&hash_table->table[i], args[0]);
		}
		
		if(one_percent > 0){
			
			if(i % one_percent == 0){
				percent = ((double)i / (double)(hash_table->number_buckets * hash_table->bucket_size)) *100;
				log_progress_bar(percent);
			} 
		}
	}
	log_progress_bar(100);
	printf("\n");
	
}

void hash_table_traverse_bucket(long long bucket, void (*f)(Element *),HashTable * hash_table){
	long long i;
	
	assert(hash_table != NULL);
	assert(f != NULL);
	assert(bucket >= 0);
	
	long long first_element = hash_table->bucket_size * bucket;
	long long last_element = first_element + hash_table->bucket_size;
	
	//printf("Last %lld < %lld\n", last_element, hash_table->number_buckets * hash_table->bucket_size);
	assert(last_element-1 < hash_table->number_buckets * hash_table->bucket_size);
	
	
	for(i = first_element; i < last_element; i++){
		if (!element_check_for_flag_ALL_OFF(&hash_table->table[i])){
			f(&hash_table->table[i]);
        }
	}
}

void hash_table_traverse_bucket_with_args(long long bucket, void (*f)(Element *, void *), void * args, HashTable * hash_table){
	long long i;
	
	assert(hash_table != NULL);
	assert(f != NULL);
	assert(bucket >= 0);
	
	long long first_element = hash_table->bucket_size * bucket;
	long long last_element = first_element + hash_table->bucket_size;
	
	//printf("Last %lld < %lld\n", last_element, hash_table->number_buckets * hash_table->bucket_size);
	assert(last_element-1 < hash_table->number_buckets * hash_table->bucket_size);
	
	
	for(i = first_element; i < last_element; i++){
		if (!element_check_for_flag_ALL_OFF(&hash_table->table[i])){
			f(&hash_table->table[i], args);
        }
	}
}

/**
 * Dumps the hash table. It stores the size of the structs so we can 
 * validate when reading that the memory do corresponds to the structure
 * that we are writing. 
 * 
 */
void hash_table_dump_memory(char * filename, HashTable * hash){
	
	FILE * fp = fopen(filename, "wb");
	char * magic = MAGIC_TEXT;
	int size_ht = sizeof(HashTable);
	int size_e = sizeof(Element);
	int magic_size = strlen(magic);
	short version = HASH_VERSION;

	
	long long number_buckets = hash->number_buckets;
	int bucket_size = hash->bucket_size; 
	long long hash_size=number_buckets * bucket_size;	
	
	//Header stuff, this is enough information to prepare the hash table.  
	fwrite(magic, sizeof(char), magic_size, fp);
	fwrite(&version, sizeof(short), 1, fp);
	fwrite(&size_ht, sizeof(int), 1,fp);
	fwrite(&size_e, sizeof(int), 1,fp);
	fwrite(&hash->kmer_size, sizeof(short), 1, fp);
	fwrite(&number_buckets, sizeof(long long), 1, fp);
	fwrite(&bucket_size, sizeof(int), 1, fp);
	fwrite(&hash->max_rehash_tries, sizeof(int), 1, fp);
	fwrite(&hash->unique_kmers, sizeof(long long), 1, fp);
	
	
	//The actual data of the hash table. We are storing everything 
	fwrite(hash->table, size_e, hash_size, fp);
	fwrite(hash->next_element, sizeof(int), number_buckets, fp);
	fwrite(hash->collisions, sizeof(long long), hash->max_rehash_tries, fp);
	log_and_screen_printf("Hash dumped to : %s \n" , filename);
	hash_table_print_stats(hash);
	fclose(fp);
	
}

static void exit_while_reading(FILE * fp, char * filename){
		log_and_screen_printf ("Error while reading file: %s\n", filename); 
		fclose(fp);
		exit (-1);
}

static void validate_read(size_t readed, size_t expected, FILE * fp, char * filename){
	if(readed != expected){
			exit_while_reading(fp, filename);
	}
}

HashTable *  hash_table_read_dumped_memory(char * filename ){
	
	HashTable * hash = calloc(1, sizeof(HashTable));
	
	FILE * fp = fopen(filename, "rb");
	int magic_size = strlen(MAGIC_TEXT);
	char * magic = calloc(magic_size, sizeof(char));
	int size_ht;
	int size_e = sizeof(Element);
	size_t readed; 
	short version;
	long long number_buckets;
	int bucket_size; 
	long long hash_size;	
	
	
	if(fp == NULL){
		exit_while_reading(fp, filename);
	}

	//Header stuff, this is enough information to prepare the hash table.  
	readed = fread(magic, sizeof(char), magic_size, fp);
	validate_read(readed, magic_size,  fp,  filename);
	if(strcmp(magic, MAGIC_TEXT) != 0){
		log_and_screen_printf( "[hash_table_read_dumped_memory] Invalid magic number!\n");
		fclose(fp);
		exit(-1);
		
	}
	//#printf("%s\n", magic);
	
	readed = fread(&version, sizeof(short), 1, fp);
	validate_read(readed, 1,  fp,  filename);
	//#printf("%d\n", version);
	if(version != HASH_VERSION){
		log_and_screen_printf( "[hash_table_read_dumped_memory] Invalid version number!\n");
		exit_while_reading(fp, filename);
		
	}
	readed = fread(&size_ht, sizeof(int), 1,fp);
	validate_read(readed, 1,  fp,  filename);
	//#printf("%d\n", size_ht);
	if(size_ht != sizeof(HashTable)){
		log_and_screen_printf( "[hash_table_read_dumped_memory] Invalid size of hash table!\n");
		exit_while_reading(fp, filename);
		
	}
	
	readed = fread(&size_e, sizeof(int), 1,fp);
	validate_read(readed, 1,  fp,  filename);
	//#printf("%d\n", size_e);
	if(size_e != sizeof(Element)){
		log_and_screen_printf( "[hash_table_read_dumped_memory] Invalid size of element!\n");
		exit_while_reading(fp, filename);
		
	}
		
	readed = fread(&hash->kmer_size, sizeof(short), 1, fp);
	//printf("kmer size %d\n", hash->kmer_size);
	validate_read(readed, 1,  fp,  filename);
	
	readed = fread(&number_buckets, sizeof(long long), 1, fp);
	validate_read(readed, 1,  fp,  filename);
	//printf("number of buckets %lld\n",number_buckets);
	
	readed = fread(&bucket_size, sizeof(int), 1, fp);
	validate_read(readed, 1,  fp,  filename);
	//printf("bucket size%d \n",bucket_size);
	
	readed = fread(&hash->max_rehash_tries, sizeof(int), 1, fp);
	validate_read(readed, 1,  fp,  filename);
	//printf("hash->max_rehash_tries %d\n", hash->max_rehash_tries);
	
	readed = fread(&hash->unique_kmers, sizeof(long long), 1, fp);
	validate_read(readed, 1,  fp,  filename);
//	printf("hash->unique_kmers %lld\n", hash->unique_kmers);
	
	hash->number_buckets = number_buckets ;
	hash->bucket_size = bucket_size; 
	
	hash_size=number_buckets * bucket_size;	
	//printf("Hash size %lld \n", hash_size);
	
	//Allocating the table according to the description of the file
	hash->table = calloc(hash_size, sizeof(Element));
	hash->next_element = calloc(number_buckets, sizeof(int));
	hash->collisions = calloc(number_buckets, sizeof(long long));
	
	if(hash->table == NULL){
		log_and_screen_printf( "Unable to create hash table\n ");
		exit_while_reading(fp, filename);
	}
	
	//Reading the actual data. 
	readed = fread(hash->table, size_e, hash_size, fp);
	validate_read(readed, hash_size,  fp,  filename);
	
	readed = fread(hash->next_element, sizeof(int), number_buckets, fp);
	validate_read(readed, number_buckets,  fp,  filename);
	
	
	readed = fread(hash->collisions, sizeof(long long), hash->max_rehash_tries, fp);
	validate_read(readed, hash->max_rehash_tries,  fp,  filename);
	
	fclose(fp);
	
	log_and_screen_printf("Hash readed from: %s \n" , filename);
	hash_table_print_stats(hash);
	
	return hash;
}

void hash_table_n_buckets_traverse(int block, int number_of_blocks, void (*f)(Element *),HashTable * hash_table){
	long long i;
	
	printf("\n");
	
	long long buckets_to_iterate = hash_table->number_buckets / number_of_blocks;
	
	assert(buckets_to_iterate  * number_of_blocks == hash_table->number_buckets);
	assert(block < number_of_blocks);
	
	
	long long one_percent = buckets_to_iterate / 100 ;
	int percent = 0;
	
	log_progress_bar(0);
	
	long long first_bucket = buckets_to_iterate * block;
	long long last_bucket = first_bucket + buckets_to_iterate;
	
	for(i=first_bucket; i < last_bucket; i++){
		
		hash_table_traverse_bucket(i, f, hash_table);
		
		
		if(one_percent > 0){
			
			if(i % one_percent == 0){
				percent = ((double)i / (double)(buckets_to_iterate)) *100;
				printf("%d:", block);
				log_progress_bar(percent);
			} 
		}
	}
	
	log_progress_bar(100);
	printf("\n");
	
}

void hash_table_n_buckets_traverse_with_args(int block, int number_of_blocks, void (*f)(Element *, void *), void * args, HashTable * hash_table){
	long long i;
	
	printf("\n");
	
	long long buckets_to_iterate = hash_table->number_buckets / number_of_blocks;
	
	assert(buckets_to_iterate  * number_of_blocks == hash_table->number_buckets);
	assert(block < number_of_blocks);
	
	
	long long one_percent = buckets_to_iterate / 100 ;
	int percent = 0;
	
	log_progress_bar(0);
	
	long long first_bucket = buckets_to_iterate * block;
	long long last_bucket = first_bucket + buckets_to_iterate;
	
	for(i=first_bucket; i < last_bucket; i++){
		
		hash_table_traverse_bucket_with_args(i, f,args, hash_table);
		
		
		if(one_percent > 0){
			
			if(i % one_percent == 0){
				percent = ((double)i / (double)(buckets_to_iterate)) *100;
				printf("%d:", block);
				log_progress_bar(percent);
			} 
		}
	}
	
	log_progress_bar(100);
	printf("\n");
	
}
#ifdef THREADS
void hash_table_threaded_traverse_with_args( void (*f)(Element *, void *), void ** args, HashTable * hash_table){
	int i;
	int no_of_threads = hash_table->number_of_threads;
	printf("\n");
	
	
	pthread_t * threads = calloc(no_of_threads, sizeof(pthread_t));
	
	void *exec_thread( void *ptr ){
		int * b_ptr = (int *) ptr;
		int block;
		block =  *b_ptr;
        
		//    printf("Executing thread %d\n", block);
		hash_table_n_buckets_traverse_with_args(block,  no_of_threads, f, &args[block],hash_table);
		return b_ptr;
	}
	int * starts = calloc(no_of_threads, sizeof(int));
	int ret;
	
	for(i = 0; i < no_of_threads; i ++){
		starts[i] = i;
		ret = pthread_create( &threads[i], NULL, exec_thread, (void*) &starts[i]);
		assert(ret==0);
	}
	
	
	printf("\n");
	
	
	for(i = 0; i < no_of_threads; i ++){
		pthread_join(threads[i], NULL);
	}
	
	free(starts);
	free(threads);	
}

void hash_table_threaded_traverse( void (*f)(Element *),HashTable * hash_table){
	int i;
	int no_of_threads = hash_table->number_of_threads;
	printf("\n");
	
	
	pthread_t * threads = calloc(no_of_threads, sizeof(pthread_t));
	
	void *exec_thread( void *ptr ){
		int * b_ptr = (int *) ptr;
		int block;
		block =  *b_ptr;
		//    printf("Executing thread %d\n", block);
		hash_table_n_buckets_traverse(block,  no_of_threads, f, hash_table);
		return b_ptr;
	}
	int * starts = calloc(no_of_threads, sizeof(int));
	int ret;
	
	for(i = 0; i < no_of_threads; i ++){
		starts[i] = i;
		ret = pthread_create( &threads[i], NULL, exec_thread, (void*) &starts[i]);
		assert(ret==0);
	}
	
	
	printf("\n");
	
	
	for(i = 0; i < no_of_threads; i ++){
		pthread_join(threads[i], NULL);
	}
	
	free(starts);
	free(threads);	
}
 
#endif

long long hash_table_array_index_of_element(Element * element, HashTable * hash_table){
	//void * ptr_element = (void *) element;
	//void * ptr_hash = (void *) hash_table;
	long long ptr_distance =  (long long)element - (long long)hash_table->table;
	long long distance = (long long)ptr_distance / sizeof(Element);
	return distance;
}





Element * hash_table_find(Key key, HashTable * hash_table)
{
	if (hash_table == NULL) 
    {
		puts("hash_table_find has been called with a NULL table! Exiting");
		exit(1);
    }
	
	Element * ret = NULL;
	long long current_pos;
	boolean overflow;
	int rehash = 0;
	boolean found; 
	
	do
    {
		found = hash_table_find_in_bucket(key,&current_pos, &overflow, hash_table,rehash);
		
		if (found) //then we know overflow is false - this is checked in find_in_bucket
		{
			ret =  &hash_table->table[current_pos];
		}
		else if (overflow)
		{ //rehash
			rehash++; 
			if (rehash>hash_table->max_rehash_tries)
			{
				fprintf(stderr,"too much rehashing!! Rehash=%d\n", rehash);
				exit(1);
			}
		}
    } while(overflow);
	
	return ret;
}


Element * hash_table_find_or_insert(Key key, boolean * found,  HashTable * hash_table){
	
	if (hash_table == NULL) {
		puts("NULL table!");
		exit(1);
	}
	
	Element element;
	Element * ret = NULL;
	int rehash = 0;
	boolean overflow; 
	
	long long current_pos;
	
	do{
		
		*found = hash_table_find_in_bucket(key,&current_pos,&overflow,hash_table,rehash);
		
		if (! *found)
		{
			if (!overflow) //it is definitely nowhere in the hashtable, so free to insert
			{
				//sanity check
				if (!element_check_for_flag_ALL_OFF(&hash_table->table[current_pos]))
				{
					fprintf(stderr,"error trying to write on an occupied element\n");
					exit(1);
				}
				
				//insert element
				//printf("Inserting element at position %qd in bucket \n", current_pos);
				element_initialise(&element,key, hash_table->kmer_size);
				
				//hash_table->table[current_pos] = element; //structure assignment
				element_assign(&(hash_table->table[current_pos]) , &element);
				//printf("Assigned at: %lld\n", current_pos);
                
				ret = &hash_table->table[current_pos];
				hash_table->unique_kmers++;
				
			}
			else
			{ //overflow -> rehashing
				
				rehash++;
				if (rehash>hash_table->max_rehash_tries)
				{
					hash_table_print_stats(hash_table);
					fprintf(stderr,"too much rehashing!! Rehash=%d\n", rehash);
					exit(1);
				}
			}
		}
		else //it is found
		{
            //printf("WOOHOO! Found\n");
			ret = &hash_table->table[current_pos];
		}
	} while (overflow);
	
	hash_table->collisions[rehash]++;
    hash_table->calculated = false;
	return ret;
}


//this methods inserts an element in the next available bucket
//it doesn't check whether another element with the same key is present in the table
//used for fast loading when it is known that all the elements in the input have different key
Element * hash_table_insert(Key key, HashTable * hash_table){
	
	if (hash_table == NULL) {
		puts("NULL table!");
		exit(1);
	}
	
	Element element;
	Element * ret = NULL;
	int rehash = 0;
	boolean inserted = false;
	do{
		//add the rehash to the final bitfield in the BinaryKmer
		BinaryKmer bkmer_with_rehash_added;
		binary_kmer_initialise_to_zero(&bkmer_with_rehash_added);
		binary_kmer_assignment_operator(bkmer_with_rehash_added, *key);
		bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] =   bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]+ (bitfield_of_64bits) rehash;
		
		int hashval = hash_value(&bkmer_with_rehash_added,(int)hash_table->number_buckets);//RHRG: Not quite sure if we want to cast here...
		
		if (hash_table->next_element[hashval] < hash_table->bucket_size)
		{ //can insert element
			long long  current_pos   = (long long) hashval * hash_table->bucket_size + (long long) hash_table->next_element[hashval] ;   //position in hash table
			
			//sanity check
			if (!element_check_for_flag_ALL_OFF(&hash_table->table[current_pos])){
				printf("Out of bounds - trying to insert new node beyond end of bucket\n");
				exit(1);
			}
			
			
			element_initialise(&element,key, hash_table->kmer_size);
			element_assign( &(hash_table->table[current_pos]),  &element); 
			hash_table->unique_kmers++;
			hash_table->next_element[hashval]++;	
			ret = &hash_table->table[current_pos];
			inserted=true;
		}
		else
		{//rehash
			rehash++;
			if (rehash>hash_table->max_rehash_tries)
			{
				fprintf(stderr,"too much rehashing!! Rehash=%d\n", rehash);
				exit(1);
			}
		}
		
	} while (! inserted);
	hash_table->calculated = false;
	return ret;
}

float hash_table_percentage_occupied(HashTable * hash_table){
float cap = (float) hash_table->unique_kmers / (float) (hash_table->bucket_size	* (float) hash_table->number_buckets) * 100;
    return cap;
}

void hash_table_print_stats(HashTable * hash_table)
{
	log_and_screen_printf("Hash:\n unique kmers: %'lld\n", hash_table->unique_kmers);
	log_and_screen_printf(" Capacity: %'lld \n", (hash_table->bucket_size * hash_table->number_buckets));
	float cap = hash_table_percentage_occupied(hash_table);
	log_and_screen_printf(" Occupied: %3.2f%%\n", cap);
	float percentage_pruned = ((float)hash_table->pruned_kmers/(float)hash_table->unique_kmers)*100;
    log_and_screen_printf(" Pruned: %'lld (%3.2f%%)\n", hash_table->pruned_kmers, percentage_pruned);
    
    int k;
	log_and_screen_printf(" Collisions:\n");
	for(k=0;k<10;k++)
    {
		if (hash_table->collisions[k] != 0) {
			log_and_screen_printf("\t tries %'d: %'qd\n",k,hash_table->collisions[k]);
		}
    }
}



long long hash_table_get_unique_kmers(HashTable * hash_table)
{
	return hash_table->unique_kmers;
}

void hash_table_set_number_of_threads(int threads, HashTable * hash_table){
	hash_table->number_of_threads = threads;
}


void hash_table_set_number_of_reads(long long read_count, short colour,  HashTable * hash_table){
    hash_table->number_of_reads[colour] = read_count;
}
long long hash_table_get_number_of_reads(short colour, HashTable * hash_table){
    return hash_table->number_of_reads[colour];
}

void hash_table_add_number_of_reads(long long read_count, short colour, HashTable * hash_table){
    hash_table->number_of_reads[colour] += read_count;
}

