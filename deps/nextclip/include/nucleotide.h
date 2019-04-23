

#ifndef nucleotide_h
#define nucleotide_h


#ifndef SOLID


typedef enum {
	Adenine = 0, Cytosine = 1, Guanine = 2, Thymine = 3, Undefined = 4
    //Never used for the binary representation of a kmer!.
} Nucleotide;
#else
typedef enum {
	Cero = 0, One = 1, Two = 2, Three = 3, Undefined = 4
    //Never used for the binary represent ation of a kmer!.
} Nucleotide;

typedef enum {
	Adenine = 0, Cytosine = 1, Guanine = 2, Thymine = 3, Undef = 4
} NucleotideBaseSpace;


NucleotideBaseSpace binary_nucleotide_base_space_get_next_base(NucleotideBaseSpace nbs, Nucleotide n);

NucleotideBaseSpace char_to_binary_nucleotide_base_space(char c);

char binary_nucleotide_base_space_to_char(NucleotideBaseSpace n);

NucleotideBaseSpace binary_kmer_get_last_nucleotide_in_base_space(BinaryKmer * kmer, Orientation o, NucleotideBaseSpace n, short kmer_size);

NucleotideBaseSpace binary_kmer_to_base_seq(BinaryKmer * kmer, Orientation o, NucleotideBaseSpace n, char * seq, short kmer_size);

#endif

#endif
