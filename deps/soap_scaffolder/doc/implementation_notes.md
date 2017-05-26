# w2rap scaffolder implementation notes

```
WARNING: these notes are meant to reflect the current status of implementation but
can be out of sync at specific times. We still provide them in the repository for
everyone's convenience.
Most of these comments refer to the latest implementations, some of them will need
to be specifically activated using the experimental flag.
```
## Loading contigs from GFA

:checkered_flag: TODO for w2rap v1 :checkered_flag:, as of now we still use `s_prepare` from a fasta file.

Contigs are stored into `CONTIG * contig_array`. The structure is defined as:

```C
typedef struct contig
{
	unsigned int from_vt;           // first kmer (vertex?)
	unsigned int to_vt;             // last kmer (vertex?)
	unsigned int length;            // what it says on the tin
	unsigned short indexInScaf;     // the index in the scaffold
	unsigned char cvg;              // coverage
	unsigned char bal_edge: 2;      // 0->FW, 2->REV or 1->PALINDROME
	unsigned char mask : 1;         // ??
	unsigned char flag : 1;         // ??
	unsigned char multi: 1;         // ??
	unsigned char inSubGraph: 1;    // ??
	unsigned char bubbleInScaff: 1; // ??
	char * seq;                     // what it says on the tin (but uses tightString)
	CONNECT * downwardConnect;      // links to other contigs (linked list)
	preARC * arcs;					// ??
	STACK * closeReads;             // ??
	size_t first_link_out;
	char * seqstr;					// seq in ascii (uppercase), for simplicity, lenght+1 allocated, with \0.
} CONTIG;
```


## Read mapping

### Using `s_map`

:warning: This behaviour is activated with the experimental flag `-x`.

### Support for external mappers

:construction: TODO after w2rap v1. :construction:

## Scaffolding

:warning: This behaviour is activated with the experimental flag `-x`.

### Read mapping to link transformation

Reads are grouped as pairs, and every pair generates a 2 links between contigs (one for the FW direction, one for the REV direction).
On this stage the links are stored in a `PAIR_LINK` array.

```C
typedef struct {
    uint32_t source;
    uint32_t source_pos;
    uint32_t dest;
    uint32_t dest_pos;
    unsigned char peGrad;
} PAIR_LINK;
```

Function to check: `PE2LinksEXP` in `linkAnalyses.c`.

### Contig breakpoint detection

:construction: TODO after w2rap v1. :construction:

Sometimes contigs obviously contradict read linkage. They should be split to avoid producing missassemblies around them.

### Contig-to-contig link analysis

Pair links between two contigs within +/- 30% of their nominal insert size (i.e. the mode) are considered satisfied.

As of now, there are 2 conditions for a link between `Contig A` and `Contig B`:

- `Contig A` has a pair link to `Contig B`, with 40% of the neighbouring links in a window of +/- 250bp have the same destination.
- A distance in the range [-1000, 10000] can be found that satisfies 75% of all pair links between `Contig A` and `Contig B`.

Function to check: `create_all_connections` on `linkAnalyses.c`.

### Scaffold construction

Function to check: `Links2ScafEXP` in `orderConting.c`.

### Scaffold Improvement

#### Repetition expansion

:checkered_flag: TODO for w2rap v1 :checkered_flag:

#### Collaborative linking

## Output generation

:checkered_flag: TODO for w2rap v1 :checkered_flag:

Need to rewrite this using the full sequence so no N->C errors are introduced and negative gaps (overlaps) produce correct output.

The current GFA output is limited and should be rewriten.

`outputScafseq` is currently walking both ways from a random start point and marking contigs as used. This is unnecessary  and cumbersome. We should start from every contig without scaffold predecessors and just walk forward, re-walking an already-used contig should be considered an error if repetitions are expanded. If proper GFA output is produced, then the function to walk lines on the GFA can be implemented in complete independence.

Function to check: `outputScafSeq` in `prlReadFillGap.c`.