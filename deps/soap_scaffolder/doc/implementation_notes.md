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

### Contig breakpoint detection

:construction: TODO after w2rap v1. :construction:

Sometimes contigs obviously contradict read linkage. They should be split to avoid producing missassemblies around them.

### Contig-to-contig link analysis



### Scaffold construction

### Scaffold Improvement

#### Repetition expansion

#### Collaborative linking

## Output generation