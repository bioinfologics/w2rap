//
// Created by Bernardo Clavijo (TGAC) on 26/05/2017.
//

#ifndef SOAP_SCAFFOLDER_LINKANALISES_H
#define SOAP_SCAFFOLDER_LINKANALISES_H

#include "stdint.h"
#include "stdinc.h"
#include "newhash.h"
#include "kmerhash.h"
#include "extfunc.h"
#include "extvab.h"

typedef struct {
    uint32_t source;
    uint32_t source_pos;
    uint32_t dest;
    uint32_t dest_pos;
    unsigned char peGrad;
} PAIR_LINK;

static inline int compare_pl (const void * a, const void * b)
{
    if ( ((PAIR_LINK*)a)->source <  ((PAIR_LINK*)b)->source ) return -1;
    if ( ((PAIR_LINK*)a)->source >  ((PAIR_LINK*)b)->source ) return 1;
    if ( ((PAIR_LINK*)a)->source_pos <  ((PAIR_LINK*)b)->source_pos ) return -1;
    if ( ((PAIR_LINK*)a)->source_pos >  ((PAIR_LINK*)b)->source_pos ) return 1;
    if ( ((PAIR_LINK*)a)->dest <  ((PAIR_LINK*)b)->dest ) return -1;
    if ( ((PAIR_LINK*)a)->dest >  ((PAIR_LINK*)b)->dest ) return 1;
    if ( ((PAIR_LINK*)a)->dest_pos <  ((PAIR_LINK*)b)->dest_pos ) return -1;
    if ( ((PAIR_LINK*)a)->dest_pos >  ((PAIR_LINK*)b)->dest_pos ) return 1;
    return 0;

}

static PAIR_LINK * pair_links;
static size_t pair_links_capacity=100000;
static size_t pair_links_size=0;

size_t create_all_connections();

int connection_prob(size_t source, size_t dest, long * dist);
uint8_t find_best_distance(size_t source, size_t dest, long * dist);
inline int max_distance(size_t linkid);
inline int min_distance(size_t linkid);
inline int proposed_distance(size_t linkid);


#endif //SOAP_SCAFFOLDER_LINKANALISES_H

