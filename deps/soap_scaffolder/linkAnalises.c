//
// Created by Bernardo Clavijo (TGAC) on 26/05/2017.
//


#include "linkAnalises.h"


int find_perfect_overlap(size_t source, size_t dest, int min, int max){

    for (long ss=max;ss>=min;--ss){
        if (ss>contig_array[source].length || ss>contig_array[dest].length || getTwinCtg(source)==source || getTwinCtg(dest)==dest) continue;
        if (memcmp(contig_array[source].seqstr+contig_array[source].length-ss,contig_array[dest].seqstr,ss)==0){
            return ss;
        }
    }
    return 0;
}

/*************************************************
Function:
    find_best_distance
Description:
    Analyses the probability of a connection between 2 contigs
Input:
    1. Source contig #
    2. Dest contig #
Output:
    3. Distance, if <0 and return >0, this is a valid overlap
Return:
    Score (probability) of a connection, form 0 (not connected) to 1000 (completely supported)
*************************************************/
//TODO: have a version of this that goes upstream to add info if needed from previous contigs
boolean find_best_distance(size_t source, size_t dest, long * dist){
    //create an array of all satisfaction intervals (ie. dist -20% insert, dist+20%insert)
    //sort array


    //count for every possible distance what the satisfaction % is.

    //take the distance that satisfies the most links, compute the mean of all satisfied links
}



/*************************************************
Function:
    connection_prob
Description:
    Analyses the probability of a connection between 2 contigs
Input:
    1. Source contig #
    2. Dest contig #
Output:
    3. Distance, if <0 and return >0, this is a valid overlap
Return:
    Score (probability) of a connection, form 0 (not connected) to 1000 (completely supported)
*************************************************/
//TODO: have a version of this that goes upstream to add info if needed from previous contigs
int connection_prob(size_t source, size_t dest, long * dist){
    size_t window_size=250;//TODO: what is a good number for this?
    long total_distance=0;
    long confirmed_links=0;
    int max_conf=0;
    for (int g=0;g<gradsCounter;++g) {
        //size_t link_dest[1000];//position of the links
        size_t link_pos[10000];
        size_t confirming_neighbours[10000];
        size_t contradicting_neighbours[10000];
        size_t link_destpos[10000];
        size_t lc = 0;
        for (size_t ci = contig_array[source].first_link_out; pair_links[ci].source == source && lc<10000; ++ci) {
            if (pair_links[ci].dest==dest && pair_links[ci].peGrad==g) {
                //add link
                link_pos[lc] = pair_links[ci].source_pos;
                link_destpos[lc] = pair_links[ci].dest_pos;
                confirming_neighbours[lc] = 0;
                contradicting_neighbours[lc] = 0;
                //go back up to window_size, counting contradicting and confirming neighbours
                for (long n = ci - 1; n >= 0 && pair_links[n].source == source &&
                                      pair_links[n].source_pos >= link_pos[lc] - window_size; --n) {
                    if (pair_links[n].peGrad == g) {
                        if (pair_links[n].dest == dest) ++confirming_neighbours[lc];
                        else ++contradicting_neighbours[lc];
                    }
                }
                //go forward up to window_size, counting contradicting and confirming neighbours
                for (long n = ci + 1; n >= pair_links_size && pair_links[n].source == source &&
                                      pair_links[n].source_pos <= link_pos[lc] + window_size; ++n) {
                    if (pair_links[n].peGrad == g) {
                        if (pair_links[n].dest == dest) ++confirming_neighbours[lc];
                        else ++contradicting_neighbours[lc];
                    }
                }

                int conf;
                //TODO: this is a stupid way to count, but as of now if window value >40%, go for it
                if (confirming_neighbours[lc] > 0) {
                    conf = confirming_neighbours[lc] * 100 /
                           (confirming_neighbours[lc] + contradicting_neighbours[lc]);

                    if (conf >= 25) {
                        ++confirmed_links;
                        total_distance += (long) pes[g].insertS - ((long) contig_array[source].length - (long) link_pos[lc]) - (long) link_destpos[lc];
                        ++lc;
                        if (conf > max_conf) max_conf = conf;
                    }
                }
            }
        }
        //TODO: as of now, 2 confirmed links call it a connection

        if (confirmed_links>1){
            *dist=(long)total_distance/(long)confirmed_links;
            if(confirmed_links<10000 && -1000<*dist && 10000>*dist) {
                if (*dist<500){
                    int ovl=find_perfect_overlap(source,dest,10,500);
                    if (ovl>0){
                        //printf ("Perfect overlap validated between %ld and %ld, %dbp\n",source,dest,ovl);
                        *dist=-ovl;
                    } else if (*dist<-500){
                        return 0;
                    } else *dist=1;
                }
                return max_conf;
            }
        }
        //find out first and last connection between source and dest for each library

        //compute proposed distance for library
        //compute expected connection mapping locations for library
        //do window-based linkage across the mapping location
    }
    //compose

    //validate overlap
    return 0;
}