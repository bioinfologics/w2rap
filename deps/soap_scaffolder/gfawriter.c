//
// Created by Bernardo Clavijo (EI) on 10/05/2017.
//



#include "gfawriter.h"

void dump_gfa(const char * filename){
    FILE * outgfa=fopen(filename,"w");
    //dump contigs as sequences (only lenghts)
    for (uint64_t i=1; i<num_ctg; ++i){
        if (contig_array[i].bal_edge==2){
            fprintf(outgfa,"S\t%llu\t*\tLN:i:%u\n",i,contig_array[i].length);
            //int col=650000;
            //outputTightStr ( outgfa, contig_array[i].seq, 0, contig_array[i].length, contig_array[i].length, 0, &col  );
            //fprintf(outgfa,"\n");
        }
    }
    //Now dump links
    for (uint64_t i=1; i<num_ctg; ++i){
        for (CONNECT * conn=contig_array[i].downwardConnect; conn ; conn=conn->next){
            if (conn->deleted || conn->mask) continue;

            int c1=i,c2=conn->contigID;
            unsigned char c1d='+',c2d='+';
            if (contig_array[c1].bal_edge==0) {
                --c1;
                c1d='-';
            }
            if (contig_array[c2].bal_edge==0) {
                --c2;
                c2d='-';
            }
            fprintf(outgfa,"L\t%d\t%c\t%d\t%c\t0M\n",c1,c1d,c2,c2d);


        }
    }

    fclose(outgfa);

    //dump

}

