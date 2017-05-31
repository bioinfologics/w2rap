//
// Created by Bernardo Clavijo (TGAC) on 26/05/2017.
//


#include "linkAnalyses.h"

/*************************************************
 Function:
    PE2LinksEXP
 Description:
    Updates connections between contigs based on alignment
    information of paired-end reads.
    EXPERIMENTAL: weights linkage through detection of unique destinations
    for pairs from a window within the read.
 Input:
    1. infile:      alignment information file of paired-end reads
 Output:
    None.
 Return:
    None.
 *************************************************/
void PE2LinksEXP ( char * infile )
{
    lineLen = 1024;
    char name[256], line[lineLen];

    const float pair_link_growth_rate=2;
    pair_links=calloc(pair_links_capacity,sizeof(PAIR_LINK));
    FILE * fp1;
    FILE * linkF;
    gzFile * fp2;
    int flag = 0;
    unsigned int j;
    printf ( "*****************************************************\nEXPERIMENTAL PE2Links starting.\n\n" );

    sprintf ( name, "%s.links", infile );
    linkF = ckopen ( name, "w" );

    if ( !pes ) { loadPEgrads ( infile ); }

    printf( "Loading all across-contigs pairs\n");

    gzFile * fp;
    sprintf ( name, "%s.readOnContig.gz", infile );
    fp = gzopen ( name, "r" );

    int current_grad=0;
    long long upper_bound=pes[0].PE_bound;
    long long readno,pre_readno=0;
    unsigned int contigno,pre_contigno;
    int pos,pre_pos;
    //unsigned char dir, pre_dir;
    long long sizedist[100001]={0};//up to 50K seems a good range of size (this goes + and - as in ABySS);
    long long wrong_direction=0;
    long long same=0;
    long long different=0;
    long long single_first=0;
    long long single_second=0;
    long long pre_bound=0;
    //TODO (in progress): load all reads, computing the size distribution as it goes
    char * r=gzgets ( fp, line, lineLen ); //discard first line
    while ( r!=NULL ){
        r = gzgets ( fp, line, lineLen );
        sscanf ( line, "%lld %d %d", &readno, &contigno, &pos );
        if (r==NULL || readno>upper_bound){
            //TODO print stats and compute distribution (save distribution to file)
            sprintf ( name, "%s.peGrad%d.hist", infile,current_grad );
            FILE * fhist=fopen(name,"w");
            size_t maxc=0;
            for (long i=1;i<100000;i++){
                if (sizedist[i]>10) fprintf(fhist,"%ld, %lld\n",i-50000,sizedist[i]);
                if (i>50000 && sizedist[i]>maxc){
                    pes[current_grad].insertS=i-50000;
                    maxc=sizedist[i];
                }
            }
            fclose(fhist);

            printf("PE Grad %d: %10lld pairs, insert size %d bp \n",current_grad, (upper_bound-pre_bound)/2, pes[current_grad].insertS);
            printf("           %10lld (%6.2f%% ) correct orientation within contig\n",same-wrong_direction,(same-wrong_direction)*200.0/(upper_bound-pre_bound));
            printf("           %10lld (%6.2f%% ) incorrect orientation within contig\n",wrong_direction,wrong_direction*200.0/(upper_bound-pre_bound));
            printf("           %10lld (%6.2f%% ) across contigs\n",different,different*200.0/(upper_bound-pre_bound));
            printf("           %10lld (%6.2f%% ) single end mapped (%lld first, %lld second)\n",single_first+single_second,(single_first+single_second)*200.0/(upper_bound-pre_bound),single_first,single_second);
            printf("           %10lld (%6.2f%% ) fully unmapped\n\n",((upper_bound-pre_bound)/2-same-different-single_first-single_second),((upper_bound-pre_bound)/2-same-different-single_first-single_second)*200.0/(upper_bound-pre_bound));

            same=wrong_direction=different=0,single_first=0,single_second=0;
            pre_bound=upper_bound;
            for (long i=1;i<100000;i++) sizedist[i]=0;
            if (r==NULL) break;
        }
        while (readno>upper_bound){ //current_grad completely read, go to next
            ++current_grad;
            upper_bound=pes[current_grad].PE_bound;
        }
        if (contig_array[contigno].bal_edge==1) continue; //skips things on palindrome contigs

        if ( readno % 2 == 1 ){
            pre_readno=readno;
            pre_contigno=contigno;
            pre_pos=pos;
            ++single_first;
        }
        else if ( pre_readno == readno - 1 ){
            --single_first;
            if (contigno==pre_contigno) {
                same++;
                wrong_direction++;
            }
            else if (contigno==getTwinCtg(pre_contigno)){
                same++;
                long long s=contig_array[pre_contigno].length-pos-pre_pos+50000;
                s=(s>0? s:0);
                s=(s<100000? s: 100000);
                ++sizedist[s];
            }
            else {
                different++;
                //insertion with block-based growth
                //TODO: this should actually be 4 instances (fw and reverse, and both combinations of src and dest)
                if (pair_links_size+2>=pair_links_capacity){
                    pair_links_capacity*=pair_link_growth_rate;
                    pair_links=realloc(pair_links,pair_links_capacity*sizeof(PAIR_LINK));
                }
                pair_links[pair_links_size].source=pre_contigno;
                pair_links[pair_links_size].source_pos=pre_pos;
                pair_links[pair_links_size].dest=getTwinCtg(contigno);
                pair_links[pair_links_size].dest_pos=contig_array[contigno].length-pos;
                pair_links[pair_links_size].peGrad=current_grad;
                ++pair_links_size;
                pair_links[pair_links_size].source=contigno;
                pair_links[pair_links_size].source_pos=pos;
                pair_links[pair_links_size].dest=getTwinCtg(pre_contigno);
                pair_links[pair_links_size].dest_pos=contig_array[pre_contigno].length-pre_pos;
                pair_links[pair_links_size].peGrad=current_grad;
                ++pair_links_size;
            }
        } else ++single_second;

    }

    printf( "%lld different-contig links for analysis, %lld bytes of memory used on link array\n", pair_links_size, pair_links_capacity*
                                                                                                                    sizeof(PAIR_LINK));



    //TODO: print reads stats on each

    printf( "Sorting all across-contigs pairs\n");

    qsort(pair_links,pair_links_size, sizeof(PAIR_LINK),compare_pl);

    printf( "Deduplicating all\n");
    size_t total_links[gradsCounter], unique_links[gradsCounter];
    for (char i=0; i<gradsCounter; ++i) total_links[i]=unique_links[i]=0;
    size_t next_w=1;
    ++total_links[pair_links[0].peGrad];
    ++unique_links[pair_links[0].peGrad];

    for (size_t i=1; i<pair_links_size; ++i){
        ++total_links[pair_links[i].peGrad];
        if (compare_pl(&pair_links[i],&pair_links[i-1])!=0 ){
            ++unique_links[pair_links[i].peGrad];
            pair_links[next_w++]=pair_links[i];

        }
    }
    pair_links_size=next_w;

    for (char i=0;i<gradsCounter;++i){
        printf ("Grad %d, %ld / %ld (%.2f%%) unique links\n",i,unique_links[i],total_links[i],unique_links[i]*100.0/total_links[i]);
    }

    gzclose ( fp2 );
    fclose ( linkF );

    //TODO fill in first-connection of a contig struct
    for (size_t i=1;i<=num_ctg;++i) contig_array[i].first_link_out=0;
    contig_array[pair_links[0].source].first_link_out=0;
    size_t linked_ctgs=1;
    for (size_t i=1; i<pair_links_size;++i) {
        if (pair_links[i].source!=pair_links[i-1].source){
            ++linked_ctgs;
            contig_array[pair_links[i].source].first_link_out=i;
        }
    }
    printf ("Grad %lld / %lld contigs have links\n",linked_ctgs,num_ctg);

    return;
}

int find_perfect_overlap(size_t source, size_t dest, int min, int max){

    for (long ss=max;ss>=min;--ss){
        if (ss>contig_array[source].length || ss>contig_array[dest].length || getTwinCtg(source)==source || getTwinCtg(dest)==dest) continue;
        if (memcmp(contig_array[source].seqstr+contig_array[source].length-ss,contig_array[dest].seqstr,ss)==0){
            return ss;
        }
    }
    return 0;
}

int proposed_distance(size_t linkid){
    return pes[pair_links[linkid].peGrad].insertS
           - (contig_array[pair_links[linkid].source].length-pair_links[linkid].source_pos)
           - pair_links[linkid].dest_pos;
}

int min_distance(size_t linkid){
    return pes[pair_links[linkid].peGrad].insertS * 7 / 10
           - (contig_array[pair_links[linkid].source].length-pair_links[linkid].source_pos)
           - pair_links[linkid].dest_pos;
}

int max_distance(size_t linkid){
    return pes[pair_links[linkid].peGrad].insertS * 13 / 10
           - (contig_array[pair_links[linkid].source].length-pair_links[linkid].source_pos)
           - pair_links[linkid].dest_pos;
}

/*************************************************
Function:
    find_best_distance
Description:
    Finds the distance that best satisfies all links between two contigs
Input:
    1. Source contig #
    2. Dest contig #
Output:
    3. Distance, if <0 this is a validated overlap, non-validated overlaps will produce dist=1
Return:
    % of links with their distances satisfied with +/- 20%
*************************************************/
//TODO: have a version of this that goes upstream to add info if needed from previous contigs
uint8_t find_best_distance(size_t source, size_t dest, long * dist) {
    unsigned satisfied[100001]; //how many satisfied links at a distances between -50K + 50K
    unsigned proposed[10001]; //how many proposed links at a distances between -50K + 50K

    size_t bs = 0, bs_start = 0, bs_end = 0;

    for (int i = 0; i <= 100000; ++i) satisfied[i] = proposed[i] = 0;
    size_t lc = 0;

    //create an array of all satisfaction intervals (ie. dist -20% insert, dist+20%insert)
    for (size_t ci = contig_array[source].first_link_out; pair_links[ci].source == source; ++ci) {
        if (pair_links[ci].dest == dest) {
            ++lc;
            int mind = min_distance(ci) + 50000;
            int maxd = max_distance(ci) + 50000;
            if (mind < 0 || maxd > 100000) continue;
            for (int d = mind; d <= maxd; ++d) ++satisfied[d];
            ++proposed[(mind+maxd)/2];
            //printf("min=%lld max=%lld proposed=%lld\n",mind,maxd,(mind+maxd)/2);
        }
    }
    if (lc==0) return 0;
    //find the interval of maximum satisfaction
    for (int i = 0; i <= 100000; ++i) {
        if (satisfied[i] > bs) {
            bs = satisfied[i];
            bs_start = i;
        }
        if (satisfied[i] == bs) bs_end = i;
    }

    //compute the mean of all satisfied links
    long long total_d = 0;
    long long sl = 0;
    for (size_t i = bs_start; i <= bs_end; ++i) {
        total_d += proposed[i] * (i - 50000);
        sl += proposed[i];
    }
    //printf("bs=%lld bs_start=%lld bs_end=%lld totald=%lld sl=%lld\n",bs, bs_start,bs_end, total_d,sl);
    if (sl>0) {
        *dist = total_d / sl;
    } else *dist=((long long) (bs_start+bs_end))/2-50000;
    //printf("%lld->%lld Distances between %lld and %lld satisfy %lld / %lld links, best distance is %ld\n",
    //source,dest, bs_start-50000, bs_end-50000, bs, lc, *dist);
    return (uint8_t) (bs*100/lc);
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

size_t create_all_connections(){
    printf( "Creating all direct connections between %d contigs:\n", num_ctg);
    size_t connections_count=0;
    //TODO for each contig, evaluate connections with each proposed linked contig in a 1-to-1 basis
    uint8_t * dests = malloc(sizeof(uint8_t)*(num_ctg+1));
    for (size_t i=1; i<=num_ctg;++i) {
        //get all links destinations
        bzero(dests,sizeof(uint8_t)*num_ctg);
        for (size_t ci = contig_array[i].first_link_out; pair_links[ci].source == i; ++ci){
            if (dests[pair_links[ci].dest]!=2) {
                if (dests[pair_links[ci].dest] == 1) {
                    long d1,d2;
                    int p1=connection_prob(i,pair_links[ci].dest,&d1);
                    int p2=connection_prob(getTwinCtg(pair_links[ci].dest),getTwinCtg(i),&d2);
                    if ( ( p1>40 || p2>40 ) ){//&& d1-d2<2000 && d2-d1<2000 ){
                        int dm;
                        if (find_best_distance(i,pair_links[ci].dest,&dm)>75) {
                            if (dm < -1000 || dm > 10000) {
                                printf("Error on d1=%ld (p1=%d) d2=%ld (p2=%d) dm=%d for connection between %lld and %lld\n",
                                       d1, p1, d2, p2, dm, i, pair_links[ci].dest);
                            }
                            add1Connect(i, pair_links[ci].dest, dm, p1, 0);
                            add1Connect(getTwinCtg(pair_links[ci].dest), getTwinCtg(i), dm, p2, 0);
                            ++connections_count;
                        }
                    }
                    dests[pair_links[ci].dest] = 2;
                } else dests[pair_links[ci].dest] = 1;
            }
        }

        //sort all destinations

        //if a destination appears twice, evaluate possible connection.


    }
    free(dests);
    printf( "%ld connections created\n", connections_count);
    return connections_count;
}