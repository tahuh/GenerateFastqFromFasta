/*
 * gen.c
 * generate fastq read from FASTA file in really fast manner
 * Compile : gcc -O2 -g -o gen gen.c -lz
 *
 * Requires Heng Li's kseq.h
 * Author : Thomas Sunghoon Heo
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <getopt.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

char *fasta = NULL;
//char *fastq = NULL;
//FILE *ofp = NULL;
kseq_t *p = NULL;
unsigned int READ_LENGTH = 150;
unsigned int SLIDE_SIZE = 1;

const char *usage = "gen [options] in.fasta\n"
"\n"
"[options]\n"
"\n"
"   -h         This message\n"
"   -s   INT   Slide size. Less than -l\n"
"   -l   INT   READ LENGTH\n"
"\n";
int main(int argc, char **argv)
{
    if(argc < 2){
		fprintf(stderr, "%s", usage);
        exit(-1);
	}
    int c;
    while( (c = (getopt(argc,argv,"hs:l:"))) != -1){
		switch(c){
			case 'h': fprintf(stderr, "%s", usage); exit(-1);
			case 's': SLIDE_SIZE = atoi(optarg); break;
			case 'l': READ_LENGTH = atoi(optarg); break;
			//case 'o': fastq = optarg;break;
		}
	}
	//if(fastq != NULL){
	//	ofp = fopen(fastq, "w");
	//}
	fasta = argv[optind];
	gzFile fasta_file = gzopen(fasta, "r");
	p = kseq_init(fasta_file);
	int l = 0;
	unsigned int i = 0;
	char *read = (char *)malloc(sizeof(char) * (READ_LENGTH + 1));
	char *qual = (char *)malloc(sizeof(char) * (READ_LENGTH + 1));
	memset(qual, 'F', sizeof(char) * (READ_LENGTH + 1));
	unsigned long reads = 0;
	unsigned long total_length = 0;
	while((l=kseq_read(p)) >= 0){
		char *name = p->name.s;
		char *seq = p->seq.s;
		total_length += p->seq.l;
		fprintf(stderr, "[gen] Producing for chromosome=%s\n", name);
		for(i = 0; i < p->seq.l - READ_LENGTH - SLIDE_SIZE + 1; i += SLIDE_SIZE){
			memset(read,0, sizeof(char) * (READ_LENGTH+1));
			memcpy(read, p->seq.s + i, sizeof(char) * READ_LENGTH);
			fprintf(stdout, "@%s:%d-%d\n%s\n+\n%s\n", name,i,i+READ_LENGTH, read, qual);
			reads++;
		
			if (reads % 1000000 == 0){
				fprintf(stderr, "[gen] Generated %ld reads from %s\n", reads, argv[optind]);
			}
		}
	}
	fprintf(stderr, "[gen] Generated %ld reads from %s\n", reads, argv[optind]);
	float cov = ((float)reads * READ_LENGTH)/((float)total_length);
	fprintf(stderr, "[gen] Expected coverage = %.3lf X\n", cov);
	free(read);
	free(qual);
	kseq_destroy(p);
	gzclose(fasta_file);
	return 0;
}
