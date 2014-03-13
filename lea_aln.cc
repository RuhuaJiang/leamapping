/*
 * lea_aln.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: jiang
 */
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include "lea_types.h"
#include "lea_seqio.h"
#include "lea_utility.h"

int lea_map_core( ReferenceInfo reference_info , TableCell *kmer_position_table,read_t *readsChunk, uint32_t n_seqs, Options *op) {
	uint32_t i = 0;
	read_t *p;

	while ((uint32_t) i < (uint32_t) n_seqs) {
		p = &readsChunk[i++]; //A read.
	}

	free_read_seq(n_seqs, readsChunk);
	return 0;
}

int lea_map(char *refFile, char *readsFile, Options opt) {
	int ret, i;
	TableCell *kmer_position_table;
	ReferenceInfo reference_info;
	char *tableFn;
	uint32_t total = 0, n_seqs = 0;
	uint32_t readsBp = 0;
	seqio_t *seq;
	read_t *readsChunk, *p;
	char *prefix = refFile;
	clock_t t, mapStart;


	reference_info.bns = bns_restore(prefix);


	opt.chunck_size =  (0x2ffff);
	opt.index_parameter.kmer_len = get_kmer_length(reference_info.bns->l_pac);

	opt.kmer_table_len =  pow(2,(opt.index_parameter.kmer_len * 4))+10;
	opt.index_parameter.index_reverse_complement = false;
	opt.index_parameter.distinct_level = 1;
	opt.index_parameter.char_size = 2;

	fprintf(stderr, "[dmap_aln] start \n");
	t = clock();
	fprintf(stderr, "[dmap_aln] retrieving index info ...\n");

    tableFn = (char*)calloc(strlen(refFile) + 10, 1);
	strcpy(tableFn, refFile);
	strcat(tableFn, ".tb");
	reference_info.pac = static_cast<uint8_t *>(calloc((reference_info.bns->l_pac) / 4 + 1, 1));

	if (reference_info.pac == 0) {
		fprintf(stderr, "[dmap_aln] insufficient memory!\n");
		exit(1);
	}

	fread(reference_info.pac, 1, (reference_info.bns->l_pac) / 4 + 1, reference_info.bns->fp_pac);
	fprintf(stderr, "[dmap_aln] number of characters of reference  %llu \n",
			reference_info.bns->l_pac);

	kmer_position_table = table_restore(tableFn,opt);
	fprintf(stderr, "[dmap_aln] retrieve success! %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);

	mapStart = clock();


	seq = seq_open(readsFile);
	/*Core loop, reads mapping. To save memory, sequentially process chunk of
	 * reads then write results to disk
	 * */
	readsChunk = (read_t*) calloc((opt.chunck_size) + 1, sizeof(read_t));

	while (kseq_read(seq->ks) >= 0) {
		total += 1;
		p = &readsChunk[n_seqs++];
		{ //Read a sequence
			p->name = strdup((const char*) seq->ks->name.s);
			p->len = seq->ks->seq.l;
			readsBp += p->len;

			//copy read sequence
			p->seq = (ubyte_t*) calloc(p->len, 1);
			for (i = 0; i != p->len; i++)
				p->seq[i] = nst_nt4_table[(int) seq->ks->seq.s[i]];

			// copy read quality
			if (seq->ks->qual.l) {
				p->qual = (ubyte_t*) strdup((char*) seq->ks->qual.s);
				p->qual_len = seq->ks->qual.l;
			}
		}

		if (n_seqs == (opt.chunck_size) || readsBp > (opt.chunck_size) * 100) {

			fprintf(stderr, "[dmap_aln] read %d sequences ...", total);
			t = clock();
			ret = lea_map_core(reference_info , kmer_position_table,readsChunk, n_seqs, &opt);
			fprintf(stderr, "%.2f sec\n",
					(float) (clock() - t) / CLOCKS_PER_SEC);
			n_seqs = 0;
			readsBp = 0;
			readsChunk = (read_t*) calloc((opt.chunck_size) + 1,
					sizeof(read_t));
		}
	}
	fprintf(stderr, "[dmap_aln] read %d sequences ...\n", total);
	ret = lea_map_core(reference_info , kmer_position_table,readsChunk, n_seqs, &opt);

	fprintf(stderr, "[dmap_aln] map time %.2f sec\n",
				(float) (clock() - mapStart) / CLOCKS_PER_SEC);



	free(tableFn);
	free(kmer_position_table);

	free(seq);
	return 0;
}
