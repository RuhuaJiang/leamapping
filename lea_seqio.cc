/*
 * lea_seqio.cc

 *
 *  Created on: Mar 12, 2014
 *      Author: jiang
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "utils.h"
#include "lea_types.h"
#include "kseq.h"

void free_read_seq(int n_seqs, read_t *seqs)
{
	int i;
	for (i = 0; i != n_seqs; ++i) {
		read_t *p = seqs + i;
		free(p->name);
		free(p->seq);  free(p->qual);
		//free(p->cigar);
	}
	free(seqs);
}

int table_dump(char *fn,const TableCell *kmer_position_table, Options opt)
{
	FILE *fp;

	fp = xopen(fn, "wb");
	fwrite(kmer_position_table,sizeof(TableCell),opt.kmer_table_len,fp );
	fclose(fp);
	return 0;


}

TableCell *table_restore(const char *fn,Options opt)
{
	TableCell *kmerPositionTable;
	FILE *fp;

	kmerPositionTable = (TableCell*)calloc(opt.kmer_table_len, sizeof(TableCell));
	if (kmerPositionTable == 0) {
			fprintf(stderr, "[dmap_aln] insufficient memory!\n");
			exit(1);
	}
	fp = xopen(fn, "rb");
	fseek(fp, 0, SEEK_SET);
	fread(kmerPositionTable, sizeof(TableCell), opt.kmer_table_len, fp);
	fclose(fp);
	return kmerPositionTable;
}

seqio_t *seq_open(const char *fn)
{
	gzFile fp;
	seqio_t *bs;
	bs = (seqio_t*)calloc(1, sizeof(seqio_t));
	fp = xzopen(fn, "r");
	bs->ks = kseq_init(fp);
	return bs;
}


