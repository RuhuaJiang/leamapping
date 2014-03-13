/*
 * lea_seqio.h
 *
 *  Created on: Mar 12, 2014
 *      Author: jiang
 */


int table_dump(char *fn,const TableCell *kmer_position_table, Options opt);
TableCell *table_restore(const char *fn,Options opt);
seqio_t *seq_open(const char *fn);
void free_read_seq(int n_seqs, read_t *seqs);
