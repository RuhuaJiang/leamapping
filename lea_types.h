/*
 * lea_types.h
 *
 *  Created on: Feb 7, 2014
 *      Author: jiang
 */

#ifndef LEA_TYPES_H_
#define LEA_TYPES_H_
#include <stdint.h>
#include "bntseq.h"
#include "utils.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
typedef struct{
	uint32_t distinct_level;
	uint32_t kmer_len;
	uint32_t char_size; // regard char_size number of character as a single character when index
	uint32_t char_int;
	bool index_reverse_complement;
}IndexParam;


typedef struct{
	IndexParam  index_parameter;
	uint64_t kmer_table_len;

	uint32_t chunck_size;

	int a; //Score of a match [1]
	int b; //Mismatch penalty [3]
	int q; //Gap open penalty [5]
	int r; //Gap extension penalty. The penalty for a contiguous gap of size k is q+k*r. [2]
	int bw; //-w INT	 Band width in the banded alignment [33]
}Options;


typedef struct{
	uint8_t *pac;
	bntseq_t *bns;
	int64_t l_pac;
	char *ref_name;
	bool contain_reverse_complement;
}ReferenceInfo;

typedef uint16_t cigar_t;
typedef struct{
	char *name;
	ubyte_t *seq,*qual;
	uint32_t aln_pos;
	cigar_t cigar;
	uint32_t len;
	uint32_t qual_len;
    /*
	int32_t *shifts; //Positions that has distinct Kmer
	int32_t *positions; //Kmer position in reference
	uint32_t kmerN; //number of distinct kmer.  Always < =  opt->firstZ
	*/

}read_t;

typedef struct{
	// for fastq input
	kseq_t *ks;
}seqio_t;
const uint32_t NOT_UNIQUE = (uint32_t)0xFFFFFFFF -1;
const int ARRAY_LEN = 1000 ;
typedef uint32_t TableCell;

typedef struct{
	int64_t position;
	uint32_t shift;
}PositionShift;

typedef struct{
	TableCell *kmer_position_table_AA;
	TableCell *kmer_position_table_CC;
	TableCell *kmer_position_table_GG;
	TableCell *kmer_position_table_TT;
}Tables;
typedef struct{
	 int64_t first;
	 int64_t second;
}pair;

#endif /* LEA_TYPES_H_ */
