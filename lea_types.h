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

typedef struct{
	uint32_t distinct_level;
	uint32_t kmer_len;
	bool index_reverse_complement;
}IndexParam;


typedef struct{
	IndexParam  index_parameter;
	uint64_t kmer_table_len;
}Options;


typedef struct{
	uint8_t *pac;
	bntseq_t *bns;
	int64_t l_pac;
	bool contain_reverse_complement;
}ReferenceInfo;

typedef uint32_t TableCell;

#endif /* LEA_TYPES_H_ */
