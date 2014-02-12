/*
 * lea_index_reglen.h
 *
 *  Created on: Feb 11, 2014
 *      Author: jiang
 */

#ifndef LEA_INDEX_REGLEN_H_
#define LEA_INDEX_REGLEN_H_
#include<list>
#include<stdint.h>
typedef struct{
	std::list<uint8_t> disances_bits;
	uint64_t last_distance;
}DistancesInfo;


void lea_index_region_length(char *indexFile, Options opt);

#endif /* LEA_INDEX_REGLEN_H_ */
