/*
 * lea_utility.h
 *
 *  Created on: Feb 11, 2014
 *      Author: jiang
 */

#ifndef LEA_UTILITY_H_
#define LEA_UTILITY_H_
#include <list>

char *int2bin(unsigned n, char *buf);
uint8_t read2bits(uint8_t byte, uint8_t pos);
uint32_t list2decimal(std::list<uint8_t> lst );
uint8_t read_position_value(uint8_t *seq, uint64_t pos, int k, int bit_byte);
int getKmerStrandIndicatorBit(uint32_t kmerMapArrayPos,const uint8_t *kmerStrandIndicator);
int setKmerStrandIndicatorBit(uint32_t kmerMapArrayPos,uint8_t *kmerStrandIndicator);
uint32_t get_kmer_length(int64_t l_pac);
uint64_t look_ahead(uint64_t start_pos, int charaters, uint8_t *seq, int k, int bit_byte);
uint32_t  vector_to_int(std::list <uint64_t> distances);
bool compareParSecondDec(std::pair<int64_t,int64_t> s1, std::pair<int64_t,int64_t> s2);
#endif /* LEA_UTILITY_H_ */
