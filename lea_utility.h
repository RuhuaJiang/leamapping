/*
 * lea_utility.h
 *
 *  Created on: Feb 11, 2014
 *      Author: jiang
 */

#ifndef LEA_UTILITY_H_
#define LEA_UTILITY_H_
char *int2bin(unsigned n, char *buf);
uint8_t read2bits(uint8_t byte, uint8_t pos);
uint32_t list2decimal(std::list<uint8_t> lst );
uint8_t read_position_value(uint8_t *seq, uint64_t pos, int k);
#endif /* LEA_UTILITY_H_ */
