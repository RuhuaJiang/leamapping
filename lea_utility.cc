/*
 * lea_utility.cc
 *
 *  Created on: Feb 11, 2014
 *      Author: jiang
 */


#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <list>
#include <stdio.h>

char *int2bin(unsigned n, char *buf) {
#define BITS (sizeof(n) * CHAR_BIT)

	static char static_buf[BITS + 1];
	int i;

	if (buf == NULL)
		buf = static_buf;

	for (i = BITS - 1; i >= 0; --i) {
		buf[i] = (n & 1) ? '1' : '0';
		n >>= 1;
	}

	buf[BITS] = '\0';
	return buf;

#undef BITS
} ////

// Returns two bits from one byte. The left most index is 0.
// The Possible values of pos are: 0, 2, 4, 6
uint8_t read2bits(uint8_t byte, uint8_t pos) {
	uint8_t twoBits;
	twoBits = ((byte) >> (6 - pos)) & ((1 << (2)) - 1);
	//fprintf(stderr, " byte %s ",int2bin(byte, NULL));
	return twoBits;
}


//Returns value at given position of given 2-bit per char presentation of sequence
uint8_t read_position_value_helper(uint8_t *seq, uint64_t pos){
	uint32_t index;
	uint8_t shift,twoBits;
	index = pos / 4 ;
	shift = (pos % 4)*2;
	//fprintf(stderr,"index:%d shift:%d",index,shift);
	twoBits= read2bits(seq[index], shift);
	return twoBits;
}
//Returns value at given position, read k charaters. Each character is in 2-bit reprensentation. k is max 4 here
uint8_t read_position_value(uint8_t *seq, uint64_t pos, int k){
	uint8_t twoBits = 0;
	if(k > 4 || k<1)
	{
		fprintf(stderr,"k not proper, should [1, 4] \n");
		 exit(1);
	}
	//fprintf(stderr,"pos%llu ",pos);
	for(int i = 0; i<k; i++){
		twoBits <<=2;
		twoBits +=read_position_value_helper(seq,pos+i);
		//fprintf(stderr,"twoBits%d ",twoBits);
	}
	return twoBits;
}

uint32_t list2decimal(std::list<uint8_t> lst ){
	uint32_t result = 0;
	std::list<uint8_t>::iterator  iter;
	for(iter = lst.begin(); iter != lst.end(); iter++ )
	{
		result <<=1;
		result += *iter;
	}
	return result;
}
