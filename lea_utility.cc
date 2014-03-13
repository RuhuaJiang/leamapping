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
uint8_t read_position_value(uint8_t *seq, uint64_t pos, int k, int bit_byte){
	uint8_t twoBits = 0;
	if(k > 4 || k<1)
	{
		fprintf(stderr,"k not proper, should [1, 4] \n");
		exit(1);
	}
	//fprintf(stderr,"pos%llu ",pos);
	for(int i = 0; i<k; i++){
		twoBits <<=2;
		//fprintf(stderr,"pos+i%d ",pos+i);
		if(bit_byte == 0)
			twoBits +=read_position_value_helper(seq,pos+i);
		else
			twoBits += seq[pos+i];
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


//Returns the kmer(seed) length, given the length
//of the reference
uint32_t get_kmer_length(int64_t l_pac)
{
	//FIXME need more complex function
	uint32_t kmer_len;
	if(l_pac < 100000 )
		kmer_len = 5;
	else
		kmer_len = 8;
	return kmer_len;
}


/*Returns the distance from start_pos to the next position that the charaters
 *occurs. charaters contain k character
 */
uint64_t look_ahead(uint64_t start_pos, int charaters, uint8_t *seq, int k, int bit_byte){
	uint64_t distance = 1;
	uint8_t key;

	while(1){
		  // fprintf(stderr,"start_pos+ distance:%d  ",start_pos+ distance);
		  key = read_position_value(seq,start_pos+distance,k, bit_byte);
		  //fprintf(stderr,"key: %d\n",key);
		  if(key == charaters ){
			  //fprintf(stderr,"dis:%d\n",distance);
			  return distance;
		  }
		  else
			  distance++;
	 }
}


//Returns a integer. Input is a size k list;
uint32_t vector_to_int(std::list <uint64_t> distances){

	if(distances.size() > 8){
		fprintf(stderr,"list to large, should smaller than or equal to 8\n");
		exit(1);
	}
	 std::list<uint64_t>::iterator  iter;
	 uint64_t val;
	 uint32_t result = 0;
	 for(iter = distances.begin(); iter !=distances.end(); iter++ )
	 {
		 result <<=4;
		 val = (*iter) > static_cast<uint64_t>(31) ? static_cast<uint64_t>(15):((*iter)/2);
		 result += val;
	 }
	 return result;
}



