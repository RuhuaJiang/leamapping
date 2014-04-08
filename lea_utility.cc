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
#include <vector>
#include <stdio.h>
#include <utility>
#include <algorithm>
#include <map>

#include "lea_types.h"
#include "lea_utility.h"

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
			twoBits += read_position_value_helper(seq,pos+i);
		else
		{
			twoBits += seq[pos+i];
			//fprintf(stderr,"pos %llu, twoBits%d ", pos+i, seq[pos+i]);
			//if(pos >600 )exit(1);
		}

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
		kmer_len = 7;
	else
		kmer_len = 7;
	return kmer_len;
}


/*Returns the distance from start_pos to the next position that the charaters
 *occurs. charaters contain k character
 */
uint64_t look_ahead(uint64_t start_pos, uint32_t charaters, uint8_t *seq, int k, int bit_byte){
	uint64_t distance = 1;
	uint8_t key;
	//fprintf(stderr,"charaters: %d\n",charaters);
	while(1){
		  // fprintf(stderr,"start_pos+ distance:%d  ",start_pos+ distance);
		  key = read_position_value(seq,start_pos+distance,k, bit_byte);
		 // fprintf(stderr,"key: %d\n",key);
		  if(static_cast<uint32_t>(key) == charaters ){
			  //fprintf(stderr,"dis:%d\n",distance);
			  return distance;
		  }
		  else
			  distance++;
	 }
}


/*
 * Returns the distance from start start_pos to the next position that 'island'(at least 2 consecutive charater)
 */
uint64_t look_ahead_island(uint64_t start_pos, uint32_t character, uint8_t *seq, uint32_t seq_len,int bit_byte)
{
	 uint64_t distance = 0;
	 uint8_t key;
	 int state= 0;// state = 0, not character state;  state = 1, character occurs once; state = 2, character occurs at least twice;
	 //state machine
	 //fprintf(stderr,"character: %d\n",character);
	 while(1)
	 {
		if(start_pos+distance < seq_len)
		{
			key = read_position_value(seq,start_pos+distance,1, bit_byte);

			distance++;
			//fprintf(stderr,"key: %d %d\n",key, distance);
			if(static_cast<uint32_t>(key) == character  && state == 0) //1
			{
				state = 1;
				continue;
			}
			if(static_cast<uint32_t>(key) == character  && state == 1)  //2
			{
				 state = 2;
				 continue;
			}
			if(static_cast<uint32_t>(key) != character  && state == 1) //3
			{
				  state = 0;
				  continue;
			}
			if(static_cast<uint32_t>(key) != character  && state == 0) //4
			{
				  state = 0;
				  continue;
			}
			if(static_cast<uint32_t>(key) != character  && state == 2) //5
			{
				  state = 0;
				  break;
			}
			if(static_cast<uint32_t>(key) == character  && state == 2) //6
			{
				  state = 2;
				  continue;
			}
		}
		else break;

	 }//while
	 return distance -1;
}

uint64_t look_ahead_dense_island(uint64_t start_pos, uint32_t character, uint8_t *seq, uint32_t seq_len,int bit_byte)
{
	 uint64_t distance = 0;
	 uint8_t key;
	 int state= 0;// state = 0, not character state;  state = 1, character occurs once; state = 2, character occurs at least twice;
	 //state machine
	 //fprintf(stderr,"character: %d\n",character);
	 while(1)
	 {
		if(start_pos+distance < seq_len)
		{
			key = read_position_value(seq,start_pos+distance,1, bit_byte);

			distance++;
			//fprintf(stderr,"key: %d %d\n",key, distance);
			if(static_cast<uint32_t>(key) == character  && state == 0) //1
			{
				state = 1;
				continue;
			}
			if(static_cast<uint32_t>(key) == character  && state == 1)  //2
			{
				 state = 2;
				 continue;
			}
			if(static_cast<uint32_t>(key) != character  && state == 1) //3
			{
				  state = 0;
				  continue;
			}
			if(static_cast<uint32_t>(key) != character  && state == 0) //4
			{
				  state = 0;
				  continue;
			}
			if(static_cast<uint32_t>(key) != character  && state == 2) //5
			{
				  state = 0;
				  break;
			}
			if(static_cast<uint32_t>(key) == character  && state == 2) //6
			{
				  state = 2;
				  continue;
			}
		}
		else break;

	 }//while

}



void get_spectrum(uint64_t start_pos,uint64_t length,uint8_t *seq,uint32_t seq_len,int bit_byte, Spectrum &spec)
{
	uint8_t key;
	for(int i=0; i<4; i++)
		spec.spec[i]=0;

	for(uint64_t i = 0; i < length; i++)
	{
		key = read_position_value(seq,start_pos+i,1, bit_byte);
		spec.spec[key] +=1;
	}
	/*
	for(int i=0; i<4; i++)
	{
		fprintf(stderr,"%d ", spec.spec[i]);
	    // spec.spec[i]= spec.spec[i]/spec.round_num;
	}
	*/


}
//Returns the max min --Idea of 2014-4-3
uint64_t process_spectrum(Spectrum spec)
{
  uint32_t maxv=spec.spec[0], minv = spec.spec[0];
  uint64_t maxi=0, mini=0;
  uint64_t result = 0;
  for(int i=1; i<4; i++){
	  if(spec.spec[i] > maxv)
	  {
		  maxv = spec.spec[i];
		  maxi = i;
	  }
  }
  for(int i=1; i<4; i++){
	  if(spec.spec[i] < minv)
  	  {
  		  minv = spec.spec[i];
  		  mini = i;
  	  }
    }
  result += maxi;
  result <<=2;
  result += mini;
  //fprintf(stderr,"max_min: %d %d\n", maxi,mini);
  return result;
}
uint32_t vector_to_int_noround(std::list <uint64_t> distances)
{

	    if(distances.size() > 8){
			fprintf(stderr,"list to large, should smaller than or equal to 8\n");
			exit(1);
		}
		 std::list<uint64_t>::iterator  iter;
		 uint64_t val;
		 uint32_t result = 0;
		 for(iter = distances.begin(); iter !=distances.end(); iter++ )
		 {
			 //fprintf(stderr,"%llu ",(*iter));
			 result <<=4;
			 val = (*iter);
			 result += val;
		 }
		 return result;

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
		 val = (*iter) > static_cast<uint64_t>(45) ? static_cast<uint64_t>(15):((*iter)/3);
		 result += val;
	 }
	 return result;
}

void supported_positions_chainning(std::vector<PositionShift> positions_shifts, uint32_t len)
{
	PositionShift position_shift;
	std::vector<std::pair<int64_t,int64_t> > posFreqDic;

	int64_t position;
	std::map<int64_t,int64_t> position_freq;
	for(int i=0; i<positions_shifts.size(); i++){
		position_shift = positions_shifts[i];
		position = (position_shift.position - static_cast<int64_t>(position_shift.shift))*position_shift.strand;
		position_freq[position/1000] = 0;
	}
	for(int i=0; i<positions_shifts.size(); i++){
		position_shift = positions_shifts[i];
		position = (position_shift.position - static_cast<int64_t>(position_shift.shift))*position_shift.strand;
		position_freq[position/1000] +=1;
	}
	std::pair<int64_t,int64_t>  _pair;
	for (std::map<int64_t,int64_t>::iterator it=position_freq.begin(); it!=position_freq.end(); ++it){
		_pair.first = it->first;
		_pair.second = it->second;
		posFreqDic.push_back(_pair);
	}
	std::sort(posFreqDic.begin(), posFreqDic.end(),&compareParSecondDec);
	if(posFreqDic.size() !=0)
	{
		if(posFreqDic[0].second >0)
		{
			//for(int i=0; i<posFreqDic.size();++i)
			//{
				fprintf(stdout,"%lld_%llu ",posFreqDic[0].first,posFreqDic[0].second);
			//}
			//fprintf(stderr,"\n");
		}
	}


	fprintf(stdout,"\n");
}

bool compareParSecondDec(std::pair<int64_t,int64_t> s1, std::pair<int64_t,int64_t> s2)
{

	 if (s1.second <= s2.second) return false;
	 if (s1.second > s2.second) return true;

}


