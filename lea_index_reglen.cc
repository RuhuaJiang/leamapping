/*
 * lea_index_reglen.cc
 *
 *  Created on: Feb 11, 2014
 *      Author: jiang
 */
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <list>
#include "utils.h"
#include "lea_types.h"
#include "lea_index_reglen.h"
#include "lea_utility.h"

inline void init_table(const Options opt,TableCell **table)
{
	(*table) =static_cast<TableCell *>(calloc(opt.kmer_table_len,sizeof(TableCell)));
	if((*table) == 0){
		fprintf(stderr,"[LEA_INDXE]Insufficient Memory!\n");
	}
}

//Returns the kmer(seed) length, given the length
//of the reference
inline uint32_t get_kmer_length(int64_t l_pac)
{
	//FIXME need more complex function
	uint32_t kmer_len;
	if(l_pac < 100000 )
		kmer_len = 20;
	else
		kmer_len = 20;
	return kmer_len;
}


/*Returns the distance from start_pos to the next position that the charaters
 *occurs. charaters contain k character
 */
static uint64_t look_ahead(uint64_t start_pos, int charaters, uint8_t *seq, int k){
	uint64_t distance = 2;
	uint8_t key;
	//fprintf(stderr,"k:%d\n",k);
	while(1){
		 // fprintf(stderr,"start_pos+ distance:%d  ",start_pos+ distance);
		  key = read_position_value(seq,start_pos+ distance,k);
		  //fprintf(stderr,"%d\n",key);

		  if(key == charaters ){
			  fprintf(stderr,"dis:%d\n",distance);

			  return distance;
		  }
		  else
			  distance += k;
	 }
}

static void  build_table(ReferenceInfo reference_info,Options opt, TableCell* kmer_position_table)
{
	init_table(opt,&kmer_position_table);
	uint64_t length;
	uint8_t key;
	uint64_t last_char2end_len[4]; // A C G T order

	char debug_array[4] = {'A','C','G','T'};

	//Only index forward sequence
	if (! opt.index_parameter.index_reverse_complement)
	{
	   length = reference_info.contain_reverse_complement ? static_cast<uint64_t>(reference_info.l_pac/2) : static_cast<uint64_t>(reference_info.l_pac);

	   //find last position that char occur, get how far it is from end
	   //FIXME new idea, 2014-2-28, charater can be two base sequence, say "AA"

	   /*
	   for(int charater = 0;charater < 4; charater++){
		  //fprintf(stderr, "working ");
		  uint64_t i = 0;
		  while(1){
			  i ++;
			  key = read_position_value(reference_info.pac, length -1 - i);
			  //fprintf(stderr, "%c",debug_array[key]);
			  if (key == charater){
				last_char2end_len[charater] = i;
				break;
			  }
		  }
	   }// end for
	    */

	   //Find first k number of character
	   int character = 0;
	   uint64_t start_pos=0, distance, debug_total=0;
	   std::vector <uint64_t> distances;
	   std::list <uint8_t> bit_list;
	   for(int i= 0; i < opt.index_parameter.kmer_len; i++){
		   distance = look_ahead(start_pos, character,reference_info.pac, opt.index_parameter.char_size);
		   start_pos = start_pos + distance;
		   distances.push_back(distance);
		   fprintf(stderr, "%d ", distance);
	   }

	   for(int i = 0; i< distances.size()-1; i++){
		   if(distances[i] <= distances[i+1])
			   bit_list.push_back(0);
		   else
			   bit_list.push_back(1);
	   }
	   fprintf(stderr, "\n");
	   //TEST

	   std::list<uint8_t>::iterator  iter;
	   for(iter = bit_list.begin(); iter !=bit_list.end(); iter++ )
	   {
		   fprintf(stderr, "%d", *iter);
	   }

	   fprintf(stderr, "\n");
	   uint64_t last_distance = distance;
	   uint32_t table_key;
	   while(start_pos < length - last_char2end_len[character]){
				 //look_ahead
		   	     if(start_pos % 50000000 == 0) fprintf(stderr, "%llu pos processed\n", start_pos);
				 distance = look_ahead(start_pos , character, reference_info.pac,2);
				 debug_total +=1;
				 start_pos = start_pos + distance;

				 bit_list.pop_front();
				 if( last_distance <= distance)
					 bit_list.push_back(0);
				 else
					 bit_list.push_back(1);
				 table_key = list2decimal(bit_list);
				 kmer_position_table[table_key] +=1;
				 last_distance = distance;
		 }


		#define  LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
		#ifdef LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
		  uint32_t countDistinct =0,countOccur = 0;
		  uint32_t sum = 0;
		  for(uint32_t i = 0; i < opt.kmer_table_len-1; ++i){
			  if(i % 500000000 == 0) fprintf(stderr,"%llu cell processed \n",i);
			  if (kmer_position_table[i] != 0) {
				  countOccur +=1;
				  sum +=kmer_position_table[i];
				  fprintf(stdout, "%d\n",kmer_position_table[i]);
				  }
			  if (kmer_position_table[i] == 1) countDistinct +=1;
		  }//for

		fprintf(stderr,"kmer %d\n",opt.index_parameter.kmer_len);
		fprintf(stderr, "total cell:%llu\n",opt.kmer_table_len-1);
		fprintf(stderr, "total occur %llu \n", debug_total);
		fprintf(stderr, "different occur %llu distinct:%llu   percentage: %f",countOccur, countDistinct,(float)countDistinct/countOccur);
		#endif

	}//end if only index forward sequence
	//fprintf(stderr, "success! \n");
}


void lea_index_region_length(char *indexFile, Options opt){
	clock_t t;
	int forward_only = 0;

	ReferenceInfo reference_info;

	//Parses reference in FASTA file
	gzFile fp = xzopen(indexFile, "r");
	reference_info.l_pac  = bns_fasta2bntseq(fp, indexFile, forward_only); // then contain reverse_complement = true
	reference_info.contain_reverse_complement = true;
	reference_info.bns = bns_restore(indexFile);
	reference_info.pac =  static_cast<uint8_t *>(calloc(reference_info.bns->l_pac / 4 + 1, 1));
	if (reference_info.pac == 0) {
				fprintf(stderr, "[LEA_INDEX] Insufficient Memory!\n");
				exit(1);
	}


	fread(reference_info.pac, 1, reference_info.bns->l_pac / 4 + 1, reference_info.bns->fp_pac);

	//Adjusts parameter based on reference length
	opt.index_parameter.kmer_len = get_kmer_length(reference_info.l_pac);
	opt.kmer_table_len =  pow(2,(opt.index_parameter.kmer_len-1))+10;
	opt.index_parameter.index_reverse_complement = false;
	opt.index_parameter.distinct_level = 1;
	opt.index_parameter.char_size = 2;

	//Builds the kmer position mapping table
	TableCell *kmer_position_table;
	build_table(reference_info,opt,kmer_position_table);

}



