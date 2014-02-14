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
		kmer_len = 31;
	return kmer_len;
}

static uint64_t look_ahead(uint64_t start_pos, int charater, uint8_t *seq){
	uint64_t distance = 1;
	uint8_t key;
	while(1){
		  key = read_position_value(seq,start_pos + distance);
		  if(key == charater ){
			  return distance;
		  }
		  else
			  distance +=1;
	 }
}

static void  build_table(ReferenceInfo reference_info,Options opt, TableCell* kmer_position_table)
{
	init_table(opt,&kmer_position_table);
	uint64_t pos, length;
	uint8_t key;
	uint32_t last_char2end_len[4]; // A C G T order

	char debug_array[4] = {'A','C','G','T'};

	//Only index forward sequence
	if (! opt.index_parameter.index_reverse_complement)
	{
	  length = reference_info.contain_reverse_complement ? static_cast<uint64_t>(reference_info.l_pac/2) : static_cast<uint64_t>(reference_info.l_pac);

	  //find last position that char occur, get how far it is from end
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


	  {
		  //FIXME actually we need build 4 table, right now only consider the table using A
		  //for each position, generate sliced kmers
		  int charater = 0;
		  uint64_t shift =0, distance = 0;
		  std::list<uint8_t> last_distances_bits;
		  std::vector<uint64_t> distances;
		  {
			  //Find firt k number of charater
			  int k = opt.index_parameter.kmer_len;
			  while(1){
				  shift++;
				  distance ++;
				  key = read_position_value(reference_info.pac, shift);
				  //fprintf(stderr, "k:%llu ",key);
				  if (key == charater){
					distances.push_back(distance-1);

					//fprintf(stderr, "d:%llu\n",distance);
					distance = 0;
					k--;
				  }
				  if(k == 0)break;
			  }
			  for(uint i=0;i< distances.size() -1;i++){
				  if( distances[i] <= distances[i+1])
					  last_distances_bits.push_back(0);
				  else
					  last_distances_bits.push_back(1);
			  }
		  }

		  std::list<uint8_t>::iterator iter;
		  for(iter = last_distances_bits.begin(); iter !=last_distances_bits.end(); iter++ )
		  {
		  	fprintf(stderr, "%d ", *iter);
		  }
		  //fprintf(stderr, "%d",distances[distances.size()-1]);
		  fprintf(stderr, "\n");

		  fprintf(stderr, "start... k %llu table_len: %llu\n",opt.index_parameter.kmer_len, opt.kmer_table_len);

		  DistancesInfo distances_info;
		  uint32_t table_key ,debug_total=0;
		  pos = shift+1;
		  uint64_t last_distance = distances[distances.size()-1];
		  while(pos < 1000/*length - last_char2end_len[charater]*/){
			  //look_ahead
			  distance = look_ahead(pos , charater, reference_info.pac);
			  debug_total +=1;
			  pos = pos + distance;
			  std::list<uint8_t>::iterator  iter;
			  for(iter = last_distances_bits.begin(); iter !=last_distances_bits.end(); iter++ )
			  {
			  		fprintf(stderr, "%d", *iter);
			  }
			  fprintf(stderr, "\nafter:\n");
			  //then remove first
			  last_distances_bits.pop_front();
			  for(iter = last_distances_bits.begin(); iter !=last_distances_bits.end(); iter++ )
			  {
			  		fprintf(stderr, "%d", *iter);
			  }
			  fprintf(stderr, "\n");
			  if(last_distance <= distance ) last_distances_bits.push_back(0);
			  else last_distances_bits.push_back(1);
			  fprintf(stderr, "%llu", distance);
			  last_distance = distance;

			  table_key = list2decimal(distances_info.disances_bits);
			  kmer_position_table[table_key] +=1;

			  fprintf(stderr, "\n");

		  }


		#define  LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
		#ifdef LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
		  uint32_t countDistinct =0,countOccur = 0;
		  uint32_t sum = 0;
		  for(uint32_t i = 0; i < opt.kmer_table_len-1; ++i){
			  if(i % 50000000 == 0) fprintf(stderr,"%llu cell processed \n",i);
			  if (kmer_position_table[i] != 0) {
				  countOccur +=1;
				  sum +=kmer_position_table[i];
				  //fprintf(stdout, "%s   %llu  %llu \n",int2bin(i, NULL), i,kmer_position_table[i]);
				  }
			  if (kmer_position_table[i] == 1) countDistinct +=1;

		  }//for

		  fprintf(stderr,"kmer_position_table[395554] %llu  \n",kmer_position_table[395554]);
		  fprintf(stderr, "total cell:%llu\n",opt.kmer_table_len-1);
		  fprintf(stderr, "total occur:%llu   %llu\n",debug_total,sum);
		  fprintf(stderr, "occur %llu distinct:%llu   percentage: %f",countOccur, countDistinct,(float)countDistinct/countOccur);
		#endif
		}
	}//end if
	//fprintf(stderr, "success! \n");
}


void lea_index_region_length(char *indexFile, Options opt){
	clock_t t;
	int forward_only = 0;

	ReferenceInfo reference_info;

	//Parses reference in FASTA file
	gzFile fp = xzopen(indexFile, "r");
	reference_info.l_pac  = bns_fasta2bntseq(fp, indexFile, forward_only);
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

	//Builds the kmer position mapping table
	TableCell *kmer_position_table;
	build_table(reference_info,opt,kmer_position_table);

}



