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
		kmer_len = 2;
	else
		kmer_len = 30;
	return kmer_len;
}
//0 3 0 1 2 0 0 0 0 0 11 2 0
//Returns a number,
static DistancesInfo get_increase_decrease(DistancesInfo distances_info, uint64_t distance){

	DistancesInfo dis_info;
	dis_info = distances_info;
	dis_info.disances_bits.pop_front();
	if( distances_info.last_distance <= distance)
		dis_info.disances_bits.push_back(0);
	else
		dis_info.disances_bits.push_back(1);
	dis_info.last_distance = distance;
	return dis_info;
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

		  DistancesInfo distances_info;
		  distances_info.last_distance = distances[distances.size()-1];
		  distances_info.disances_bits =last_distances_bits;
		  for(int i=0;i< distances.size() ;i++) fprintf(stderr, "%d ",distances[i]);
		  //fprintf(stderr, "last: %llu ",  distances_info.last_distance);
		  int _distance = 0;
		  for(pos = shift+1; pos < 1000/*length - last_char2end_len[charater]*/; pos++){
			key = read_position_value(reference_info.pac,pos);
			fprintf(stderr, "%d\n",key);
			if (key == charater){
				distances_info = get_increase_decrease(distances_info,_distance);
				//fprintf(stderr, "%d\n",_distance);
				_distance = 0;
			}
			_distance++;
		  }
		//#define  LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
		#ifdef LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG

		#endif
		}
	}//end if
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
	opt.kmer_table_len =  pow(2,opt.index_parameter.kmer_len);
	opt.index_parameter.index_reverse_complement = false;
	opt.index_parameter.distinct_level = 1;

	//Builds the kmer position mapping table
	TableCell *kmer_position_table;
	build_table(reference_info,opt,kmer_position_table);

}



