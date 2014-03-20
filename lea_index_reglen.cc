/*
 * lea_index_reglen.cc
 *
 *  Created on: Feb 11, 2014
 *      Author: jiang
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
#include <time.h>
#include <string.h>
#include "utils.h"
#include "lea_types.h"
#include "lea_index_reglen.h"
#include "lea_utility.h"
#include "lea_seqio.h"

inline void init_table(const Options opt,TableCell **table)
{
	(*table) =static_cast<TableCell *>(calloc(opt.kmer_table_len,sizeof(TableCell)));
	if((*table) == 0){
		fprintf(stderr,"[LEA_INDXE]Insufficient Memory!\n");
	}
}


static void  build_table(ReferenceInfo reference_info,Options opt, TableCell* kmer_position_table)
{
   init_table(opt,&kmer_position_table);
   uint64_t length;
   uint8_t key;
   uint64_t last_char2end_len[4]; // A C G T order
   clock_t t;

	//Only index forward sequence


   length = reference_info.contain_reverse_complement ? static_cast<uint64_t>(reference_info.l_pac/2) : static_cast<uint64_t>(reference_info.l_pac);

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
   int character = opt.index_parameter.char_int;
   std::list <uint64_t> distances;
   int step = opt.index_parameter.char_size;
   uint64_t start_pos=0, distance, debug_total=0;
   for(int i= 0; i < opt.index_parameter.kmer_len; ++i){
	   //distance = look_ahead(start_pos, character,reference_info.pac, opt.index_parameter.char_size,0);
	   distance =  distance = look_ahead_island(start_pos,character,reference_info.pac,reference_info.l_pac,0);

	   start_pos = start_pos + distance;
	   distances.push_back(distance);
	   fprintf(stderr, "%d ", distance);

   }
   fprintf(stderr, "\n");


   std::list<uint64_t>::iterator  iter;

   uint64_t last_distance = distance, test_sum = 0, test_count =0 ;
   uint32_t table_key;


   while(start_pos < length - 20){
			 //look_ahead
			 if(start_pos % 50000000 == 0) fprintf(stderr, "%llu pos processed\n", start_pos);
			 //distance = look_ahead(start_pos , character, reference_info.pac,opt.index_parameter.char_size,0);
			 distance = look_ahead_island(start_pos,character,reference_info.pac,reference_info.l_pac,0);
			 test_count +=1;
			 test_sum +=distance;
			 start_pos = start_pos + distance;

			 //fprintf(stderr, "%d ", start_pos);
			 distances.push_back(distance);
			 distances.pop_front();
			 /*
			 for(iter = distances.begin(); iter !=distances.end(); iter++ )
			 {
			   fprintf(stderr, "%d ", *iter);
			 }
			 fprintf(stderr, "\n");
			 */

			 vector_to_int(distances);

			 table_key = vector_to_int(distances);
			 //fprintf(stderr, "%d ", *iter);

			 //kmer_position_table[table_key] +=1;


			 if ( kmer_position_table[table_key] == 0 ) // state 0 -> state 1
				  kmer_position_table[table_key] =  static_cast<uint32_t>(start_pos);
			 else if (kmer_position_table[table_key] == NOT_UNIQUE)  // state 2 -> state 2
			 {
			 }
			 else if(kmer_position_table[table_key] != 0 && kmer_position_table[table_key] != NOT_UNIQUE)  // state 1 -> state 2
			 {
				 kmer_position_table[table_key] = NOT_UNIQUE;
			 }
			 else
			 {
				 fprintf(stderr,"[Index Ref] ERRO! \n");
				 exit(1);
			 }

	 }
    fprintf(stderr, "testsum:%llu % llu %d\n",test_sum, test_count,test_sum/test_count);
    //exit(1);
	//dump
	t = clock();
	char *tableFn;
	tableFn = (char*)calloc(strlen(reference_info.ref_name) + 10, 1);
	fprintf(stderr,"%s",reference_info.ref_name);
	char *name_suffix;
	name_suffix = (char *) malloc(sizeof(char) * (10));
	sprintf(name_suffix, "%d", opt.index_parameter.char_int);
	strcpy(tableFn , reference_info.ref_name);strcat(tableFn,".tb");strcat(tableFn,name_suffix);


	fprintf(stderr, "[dmap index] write index results to file...\n" );
	if (table_dump(tableFn,kmer_position_table, opt) == 0)
	{
		fprintf(stderr, "[dmap index]  write success! %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	#define  LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
	#ifdef LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
	  uint32_t countDistinct =0,countOccur = 0;
	  uint32_t sum = 0;
	  for(uint64_t i = 0; i < static_cast<uint64_t>(opt.kmer_table_len-1); ++i){
		  if(i % 500000000 == 0) fprintf(stderr,"%llu cell processed \n",i);
		  if (kmer_position_table[i] != 0) {
			  countOccur +=1;
			  sum +=kmer_position_table[i];
			  //fprintf(stdout, "%d\n",kmer_position_table[i]);
			  }
		  if (kmer_position_table[i] != 0 && kmer_position_table[i] != NOT_UNIQUE) countDistinct +=1;
	  }//for
	fprintf(stderr,"reference_info.l_pac %d\n",reference_info.l_pac);
	fprintf(stderr,"kmer %d\n",opt.index_parameter.kmer_len);
	fprintf(stderr, "total cell:%llu\n",opt.kmer_table_len-1);
	fprintf(stderr, "total occur %llu \n", debug_total);
	fprintf(stderr, "different occur %llu distinct:%llu   percentage: %f",countOccur, countDistinct,(float)countDistinct/countOccur);
	#endif


}


void lea_index_region_length(char *indexFile, Options opt){
	clock_t t;
	int forward_only = 0;

	ReferenceInfo reference_info;

	//Parses reference in FASTA file
	gzFile fp = xzopen(indexFile, "r");
	reference_info.ref_name = indexFile;
	reference_info.l_pac  = bns_fasta2bntseq(fp, indexFile, forward_only); //if forward_only =0, then contain reverse_complement = true
	reference_info.contain_reverse_complement = true;
	reference_info.bns = bns_restore(indexFile);
	reference_info.pac =  static_cast<uint8_t *>(calloc(reference_info.bns->l_pac / 4 + 1, 1));
	if (reference_info.pac == 0) {
				fprintf(stderr, "[LEA_INDEX] Insufficient Memory!\n");
				exit(1);
	}


	fread(reference_info.pac, 1, reference_info.bns->l_pac/4 + 1, reference_info.bns->fp_pac);

	//Adjusts parameter based on reference length
	opt.index_parameter.kmer_len = get_kmer_length(reference_info.l_pac);
	opt.kmer_table_len =  pow(2,(opt.index_parameter.kmer_len * 4))+10;
	opt.index_parameter.index_reverse_complement = false;
	opt.index_parameter.distinct_level = 1;
	opt.index_parameter.char_size = 2;



	//Builds the kmer position mapping table
	TableCell *kmer_position_table_AA,*kmer_position_table_CC,*kmer_position_table_GG,*kmer_position_table_TT;
	opt.index_parameter.char_int =0;
	build_table(reference_info,opt,kmer_position_table_AA);

	opt.index_parameter.char_int =1;
	build_table(reference_info,opt,kmer_position_table_CC);
	opt.index_parameter.char_int =2;
	build_table(reference_info,opt,kmer_position_table_GG);
	opt.index_parameter.char_int =3;
	build_table(reference_info,opt,kmer_position_table_TT);


}



