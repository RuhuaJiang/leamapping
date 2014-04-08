/*
 * lea_aln.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: jiang
 */
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>
#include "lea_types.h"
#include "lea_seqio.h"
#include "lea_utility.h"

#include <inttypes.h>

static int count_supported_positions_shifts(ReferenceInfo reference_info , TableCell *kmer_position_table,uint8_t *read, uint32_t read_len,Options opt,bool is_rvc, std::vector<PositionShift> &positions_shifts){

	uint64_t last_distance = opt.index_parameter.char_size;
	uint8_t key;
	uint32_t character = opt.index_parameter.char_int;


	while(1){
		      key =read_position_value(read,read_len-last_distance,opt.index_parameter.char_size,1);
			  if(static_cast<uint32_t>(key) == character ){
				  break;
			  }
			  else
				  last_distance++;
	}

    std::list <uint64_t> distances;
    std::list <uint64_t> max_mins;
    int step = opt.index_parameter.char_size;
    uint64_t start_pos=0, distance, debug_total=0;
    /*
    for(int i=0;i< read_len; i++)
    	fprintf(stderr,"%d",read[i]);
    fprintf(stderr,"%\n");
    */

    for(int i= 0; i < opt.index_parameter.kmer_len; ++i){
	   //distance = look_ahead(start_pos, character,read, opt.index_parameter.char_size,1);
	   distance = distance = look_ahead_island(start_pos,character,read,read_len,1);
	   {
	  	   Spectrum spec;
	  	   uint64_t max_min;
	  	   get_spectrum(start_pos,distance,read,read_len,1,spec);
	  	   max_min = process_spectrum(spec);
	  	   max_mins.push_back(max_min);
	   }
	   start_pos = start_pos + distance;
	   distances.push_back(distance);
	   //fprintf(stderr, "%d ", distance);
    }

    std::list<uint64_t>::iterator  iter;
    /*
    fprintf(stderr, "start:\n");
    for(iter = distances.begin(); iter !=distances.end(); iter++ )
	{
	   fprintf(stderr, "%d ", *iter);
	}
	fprintf(stderr, "\n"); */

	uint64_t test_sum = 0, test_count =0 ;
	uint32_t table_key, position;

	PositionShift _postion_shift;
	//fprintf(stderr,"%llu, %llu",start_pos,static_cast<uint64_t>(read_len) - last_distance );
	//fprintf(stderr,"int:%d\n",character);
	while(start_pos < static_cast<uint64_t>(read_len) - 100 ){
		//distance = look_ahead(start_pos , character, read,opt.index_parameter.char_size,1);
		distance = look_ahead_island(start_pos,character,read,read_len,1);
		{
			 Spectrum spec;
			 uint64_t max_min;
			 get_spectrum(start_pos,distance,read,read_len,1, spec);
			 max_min = process_spectrum(spec);
			 max_mins.push_back(max_min);
			 max_mins.pop_front();
			 table_key = vector_to_int_noround(max_mins);
			 //fprintf(stderr, "key: %llu  read_len%llu\n", table_key,read_len);
			 //exit(1);
		}
		start_pos = start_pos + distance;

		#ifdef USE_DISTANCE
		{
			distances.push_back(distance);
			distances.pop_front();
			table_key = vector_to_int(distances);
		}
		#endif

		/*
		for(iter = distances.begin(); iter !=distances.end(); iter++ )
		{
		   fprintf(stderr, "%d ", *iter);
		}
		fprintf(stderr, "\n");
		*/
		//exit(1);

		position = kmer_position_table[table_key];
		//fprintf(stderr,"%llu\t%llu\t", table_key,position);
		if(position !=0 && position != NOT_UNIQUE)
		{
			//fprintf(stderr,"distinct\n");
			if(is_rvc)
				_postion_shift.strand = -1;
			else
				_postion_shift.strand = 1;
			_postion_shift.position = static_cast<int64_t>(position);
			_postion_shift.shift = static_cast<uint64_t>(start_pos);
			//fprintf(stderr,"%llu,%llu ",table_key,_postion_shift.position);
			//fprintf(stderr,"%llu",_postion_shift.shift);
			positions_shifts.push_back(_postion_shift);
			//exit(1);
		}
	}

	//fprintf(stderr,"\n");
	return 0;
}

static int lea_map_single_read(ReferenceInfo reference_info , Tables kmer_position_tables,read_t *read,Options opt){

	uint8_t *rvc_read_seq;
	std::vector<PositionShift> positions_shifts;
	std::vector<std::pair<int64_t,int64_t> > posFreqDic;
	uint32_t posFreqDicCount = 0;
	rvc_read_seq = (uint8_t*)calloc(read->len, 1);
	for (int _i = 0; _i < read->len; ++_i)
			rvc_read_seq[read->len - _i - 1] = 3 - read->seq[_i];

	opt.index_parameter.char_int = 0;
	count_supported_positions_shifts(reference_info , kmer_position_tables.kmer_position_table_AA,read->seq,read->len,opt, false, positions_shifts);
	count_supported_positions_shifts(reference_info , kmer_position_tables.kmer_position_table_AA,rvc_read_seq,read->len,opt, true, positions_shifts);


	opt.index_parameter.char_int = 1;
	count_supported_positions_shifts(reference_info , kmer_position_tables.kmer_position_table_CC,read->seq,read->len,opt, false, positions_shifts);
	count_supported_positions_shifts(reference_info , kmer_position_tables.kmer_position_table_CC,rvc_read_seq,read->len,opt, true, positions_shifts);


	opt.index_parameter.char_int = 2;
	count_supported_positions_shifts(reference_info , kmer_position_tables.kmer_position_table_GG,read->seq,read->len,opt, false, positions_shifts);
	count_supported_positions_shifts(reference_info , kmer_position_tables.kmer_position_table_GG,rvc_read_seq,read->len,opt, true, positions_shifts);

	opt.index_parameter.char_int = 3;
	count_supported_positions_shifts(reference_info , kmer_position_tables.kmer_position_table_TT,read->seq,read->len,opt, false, positions_shifts);
	count_supported_positions_shifts(reference_info , kmer_position_tables.kmer_position_table_TT,rvc_read_seq,read->len,opt, true, positions_shifts);

	//2014-4-3, we should use global chaining to tolerate error here
	fprintf(stdout,"%s\t",read->name);
	supported_positions_chainning(positions_shifts,read->len);
    #ifdef  OLD_CHAIN
	int count_debug=0;

	if(positions_shifts.size() > 1)
	{	count_debug ++;
		/*
		for(int i=0; i< positions_shifts.size(); ++i)
		{
			fprintf(stderr,"%c%llu ",positions_shifts[i].strand,positions_shifts[i].position-positions_shifts[i].shift);
		}
		fprintf(stderr,"\n");
		*/
		std::map<int64_t,int64_t> position_freq;
		for(int i=0; i<positions_shifts.size(); i++){
			position_freq[(positions_shifts[i].position - static_cast<int64_t>(positions_shifts[i].shift))*positions_shifts[i].strand] = 0;
		}
		for(int i=0; i<positions_shifts.size(); i++){
			position_freq[(positions_shifts[i].position - static_cast<int64_t>(positions_shifts[i].shift))*positions_shifts[i].strand] +=1;
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
				fprintf(stderr,"%s\t",read->name);
				for(int i=0; i<posFreqDic.size();++i)
				{
					fprintf(stderr,"%lld_%llu ",posFreqDic[i].first,posFreqDic[i].second);
				}
				fprintf(stderr,"\n");
			}
		}
	}
	else{
		fprintf(stderr,"%s\n",read->name);
	}
#endif

}

static int lea_map_core(ReferenceInfo reference_info , Tables kmer_position_tables,read_t *readsChunk, uint32_t n_seqs, Options opt) {
	uint32_t i = 0;
	read_t  *read;

	while ((uint32_t) i < (uint32_t) n_seqs) {
		read = &readsChunk[i++]; //A read.
		lea_map_single_read(reference_info ,kmer_position_tables,read,opt);
	}
	free_read_seq(n_seqs, readsChunk);
	return 0;
}

int lea_map(char *refFile, char *readsFile, Options opt) {
	int ret, i;
	Tables kmer_position_tables;
	ReferenceInfo reference_info;
	char *tableFn;
	uint32_t total = 0, n_seqs = 0;
	uint32_t readsBp = 0;
	seqio_t *seq;
	read_t *readsChunk, *p;
	char *prefix = refFile;
	clock_t t, mapStart;

	reference_info.bns = bns_restore(prefix);

	opt.chunck_size =  (0x2ffff);
	opt.index_parameter.kmer_len = get_kmer_length(reference_info.bns->l_pac);

	opt.kmer_table_len =  pow(2,(opt.index_parameter.kmer_len * 4))+10;
	opt.index_parameter.index_reverse_complement = false;
	opt.index_parameter.distinct_level = 1;
	opt.index_parameter.char_size = 2;
	opt.index_parameter.char_int = 0;


	fprintf(stderr, "[dmap_aln] start \n");
	t = clock();
	fprintf(stderr, "[dmap_aln] retrieving index info ...\n");


	reference_info.pac = static_cast<uint8_t *>(calloc((reference_info.bns->l_pac) / 4 + 1, 1));
	reference_info.l_pac = reference_info.bns->l_pac;

	if (reference_info.pac == 0) {
		fprintf(stderr, "[dmap_aln] insufficient memory!\n");
		exit(1);
	}

	fread(reference_info.pac, 1, (reference_info.bns->l_pac) / 4 + 1, reference_info.bns->fp_pac);
	fprintf(stderr, "[dmap_aln] number of characters of reference  %llu \n",
			reference_info.bns->l_pac);


	tableFn = (char*)calloc(strlen(refFile) + 10, 1);
	strcpy(tableFn, refFile);strcat(tableFn, ".tb");strcat(tableFn,"0");
	kmer_position_tables.kmer_position_table_AA = table_restore(tableFn,opt);


	tableFn = (char*)calloc(strlen(refFile) + 10, 1);
	strcpy(tableFn, refFile);strcat(tableFn, ".tb");strcat(tableFn,"1");
	kmer_position_tables.kmer_position_table_CC = table_restore(tableFn,opt);

	tableFn = (char*)calloc(strlen(refFile) + 10, 1);
	strcpy(tableFn, refFile);strcat(tableFn, ".tb");strcat(tableFn,"2");
	kmer_position_tables.kmer_position_table_GG = table_restore(tableFn,opt);

	tableFn = (char*)calloc(strlen(refFile) + 10, 1);
	strcpy(tableFn, refFile);strcat(tableFn, ".tb");strcat(tableFn,"3");
	kmer_position_tables.kmer_position_table_TT = table_restore(tableFn,opt);




	fprintf(stderr, "[dmap_aln] retrieve success! %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);

	//#define  LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
	#ifdef LEA_INDEX_REGLEN_BUILD_TABLE_DEBUG
	  uint32_t countDistinct =0,countOccur = 0;
	  uint32_t sum = 0;
	  for(uint64_t i = 0; i < static_cast<uint64_t>(opt.kmer_table_len-1); ++i){
		  if (kmer_position_table[i] != 0) {
			  countOccur +=1;
			  sum +=kmer_position_table[i];
			  fprintf(stderr, "%llu\t%llu\n",i,kmer_position_table[i]);
			  }
		  if (kmer_position_table[i]!=0 && kmer_position_table[i] != NOT_UNIQUE) countDistinct +=1;
	  }//for

	fprintf(stderr, "reference_info.l_pac %d\n",reference_info.l_pac);
	fprintf(stderr, "kmer %d\n",opt.index_parameter.kmer_len);
	fprintf(stderr, "total cell:%llu\n",opt.kmer_table_len-1);
	fprintf(stderr, "different occur %llu distinct:%llu   percentage: %f",countOccur, countDistinct,(float)countDistinct/countOccur);
	exit(1);
	#endif



	mapStart = clock();

	seq = seq_open(readsFile);
	/*Core loop, reads mapping. To save memory, sequentially process chunk of
	 * reads then write results to disk
	 * */
	readsChunk = (read_t*) calloc((opt.chunck_size) + 1, sizeof(read_t));

	while (kseq_read(seq->ks) >= 0) {
		total += 1;
		p = &readsChunk[n_seqs++];
		{
			//Read a sequence
			p->name = strdup((const char*) seq->ks->name.s);
			p->len = seq->ks->seq.l;
			readsBp += p->len;

			//copy read sequence
			p->seq = (ubyte_t*) calloc(p->len, 1);
			for (i = 0; i != p->len; i++)
				p->seq[i] = nst_nt4_table[(int) seq->ks->seq.s[i]];

			// copy read quality
			if (seq->ks->qual.l) {
				p->qual = (ubyte_t*) strdup((char*) seq->ks->qual.s);
				p->qual_len = seq->ks->qual.l;
			}
		}

		if (n_seqs == (opt.chunck_size) || readsBp > (opt.chunck_size) * 100) {

			fprintf(stderr, "[dmap_aln] read %d sequences ...", total);
			t = clock();
			ret = lea_map_core(reference_info , kmer_position_tables,readsChunk, n_seqs, opt);
			fprintf(stderr, "%.2f sec\n",
					(float) (clock() - t) / CLOCKS_PER_SEC);
			n_seqs = 0;
			readsBp = 0;
			readsChunk = (read_t*) calloc((opt.chunck_size) + 1,
					sizeof(read_t));
		}
	}
	fprintf(stderr, "[dmap_aln] read %d sequences ...\n", total);
	ret = lea_map_core(reference_info , kmer_position_tables,readsChunk, n_seqs, opt);

	fprintf(stderr, "[dmap_aln] map time %.2f sec\n",
				(float) (clock() - mapStart) / CLOCKS_PER_SEC);

	free(tableFn);
	free(kmer_position_tables.kmer_position_table_AA);
	free(kmer_position_tables.kmer_position_table_CC);
	free(kmer_position_tables.kmer_position_table_GG);
	free(kmer_position_tables.kmer_position_table_TT);

	free(seq);
	return 0;
}
