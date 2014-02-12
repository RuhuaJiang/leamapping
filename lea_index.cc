/*
 * lea_index.cc
 *
 *  Created on: Feb 7, 2014
 *      Author: jiang
 */

#include<stdio.h>
#include<string.h>
#include<unistd.h>



#include "lea_types.h"
#include "lea_index_reglen.h"

//Parsing and setting options
static void get_options(int argc, char *argv[],Options *opt)
{

	return;
}

//Returns success or not. Parsing and setting argument, decide which
//index algorithm to use.
int lea_index(int argc, char *argv[])
{
	int alg_type;
	char *indexFile;
	Options opt;

	//parsing arguments
	alg_type = 1;
	indexFile = argv[1];

	get_options(argc, argv,&opt);

	if(alg_type == 1){
		//Increase_and Decrease
		lea_index_region_length(indexFile,opt );

	}
	return 0;
}

