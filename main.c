/*
 *
 * Program SpeeDB: an effective filter for IBD detection in large data sets
 * by Lin Huang <linhuang@cs.stanford.edu> 2014
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Lin Huang
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "datastruct.h"
#include "prepDB.h"
#include "filt.h"
#include "query.h"
#include "common.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.0"
#endif

static int usage() 
{
	printf("Program: speedb\n");
	printf("Version: %s\n", PACKAGE_VERSION);
	printf("Contact: Lin Huang <linhuang@cs.stanford.edu>\n\n");
	printf("Usage:   speedb command [options] \n");
	printf("Command: tped2dbs    	produce DBS file\n");
	printf("         dbs2query      produce QUERY file\n");
	printf("         genQ       	generate queries\n");
	printf("         filter     	find IBD candidate segments\n");
	printf("\n");
	return 1;
}

static int tped2db_usage() 
{
	printf("Usage:   speedb tped2dbs <input.tped> <output.dbs> \n");
	printf("\n");
	return 1;
}

static int dbs2query_usage()
{
        printf("Usage:   speedb dbs2query <input.dbs> <output.query> \n");
        printf("\n");
        return 1;
}

static int filter_usage() 
{
        printf("Usage:   speedb filter [options] <input.dbs> <input.query> <output.candidate> \n");
	printf("Options: b    block size [default: 4]\n");
	printf("         s    step size [default: 0.5]\n");
	printf("         e    etaJ [default: 0.9]\n");
	printf("         1    filter 1's p threshold: 1e-x [options: 0, 1, 2; default: 1]\n");
	printf("         2    filter 2's p threshold: 1e-x [options: 0, 1, 2; default: 1]\n");
	printf("         p    read filter 1's p threshold from file\n");
	printf("         q    read filter 2's p threshold from file\n");
	printf("         w    pairwise filtering [default: queries are not from the database]\n");
        printf("\n");
        return 1;
}

static int genQ_usage() 
{
        printf("Usage:   speedb genQ [options] <input.dbs> <output.query> \n");
	printf("Options: b    block size [default: 4]\n");
	printf("         g    segment size [default: 0.2]\n");
	printf("         t    theta [default: 0.005]\n");
	printf("         c    query count [default: 1000]\n");
	printf("         N    no artificial IBD segments [default: having IBD segments]\n");
        printf("\n");
        return 1;
}

void set_default_filter_pars(filter_pars_t* pars)
{
	pars->block_size = 4;
	pars->step_size = 0.5;
	pars->eta_J = 0.9;
	pars->p_th_fw = 1;
	pars->p_th_cp = 1;
	pars->p_filename = NULL;
	pars->q_filename = NULL;
	pars->pairwise = 0;
}

void set_default_query_pars(query_pars_t* pars)
{
	pars->block_size = 4;
	pars->segment_size = 0.2;
	pars->theta = 0.005;
	pars->query_count = 1000;
	pars->ibd = 1;
}

void printFilterPars(filter_pars_t* pars)
{
	printf("block size: %f\tstep size: %f\tetaJ: %f\tstep count: %d\tfilter 1's threshold: 1e-%d\tfilter 2's threshold: 1e-%d\n", pars->block_size, pars->step_size, pars->eta_J, pars->step_count, pars->p_th_fw, pars->p_th_cp);
	if(pars->pairwise) printf("Pairwise inferrence mode\n");
}

void printQueryPars(query_pars_t* pars)
{
	if(pars->ibd)
	{
		printf("[IBD]\t");
	}
	else
	{
		printf("[non-IBD]\t");
	}
	printf("block size: %f\tsegment size: %f\ttheta: %f\tquery count: %d\n", pars->block_size, pars->segment_size, pars->theta, pars->query_count);
}

int main(int argc, char *argv[]) 
{
	if(argc < 2) return usage();
	if(strcmp(argv[1], "tped2dbs") == 0) 
	{
		if(argc < 4)
		{
			tped2db_usage();
			exit(1);
		}
		genDB(argv[2], argv[3]);
	}
	else if(strcmp(argv[1], "dbs2query") == 0)
        {
                if(argc < 4)
                {
                        dbs2query_usage();
                        exit(1);
                }
                genQUERY(argv[2], argv[3]);
        }
	else if(strcmp(argv[1], "filter") == 0) 
	{
		if(argc < 5)
                {
                        filter_usage();
                        exit(1);
                }
		filter_pars_t* pars = (filter_pars_t*) calloc(1, sizeof(filter_pars_t));
		set_default_filter_pars(pars);

		int c;
		while ((c = getopt(argc-1, argv+1, "b:s:e:1:2:p:q:w")) >= 0) {
				switch (c) {
					case 'b': pars->block_size = atof(optarg); break;
					case 's': pars->step_size = atof(optarg); break;
					case 'e': pars->eta_J = atof(optarg); break;
					case '1': pars->p_th_fw = atoi(optarg); break;
					case '2': pars->p_th_cp = atoi(optarg); break;
					case 'p': pars->p_filename = (char*) malloc((strlen(optarg) + 1) * sizeof(char)); strcpy(pars->p_filename, optarg); break;
					case 'q': pars->q_filename = (char*) malloc((strlen(optarg) + 1) * sizeof(char)); strcpy(pars->q_filename, optarg); break;
					case 'w': pars->pairwise = 1; break;
					case '?': filter_usage(); return 1;
					default: return 1;
				}
		}
		pars->step_count = ceil(pars->block_size / pars->step_size) - 1;
		printFilterPars(pars);

		filterDB(argv[optind+1], argv[optind+2], argv[optind+3], pars);
		free(pars);
	}
	else if(strcmp(argv[1], "genQ") == 0)
        {
                if(argc < 4)
                {
                        genQ_usage();
                        exit(1);
                }
                query_pars_t* pars = (query_pars_t*)calloc(1, sizeof(query_pars_t));
		set_default_query_pars(pars);

		int c;
                while ((c = getopt(argc-1, argv+1, "b:g:t:c:N")) >= 0) {
                                switch (c) {
                                        case 'b': pars->block_size = atof(optarg); break;
					case 'g': pars->segment_size = atof(optarg); break;
                                        case 't': pars->theta = atof(optarg); break;
					case 'c': pars->query_count = atoi(optarg); break;
					case 'N': pars->ibd = 0; break;
                                        case '?': filter_usage(); return 1;
                                        default: return 1;
                                }
                }
                printQueryPars(pars);

		genQ(argv[optind+1], argv[optind+2], pars);
		free(pars);
        }
	else
	{
		printf("Error: Unknown command '%s'\n", argv[1]);
		usage();
	}
	return 0;
}
