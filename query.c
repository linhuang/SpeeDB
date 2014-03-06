/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "datastruct.h"
#include "common.h"
#include "query.h"
#include "dbs.h"

segments_t* genSegments(dbs_t* dbs, query_pars_t* pars)
{
	segments_t* seg = (segments_t*) calloc(1, sizeof(segments_t));
	if(seg == NULL) printf("malloc error: seg\n");
	if((seg->segments = (int64_t*) malloc(dbs->SNP_count * sizeof(int64_t))) == NULL) printf("malloc error: seg\n");
	int64_t p = 0;
	double start = dbs->genetic_position[0], end;
        while(p < dbs->SNP_count)
        {
                end = start + pars->segment_size;
                if(dbs->genetic_position[p] > end)
                {
                        seg->segments[seg->segments_size] = p - 1;
                        seg->segments_size++;
                        start = end;
                }
                p++;
        }
        seg->segments[seg->segments_size] = p - 1;
        seg->segments_size++;
	return seg;
}

void freeSegments(segments_t* seg)
{
	free(seg->segments);
	free(seg);
}

void printSegments(segments_t* seg, dbs_t* dbs)
{
	int64_t i;
	for(i = 0; i < seg->segments_size; i++)
	{
		printf("%" PRIu64 "\t%f\n", seg->segments[i], dbs->genetic_position[seg->segments[i]]);
	}
}

queries_t* newQuery(dbs_t* dbs, int query_count)
{
	queries_t* queries = (queries_t*) calloc(1, sizeof(queries_t));
	if(queries == NULL) printf("malloc error: queries\n");
	queries->query_count = query_count;
	int64_t i;
	if((queries->query = (tuple_t*) malloc(queries->query_count * sizeof(tuple_t))) == NULL) printf("malloc error: queries\n");
        for(i = 0; i < queries->query_count; i++)
        {
                if((queries->query[i].genotype = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: queries\n");
        }
        if((queries->db_index = (int64_t*) malloc(queries->query_count * sizeof(int64_t))) == NULL) printf("malloc error: queries\n");
        if((queries->ibd_start = (double*) calloc(queries->query_count, sizeof(double))) == NULL) printf("malloc error: queries\n");
        if((queries->ibd_end = (double*) calloc(queries->query_count, sizeof(double))) == NULL) printf("malloc error: queries\n");
        if((queries->acc_start = (double*) calloc(queries->query_count, sizeof(double))) == NULL) printf("malloc error: queries\n");
        if((queries->acc_end = (double*) calloc(queries->query_count, sizeof(double))) == NULL) printf("malloc error: queries\n");
	return queries;
}

void freeQuery(queries_t* queries)
{
	int64_t i;
	free(queries->db_index);
	free(queries->ibd_start);
	free(queries->ibd_end);
	free(queries->acc_start);
	free(queries->acc_end);
	for(i = 0; i < queries->query_count; i++)
	{
		free(queries->query[i].genotype);
	}
	free(queries->query);
	free(queries);
}

void printQuery(dbs_t* dbs, queries_t* queries)
{
        int64_t i, j;
        for(i = 0; i < queries->query_count; i++)
        {
		printf("%" PRIu64 "\t%f\t%f\t%f\t%f\t", queries->db_index[i], queries->ibd_start[i], queries->ibd_end[i], queries->acc_start[i], queries->acc_end[i]);
                for(j = 0; j < dbs->SNP_count; j++)
                {
			printf("%d", queries->query[i].genotype[j]);
                }
		printf("\n");
        }
}

void genQUERY(char* input_filename, char* output_filename)
{
	dbs_t* dbs = loadDBS(input_filename);
	queries_t* queries = newQuery(dbs, dbs->tuple_count);

	int64_t i, j;
	for(i = 0; i < dbs->SNP_count; i++)
	{
                for(j = 0; j < dbs->tuple_count; j++)
		{
			queries->query[j].genotype[i] = dbs->tuple[j].genotype[i];
		}
	}

	storeQuery(output_filename, dbs, queries);

	char* meta_filename = (char*)malloc(strlen(output_filename) + 6);
        sprintf(meta_filename, "%s.meta", output_filename);
        FILE* meta = fileOpenWB(meta_filename);

	fileClose(meta);
	freeQuery(queries);
	freeDBS(dbs);
}

void genQ(char* input_filename, char* output_filename, query_pars_t* pars)
{
	dbs_t* dbs = loadDBS(input_filename);
	segments_t* seg = genSegments(dbs, pars);
	queries_t* queries = newQuery(dbs, pars->query_count);

	int64_t i, j, k;
        int64_t** segment_origin = (int64_t**) malloc(pars->query_count * 2 * sizeof(int64_t*));
	if(segment_origin == NULL) printf("malloc error: segment_origin\n");
        for(i = 0; i < pars->query_count * 2; i++)
        {
                if((segment_origin[i] = (int64_t*) malloc(dbs->SNP_count * sizeof(int64_t))) == NULL) printf("malloc error: segment_origin\n");
        }
        int64_t origin;
        for(i = 0; i < pars->query_count * 2; i++)
        {
                for(j = 0; j < seg->segments_size; j++)
                {
                        origin = rand() % dbs->haplotype_count;
                        if(j == 0)
                        {
                                for(k = 0; k <= seg->segments[j]; k++) segment_origin[i][k] = origin;
                        }
                        else
                        {
                                for(k = seg->segments[j - 1] + 1; k <= seg->segments[j]; k++) segment_origin[i][k] = origin;
                        }
                }
        }

	char* meta_filename = (char*)malloc(strlen(output_filename) + 6);
	sprintf(meta_filename, "%s.meta", output_filename);
	FILE* meta = fileOpenWB(meta_filename);

	if(pars->ibd)
	{
		double start, end;
		double genetic_start, genetic_end;
		for(i = 0; i < pars->query_count * 2; i += 2)
		{
			origin = rand() % dbs->haplotype_count;
			start = (rand() % (int)(dbs->genetic_position[dbs->SNP_count - 1] - dbs->genetic_position[0] - pars->block_size)) + dbs->genetic_position[0];
			end = start + pars->block_size;
			genetic_start = -1, genetic_end = -1;
			for(j = 0; j < dbs->SNP_count; j++)
			{
				if(dbs->genetic_position[j] >= start && dbs->genetic_position[j] <= end)
				{
					segment_origin[i][j] = origin;
					if(genetic_start < 0) genetic_start = dbs->genetic_position[j];
					genetic_end = dbs->genetic_position[j];
				}
			}
			
			int64_t value = origin / 2;
			fwrite(&value, sizeof(int64_t), 1, meta);
			fwrite(&start, sizeof(double), 1, meta);
			fwrite(&end, sizeof(double), 1, meta);
			fwrite(&genetic_start, sizeof(double), 1, meta);
                        fwrite(&genetic_end, sizeof(double), 1, meta);
		}
	}

	fileClose(meta);
	free(meta_filename);

	int genotype;
	int64_t k1, k2;
	char c1, c2;
	for(i = 0; i < dbs->SNP_count; i++)
        {
                for(j = 0; j < pars->query_count; j++)
                {
                        genotype = 0;
                        k1 = segment_origin[j * 2][i];
			c1 = dbs->haplotype[k1].allele[i];
			k2 = segment_origin[j * 2 + 1][i];
                        c2 = dbs->haplotype[k2].allele[i];
			if(isMissing(c1) || isMissing(c2))
			{
				queries->query[j].genotype[i] = 1;
				continue;
			}

			if(c1 != dbs->minor_allele[i])
                        {
                                if((double)(rand()) / (double)(RAND_MAX) < pars->theta)
                                {
                                }
                                else
                                {
                                        genotype++;
                                }
                        }
                        else
                        {
                                if((double)(rand()) / (double)(RAND_MAX) < pars->theta)
                                {
                                        genotype++;
                                }
                                else
                                {
                                }
                        }
			if(c2 != dbs->minor_allele[i])
                        {
                                if((double)(rand()) / (double)(RAND_MAX) < pars->theta)
                                {
                                }
                                else
                                {
                                        genotype++;
                                }
                        }
                        else
                        {
                                if((double)(rand()) / (double)(RAND_MAX) < pars->theta)
                                {
                                        genotype++;
                                }
                                else
                                {
                                }
                        }
			queries->query[j].genotype[i] = genotype;
                }
        }
	storeQuery(output_filename, dbs, queries);

	freeSegments(seg);
	for(i = 0; i < pars->query_count * 2; i++)
	{
		free(segment_origin[i]);
	}
	free(segment_origin);
	freeQuery(queries);
	freeDBS(dbs);
}

void storeQuery(char* output_filename, dbs_t* dbs, queries_t* queries)
{
	int64_t i, j;
	FILE* output = fileOpenWB(output_filename);
	fwrite(&queries->query_count, sizeof(int), 1, output);
	for(i = 0; i < dbs->SNP_count; i++)
        {
                for(j = 0; j < queries->query_count; j++)
                {
			fwrite(&queries->query[j].genotype[i], sizeof(int), 1, output);
		}
	}
	fileClose(output);
}

queries_t* loadQuery(char* query_filename, dbs_t* dbs, filter_pars_t* pars)
{
	int64_t i, j;
        FILE* input = fileOpenRB(query_filename);
	int ret, query_count;
        ret = fread(&query_count, sizeof(int), 1, input);
	queries_t* queries = newQuery(dbs, query_count);
	queries->query_count = query_count;
        for(i = 0; i < dbs->SNP_count; i++)
        {
                for(j = 0; j < queries->query_count; j++)
                {
                        ret = fread(&queries->query[j].genotype[i], sizeof(int), 1, input);
			if(ret != 1)
			{
				printf("Loading error\n");
				exit(1);
			}
                }
        }
        fileClose(input);

	char* meta_filename = (char*)malloc(strlen(query_filename) + 6);
        sprintf(meta_filename, "%s.meta", query_filename);
        FILE* meta = fileOpenRB(meta_filename);
	for(i = 0; i < queries->query_count; i++)
	{
		ret = fread(&queries->db_index[i], sizeof(int64_t), 1, meta);
		ret = fread(&queries->ibd_start[i], sizeof(double), 1, meta);
		ret = fread(&queries->ibd_end[i], sizeof(double), 1, meta);
		ret = fread(&queries->acc_start[i], sizeof(double), 1, meta);
                ret = fread(&queries->acc_end[i], sizeof(double), 1, meta);
	}
	fileClose(meta);
	free(meta_filename);
	return queries;
}
