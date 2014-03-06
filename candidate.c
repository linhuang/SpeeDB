/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "datastruct.h"
#include "candidate.h"

#define candidatesStart steps->start_genetic_position[candidates->start_step[i]]
#define candidatesEnd getMinDouble(steps->end_genetic_position[candidates->end_step[i]], dbs->genetic_position[dbs->SNP_count - 1])

candidates_t* newCandidates()
{
        candidates_t* candidates = (candidates_t*) malloc(1 * sizeof(candidates_t));
        if(candidates == NULL) printf("malloc error: candidates\n");
        candidates->set_size = 0;
        candidates->max_set_size = MAX_CANDIDATE_SIZE;
        if((candidates->query_index = (int64_t*) malloc(MAX_CANDIDATE_SIZE * sizeof(int64_t))) == NULL) printf("malloc error: candidates\n");
        if((candidates->tuple_index = (int64_t*) malloc(MAX_CANDIDATE_SIZE * sizeof(int64_t))) == NULL) printf("malloc error: candidates\n");
	if((candidates->start_step = (int64_t*) malloc(MAX_CANDIDATE_SIZE * sizeof(int64_t))) == NULL) printf("malloc error: candidates\n");
        if((candidates->end_step = (int64_t*) malloc(MAX_CANDIDATE_SIZE * sizeof(int64_t))) == NULL) printf("malloc error: candidates\n");
        return candidates;
}

void reallocCandidates(candidates_t* candidates)
{
	if((candidates->query_index = (int64_t*) realloc(candidates->query_index, candidates->max_set_size * sizeof(int64_t))) == NULL) printf("realloc error: candidates\n");
        if((candidates->tuple_index = (int64_t*) realloc(candidates->tuple_index, candidates->max_set_size * sizeof(int64_t))) == NULL) printf("realloc error: candidates\n");
	if((candidates->start_step = (int64_t*) realloc(candidates->start_step, candidates->max_set_size * sizeof(int64_t))) == NULL) printf("realloc error: candidates\n");
        if((candidates->end_step = (int64_t*) realloc(candidates->end_step, candidates->max_set_size * sizeof(int64_t))) == NULL) printf("realloc error: candidates\n");
}

void storeCandidates(char* output_filename, dbs_t* dbs, steps_t* steps, candidates_t* candidates)
{
        FILE* output = fileOpenW(output_filename);
        int64_t i;
        fprintf (output, "Query\tTuple\tStart\tSNP_index\tEnd\tSNP_index\n");
        for(i = 0; i < candidates->set_size; i++)
        {
                fprintf (output, "%" PRIu64 "\t%" PRIu64 "\t%f\t%" PRIu64 "\t%f\t%" PRIu64 "\n", candidates->query_index[i], candidates->tuple_index[i], candidatesStart, steps->start_SNP[candidates->start_step[i]], candidatesEnd, steps->end_SNP[candidates->end_step[i]]);
        }
        fileClose(output);
}

void freeCandidates(candidates_t* candidates)
{
	free(candidates->query_index);
	free(candidates->tuple_index);
	free(candidates->start_step);
	free(candidates->end_step);
	free(candidates);
}

void evalCandidates(candidates_t* candidates, dbs_t* dbs, steps_t* steps, queries_t* queries, filter_pars_t* pars)
{
        int64_t i;
        double total_candidate_length = 0.0, pass_ratio;
        for(i = 0; i < candidates->set_size; i++)
        {
                total_candidate_length += (candidatesEnd - candidatesStart);
        }
	double total_DB_length = (dbs->genetic_position[dbs->SNP_count - 1] - dbs->genetic_position[0]) * queries->query_count * dbs->tuple_count;
        pass_ratio = total_candidate_length / total_DB_length;
        printf("pass ratio = %f\n", pass_ratio);

	double total_IBD_length = 0.0;
	for(i = 0; i < queries->query_count; i++)
	{
		total_IBD_length += (queries->ibd_end[i] - queries->ibd_start[i]);
	}

	double total_overlap_length = 0.0;
	for(i = 0; i < candidates->set_size; i++)
        {
		if(candidates->tuple_index[i] == queries->db_index[candidates->query_index[i]])
		{
			double c_start = candidatesStart;
			double q_start = queries->ibd_start[candidates->query_index[i]];
			double c_end = candidatesEnd;
			double q_end = queries->ibd_end[candidates->query_index[i]];

			double o_start = getMaxDouble(c_start, q_start);
			double o_end = getMinDouble(c_end, q_end);
			if(o_end > o_start) total_overlap_length += (o_end - o_start);
		}
	}

	double partial_credit_sensitivity = total_overlap_length / total_IBD_length;
	printf("partial credit sensitivity = %f\n", partial_credit_sensitivity);

	double partial_credit_specificity = (total_DB_length - (total_candidate_length + total_IBD_length - total_overlap_length)) / (total_DB_length - total_IBD_length);
	printf("specificity = %f\n", partial_credit_specificity);

        int64_t sensitivity = 0;
        int found = 0;
        for(i = 0; i < candidates->set_size; i++)
        {
                if(!i || candidates->query_index[i] != candidates->query_index[i - 1])
                {
                        if(found) 
			{
				sensitivity++;
			}
			else
			{
			}
                        found = 0;
                        if(candidates->tuple_index[i] == queries->db_index[candidates->query_index[i]] && candidatesStart <= queries->ibd_start[candidates->query_index[i]] && candidatesEnd >= queries->ibd_end[candidates->query_index[i]]) 
			{
				found = 1;
			}
			else
			{
			}
                }
                else
                {
                        if(!found && candidates->tuple_index[i] == queries->db_index[candidates->query_index[i]] && candidatesStart <= queries->ibd_start[candidates->query_index[i]] && candidatesEnd >= queries->ibd_end[candidates->query_index[i]]) 
			{
				found = 1;
			}
			else
                        {
                        }
                }
        }
        if(found) 
	{
		sensitivity++;
	}
	else
        {
        }
        printf("sensitivity = %" PRIu64 "\n", sensitivity);
}
