/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "datastruct.h"
#include "prepDB.h"
#include "filt.h"
#include "dbs.h"
#include "precalc.h"
#include "common.h"
#include "query.h"
#include "step.h"
#include "occ.h"
#include "selectmarker.h"
#include "candidate.h"
#include "tol.h"

#define DIFF (((marker->marker_fw_tail - marker->marker_fw_head) <= 100) ? tol[marker->marker_fw_tail - marker->marker_fw_head] : tol[100])
#define CompDIFF (((marker->marker_cp_tail - marker->marker_cp_head) <= 100) ? comp_tol[marker->marker_cp_tail - marker->marker_cp_head] : comp_tol[100])

steps_t* setStepDBS(dbs_t* dbs, filter_pars_t* pars)
{
	setBlockSNP(dbs, pars);
	setStepSNP(dbs, pars);
	int64_t steps_count = countSteps(dbs, pars);
	steps_t* steps = newSteps(dbs, steps_count);
	setSteps(dbs, steps, pars);
	return steps;
}

precalc_t* precalculation(char* dbs_filename, dbs_t* dbs, filter_pars_t* pars)
{
	char* precalc_filename = (char*)malloc(strlen(dbs_filename) + 6);
	sprintf(precalc_filename, "%s.prec", dbs_filename);	

	if(fileExists(precalc_filename))
	{
		printf("load %s file ...\n", precalc_filename);
		return loadPrecalc(dbs, precalc_filename);
	}
	else
	{
		printf("generate %s file ...\n", precalc_filename);
		return precalcDBS(dbs, precalc_filename, pars);
	}
}

void mergeSet(selected_marker_t* marker, precalc_t* precalc, occ_t* occ, dbs_t* dbs, bitmap_t* bitmap, int step_index, int* tol, int* comp_tol)
{
	int64_t i, j;

	for(i = marker->marker_fw_head; i < marker->marker_fw_tail; i++)
	{
		for(j = 0; j < precalc->set_size[marker->marker_fw[i]]; j++)
		{
			occ->occ_fw[precalc->set[marker->marker_fw[i]][j]]++;
		}
	}

	for(i = marker->marker_cp_head; i < marker->marker_cp_tail; i++)
        {
                for(j = 0; j < precalc->complement_set_size[marker->marker_cp[i]]; j++)
                {
                        occ->occ_cp[precalc->complement_set[marker->marker_cp[i]][j]]++;
                }
        }

	for(i = 0; i < dbs->tuple_count; i++)
	{
		if(occ->occ_fw[i] <= DIFF && occ->occ_cp[i] <= CompDIFF)
		{
			bitmap->res[i][step_index] = 1;
		}
	}
}

void replaceSet(selected_marker_t* marker, precalc_t* precalc, occ_t* occ, dbs_t* dbs, bitmap_t* bitmap, int step_index, int* tol, int* comp_tol)
{
	int64_t i, j;

	for(i = marker->marker_fw_head; i < marker->marker_fw_head + marker->marker_count_fw[marker->marker_count_fw_head]; i++)
	{
		for(j = 0; j < precalc->set_size[marker->marker_fw[i]]; j++)
                {
                        occ->occ_fw[precalc->set[marker->marker_fw[i]][j]]--;
                }
	}

	for(i = marker->marker_fw_tail - marker->marker_count_fw[marker->marker_count_fw_tail - 1]; i < marker->marker_fw_tail; i++)
	{
		for(j = 0; j < precalc->set_size[marker->marker_fw[i]]; j++)
                {
                        occ->occ_fw[precalc->set[marker->marker_fw[i]][j]]++;
                }
	}

	marker->marker_fw_head += marker->marker_count_fw[marker->marker_count_fw_head];
	marker->marker_count_fw_head++;

	for(i = marker->marker_cp_head; i < marker->marker_cp_head + marker->marker_count_cp[marker->marker_count_cp_head]; i++)
        {
                for(j = 0; j < precalc->complement_set_size[marker->marker_cp[i]]; j++)
                {
                        occ->occ_cp[precalc->complement_set[marker->marker_cp[i]][j]]--;
                }
        }

        for(i = marker->marker_cp_tail - marker->marker_count_cp[marker->marker_count_cp_tail - 1]; i < marker->marker_cp_tail; i++)
        {
                for(j = 0; j < precalc->complement_set_size[marker->marker_cp[i]]; j++)
                {
                        occ->occ_cp[precalc->complement_set[marker->marker_cp[i]][j]]++;
                }
        }

        marker->marker_cp_head += marker->marker_count_cp[marker->marker_count_cp_head];
        marker->marker_count_cp_head++;

	for(i = 0; i < dbs->tuple_count; i++)
        {
                if(occ->occ_fw[i] <= DIFF && occ->occ_cp[i] <= CompDIFF)
                {
			bitmap->res[i][step_index] = 1;
                }
        }
}

bitmap_t* newBitmap(steps_t* steps, dbs_t* dbs)
{
	int64_t i;
	bitmap_t* bitmap = (bitmap_t*) malloc(1 * sizeof(bitmap_t));
        if(bitmap == NULL) printf("malloc error: bitmap\n");
        if((bitmap->res = (int8_t**) malloc(dbs->tuple_count * sizeof(int8_t*))) == NULL) printf("malloc error: bitmap\n");
	for(i = 0; i < dbs->tuple_count; i++)
	{
		if((bitmap->res[i] = (int8_t*) malloc(steps->steps_count * sizeof(int8_t))) == NULL) printf("malloc error: bitmap\n");
	}
	return bitmap;
}

void clearBitmap(bitmap_t* bitmap, dbs_t* dbs, steps_t* steps)
{
	int64_t i;
	for(i = 0; i < dbs->tuple_count; i++)
	{
		memset(bitmap->res[i], 0, steps->steps_count * sizeof(int8_t));
	}
}

void freeBitmap(dbs_t* dbs, bitmap_t* bitmap)
{
	int64_t i;
	for(i = 0; i < dbs->tuple_count; i++)
	{
		free(bitmap->res[i]);
	}
	free(bitmap->res);
	free(bitmap);
}

void computeIBD(bitmap_t* bitmap, dbs_t* dbs, steps_t* steps, filter_pars_t* pars, queries_t* queries, int64_t query_index, candidates_t* candidates)
{
	int64_t i, j;
	double start = -1, end = -1;
	int64_t start_step = -1, end_step = -1;
	int flag;
	for(i = 0; i < dbs->tuple_count; i++)
	{
		flag = 0;
		for(j = 0; j < steps->steps_count; j++)
		{
			if(bitmap->res[i][j])
			{
				if(!bitmap->res[i][j - 1])
				{
					if(flag)
					{
						if(steps->start_genetic_position[j - pars->step_count + 1] - pars->step_size <= end)
						{
							end = steps->end_genetic_position[j] + pars->step_size;
							end_step = j + 1;
						}
						else
						{
							if(pars->pairwise == 1 && query_index == i)
							{
							}
							else
							{
								if(start < dbs->genetic_position[0]) start = dbs->genetic_position[0];
								if(end > dbs->genetic_position[dbs->SNP_count - 1]) end = dbs->genetic_position[dbs->SNP_count - 1];
								if(start_step < 0) start_step = 0;
								if(end_step >= steps->steps_count) end_step = steps->steps_count - 1;

								candidates->query_index[candidates->set_size] = query_index;
								candidates->tuple_index[candidates->set_size] = i;
								candidates->start_step[candidates->set_size] = start_step;
        	                                                candidates->end_step[candidates->set_size] = end_step;
								candidates->set_size++;
								if(candidates->set_size == candidates->max_set_size)
								{
									candidates->max_set_size *= 2;
									printf("Note: candidate list size is updated: %" PRIu64 "\n", candidates->max_set_size);
									reallocCandidates(candidates);
								}

								start = steps->start_genetic_position[j - pars->step_count + 1] - pars->step_size;
								end = steps->end_genetic_position[j] + pars->step_size;
								start_step = j - pars->step_count;
								end_step = j + 1;
							}
						}
					}
					else
					{
						start = steps->start_genetic_position[j - pars->step_count + 1] - pars->step_size;
						end = steps->end_genetic_position[j] + pars->step_size;
						start_step = j - pars->step_count;
						end_step = j + 1;
						flag = 1;
					}
				}
				else
				{
					end = steps->end_genetic_position[j] + pars->step_size;
					end_step = j + 1;
				}
			}
		}
		if(flag)
		{
			if(pars->pairwise == 1 && query_index == i)
			{
			}
			else
			{
				if(start < dbs->genetic_position[0]) start = dbs->genetic_position[0];
				if(end > dbs->genetic_position[dbs->SNP_count - 1]) end = dbs->genetic_position[dbs->SNP_count - 1];
				if(start_step < 0) start_step = 0;
                        	if(end_step >= steps->steps_count) end_step = steps->steps_count - 1;

				candidates->query_index[candidates->set_size] = query_index;
        	                candidates->tuple_index[candidates->set_size] = i;
				candidates->start_step[candidates->set_size] = start_step;
        	                candidates->end_step[candidates->set_size] = end_step;
                	        candidates->set_size++;
				if(candidates->set_size == candidates->max_set_size)
	                        {
        	                        candidates->max_set_size *= 2;
					printf("Note: candidate list size is updated: %" PRIu64 "\n", candidates->max_set_size);
                        	        reallocCandidates(candidates);
                        	}
			}
		}
	}
}

candidates_t* runQuery(dbs_t* dbs, queries_t* queries, precalc_t* precalc, steps_t* steps, filter_pars_t* pars, int* tol, int* comp_tol)
{
	selected_marker_t* marker = newSelectedMarker(dbs);
	occ_t* occ = newOcc(dbs);
	bitmap_t* bitmap = newBitmap(steps, dbs);
	candidates_t* candidates = newCandidates();

	int64_t i, j;
	for(i = 0; i < queries->query_count; i++)
	{
		clearSelectedMarker(marker, dbs);
		clearOcc(occ, dbs);
		clearBitmap(bitmap, dbs, steps);

		int64_t step;
		int64_t start_position, end_position;
		for(step = 0; step < steps->steps_count; step++)
		{
			start_position = steps->start_SNP[step];
			end_position = steps->end_SNP[step];

			int marker_count = 0;
			for(j = start_position; j <= end_position; j++)
			{
				if(dbs->minor_allele_frequency[j] < 0.2) continue; 
				if(precalc->cluster[j] == -1) continue;
				if(marker->cluster_fw[precalc->cluster[j]] == 1) continue;
				if(queries->query[i].genotype[j] == 0)	// query i has two minor alleles in marker j
				{
					marker_count++;
					marker->marker_fw[marker->marker_fw_tail] = j;
					marker->marker_fw_tail++;
					marker->cluster_fw[precalc->cluster[j]] = 1;
				}
			}
			marker->marker_count_fw[marker->marker_count_fw_tail] = marker_count;
			marker->marker_count_fw_tail++;

			marker_count = 0;
			for(j = start_position; j <= end_position; j++)
                        {
                                if(precalc->complement_cluster[j] == -1) continue;
                                if(marker->cluster_cp[precalc->complement_cluster[j]] == 1) continue;
                                if(queries->query[i].genotype[j] == 2)  // query i has two major alleles in marker j
                                {
                                        marker_count++;
                                        marker->marker_cp[marker->marker_cp_tail] = j;
                                        marker->marker_cp_tail++;
                                        marker->cluster_cp[precalc->complement_cluster[j]] = 1;
                                }
                        }
                        marker->marker_count_cp[marker->marker_count_cp_tail] = marker_count;
                        marker->marker_count_cp_tail++;

			if(step > pars->step_count - 1)
			{
				replaceSet(marker, precalc, occ, dbs, bitmap, step, tol, comp_tol);
			}
			else if(step == pars->step_count - 1)
			{
				mergeSet(marker, precalc, occ, dbs, bitmap, step, tol, comp_tol);
			}
		}
		computeIBD(bitmap, dbs, steps, pars, queries, i, candidates);
	}

	freeSelectedMarker(marker);
	freeOcc(occ);
	freeBitmap(dbs, bitmap);	

	return candidates;
}

int comp(candidates_t* candidates, int64_t i, int64_t j)
{
	if(candidates->query_index[i] > candidates->query_index[j]) return 1;
	else if(candidates->query_index[i] < candidates->query_index[j]) return -1;

	if(candidates->tuple_index[i] > candidates->tuple_index[j]) return 1;
	else if(candidates->tuple_index[i] < candidates->tuple_index[j]) return -1;

	if(candidates->start_step[i] > candidates->start_step[j]) return 1;
	else if(candidates->start_step[i] < candidates->start_step[j]) return -1;

	if(candidates->end_step[i] > candidates->end_step[j]) return 1;
        else if(candidates->end_step[i] < candidates->end_step[j]) return -1;

	return 0;
}

void swap(candidates_t* candidates, int64_t i, int64_t j)
{
	int64_t temp;
	temp = candidates->query_index[i];
	candidates->query_index[i] = candidates->query_index[j];
	candidates->query_index[j] = temp;

	temp = candidates->tuple_index[i];
	candidates->tuple_index[i] = candidates->tuple_index[j];
	candidates->tuple_index[j] = temp;

	temp = candidates->start_step[i];
	candidates->start_step[i] = candidates->start_step[j];
	candidates->start_step[j] = temp;

	temp = candidates->end_step[i];
	candidates->end_step[i] = candidates->end_step[j];
	candidates->end_step[j] = temp;
}

candidates_t* dedupIBD(candidates_t* candidates)
{
	candidates_t* dedup_candidates = newCandidates();
	if(candidates->set_size == 0) return dedup_candidates;

	int64_t i, j, temp;
        for(i = 0; i < candidates->set_size; i++)
	{
		if(candidates->query_index[i] > candidates->tuple_index[i])
		{
			temp = candidates->query_index[i];
			candidates->query_index[i] = candidates->tuple_index[i];
			candidates->tuple_index[i] = temp;
		}
	}

	for(i = 0; i < candidates->set_size - 1; i++)
	{
		for(j = i + 1; j < candidates->set_size; j++)
		{
			if(comp(candidates, i, j) > 0)
			{
				swap(candidates, i, j);
			}
		}
	}

	dedup_candidates->query_index[0] = candidates->query_index[0];
        dedup_candidates->tuple_index[0] = candidates->tuple_index[0];
        dedup_candidates->start_step[0] = candidates->start_step[0];
        dedup_candidates->end_step[0] = candidates->end_step[0];
        dedup_candidates->set_size++;

	for(i = 1; i < candidates->set_size; i++)
	{
		if(candidates->query_index[i] == dedup_candidates->query_index[dedup_candidates->set_size - 1] && candidates->tuple_index[i] == dedup_candidates->tuple_index[dedup_candidates->set_size - 1] && candidates->start_step[i] <= dedup_candidates->end_step[dedup_candidates->set_size - 1] + 1)
		{
			dedup_candidates->end_step[dedup_candidates->set_size - 1] = getMaxInt64(dedup_candidates->end_step[dedup_candidates->set_size - 1], candidates->end_step[i]);
		}
		else
		{
			dedup_candidates->query_index[dedup_candidates->set_size] = candidates->query_index[i];
		        dedup_candidates->tuple_index[dedup_candidates->set_size] = candidates->tuple_index[i];
		        dedup_candidates->start_step[dedup_candidates->set_size] = candidates->start_step[i];
		        dedup_candidates->end_step[dedup_candidates->set_size] = candidates->end_step[i];
		        dedup_candidates->set_size++;

			if(dedup_candidates->set_size == dedup_candidates->max_set_size)
                        {
                                dedup_candidates->max_set_size *= 2;
                                printf("Note: candidate list size is updated: %" PRIu64 "\n", dedup_candidates->max_set_size);
                                reallocCandidates(dedup_candidates);
                        }
		}
	}

	return dedup_candidates;
}

void printDummyOutput(char* output_filename, dbs_t* dbs, queries_t* queries)
{
	FILE* output = fileOpenW(output_filename);
        int64_t i, j;
        fprintf (output, "Query\tTuple\tStart\tSNP_index\tEnd\tSNP_index\n");
	for(i = 0; i < queries->query_count; i++)
        {
		for(j = 0; j < dbs->tuple_count; j++)
		{
	                fprintf (output, "%" PRIu64 "\t%" PRIu64 "\t%f\t%" PRIu64 "\t%f\t%" PRIu64 "\n", i, j, dbs->genetic_position[0], 0, dbs->genetic_position[dbs->SNP_count - 1], dbs->SNP_count - 1);
		}
        }
        fileClose(output);
}

int filterDB(char* dbs_filename, char* query_filename, char* output_filename, filter_pars_t* pars)
{
	clock_t t;
        t = clock();
	printf("Loading DBS ...\n");
        dbs_t* dbs = loadDBS(dbs_filename);
	steps_t* steps = setStepDBS(dbs, pars);
	
	printf("Loading queries ...\n");
        queries_t* queries = loadQuery(query_filename, dbs, pars);
	printf("Loaded %d queries\n", queries->query_count);
        printf("Total loading time: %.2f sec\n\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	if(pars->dummy == 1)
	{
		printDummyOutput(output_filename, dbs, queries);
		return 1;
	}

	t = clock();
        printf("Precalculating ...\n");
        precalc_t* precalc = precalculation(dbs_filename, dbs, pars);
        printf("Total precalculation time: %.2f sec\n\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	int* tol = NULL;
	int* comp_tol = NULL;
	if(pars->p_filename == NULL)
	{
		tol = selectTol(pars);
	}
	else
	{
		tol = loadTol(pars);
	}
	if(pars->q_filename == NULL)
	{
		comp_tol = selectCompTol(pars);
	}
	else
	{
		comp_tol = loadCompTol(pars);
	}
	if(tol == NULL || comp_tol == NULL) printf("error: bad p threshold\n");

	t = clock();
        candidates_t* candidates = runQuery(dbs, queries, precalc, steps, pars, tol, comp_tol);
        printf("Total query time: %.2f sec\n\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	if(pars->pairwise == 1)
	{
		t = clock();
		printf("Merging ...\n");
		candidates_t* dedup_candidates = dedupIBD(candidates);
		printf("Total merging time: %.2f sec\n\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		storeCandidates(output_filename, dbs, steps, dedup_candidates);
		freeCandidates(dedup_candidates);
	}
	else
	{
		storeCandidates(output_filename, dbs, steps, candidates);
		evalCandidates(candidates, dbs, steps, queries, pars);
	}

	freeSteps(steps);
	freeQuery(queries);
	freePrecalc(dbs, precalc);
	freeCandidates(candidates);
	freeDBS(dbs);

	return 0;
}
