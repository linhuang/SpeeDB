/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _CANDIDATE_H_
#define _CANDIDATE_H_

#define MAX_CANDIDATE_SIZE 10000000

candidates_t* newCandidates();
void reallocCandidates(candidates_t* candidates);
void freeCandidates(candidates_t* candidates);
void storeCandidates(char* output_filename, dbs_t* dbs, steps_t* steps, candidates_t* candidates);
void evalCandidates(candidates_t* candidates, dbs_t* dbs, steps_t* steps, queries_t* queries, filter_pars_t* pars);

#endif
