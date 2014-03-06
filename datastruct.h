/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _DATA_STRUCTURE_H_
#define _DATA_STRUCTURE_H_

#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

typedef struct
{
        int64_t SNP_count;
        int64_t haplotype_count;
} datasize_t;

typedef struct
{
        int *genotype;
} tuple_t;

typedef struct
{
        char *allele;
} haplotype_t;

typedef struct
{
        int64_t SNP_count;
        int64_t haplotype_count;
        int64_t tuple_count;
        double *genetic_position;
        int *block_SNP;
        int *step_SNP;
        char *minor_allele;
	char *major_allele;
        double *minor_allele_frequency;
        tuple_t *tuple;
	haplotype_t *haplotype;
} dbs_t;

typedef struct
{
        int query_count;
        tuple_t *query;
        int64_t *db_index;
        double *ibd_start;
        double *ibd_end;
        double *acc_start;
        double *acc_end;
} queries_t;

typedef struct
{
        int *set_size;  // query g = 2
        int **set;      // query g = 2
        int *cluster;   // query g = 2
        int *complement_set_size;       // query g = 0
        int **complement_set;           // query g = 0
        int *complement_cluster;        // query g = 0
} precalc_t;

typedef struct
{
	double block_size;
	double step_size;
	double eta_J;
	int step_count;
	int p_th_fw;
	int p_th_cp;
	char* p_filename;
	char* q_filename;
	int pairwise;
} filter_pars_t;

typedef struct
{
	double block_size;
	double segment_size;
	double theta;
	int query_count;
	int ibd;
} query_pars_t;

typedef struct
{
        int64_t segments_size;
        int64_t* segments;
} segments_t;

typedef struct
{
        int64_t steps_count;
        double* start_genetic_position;
        double* end_genetic_position;
	int64_t* start_SNP;
	int64_t* end_SNP;
} steps_t;

typedef struct
{
        int* occ_fw;
        int* occ_cp;
} occ_t;

typedef struct
{
        int64_t* marker_fw;
        int64_t* marker_count_fw;
        int64_t marker_fw_head, marker_fw_tail;
        int64_t marker_count_fw_head, marker_count_fw_tail;
        int* cluster_fw;

        int64_t* marker_cp;
        int64_t* marker_count_cp;
        int64_t marker_cp_head, marker_cp_tail;
        int64_t marker_count_cp_head, marker_count_cp_tail;
        int* cluster_cp;
} selected_marker_t;

typedef struct
{
        int8_t** res;
} bitmap_t;

typedef struct
{
	int64_t set_size;
	int64_t max_set_size;
	int64_t* query_index;
	int64_t* tuple_index;
	int64_t* start_step;
	int64_t* end_step;
} candidates_t;

#endif
