/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "selectmarker.h"

selected_marker_t* newSelectedMarker(dbs_t* dbs)
{
        selected_marker_t* marker = (selected_marker_t*) malloc(1 * sizeof(selected_marker_t));
        if(marker == NULL) printf("malloc error: marker\n");
        if((marker->marker_fw = (int64_t*) malloc(dbs->SNP_count * sizeof(int64_t))) == NULL) printf("malloc error: marker\n");
        if((marker->marker_count_fw = (int64_t*) malloc(dbs->SNP_count * sizeof(int64_t))) == NULL) printf("malloc error: marker\n");
        if((marker->cluster_fw = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: marker\n");
        if((marker->marker_cp = (int64_t*) malloc(dbs->SNP_count * sizeof(int64_t))) == NULL) printf("malloc error: marker\n");
        if((marker->marker_count_cp = (int64_t*) malloc(dbs->SNP_count * sizeof(int64_t))) == NULL) printf("malloc error: marker\n");
        if((marker->cluster_cp = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: marker\n");
        return marker;
}

void clearSelectedMarker(selected_marker_t* marker, dbs_t* dbs)
{
        marker->marker_fw_head = 0;
        marker->marker_fw_tail = 0;
        marker->marker_count_fw_head = 0;
        marker->marker_count_fw_tail = 0;
        memset(marker->cluster_fw, 0, dbs->SNP_count * sizeof(int));

        marker->marker_cp_head = 0;
        marker->marker_cp_tail = 0;
        marker->marker_count_cp_head = 0;
        marker->marker_count_cp_tail = 0;
        memset(marker->cluster_cp, 0, dbs->SNP_count * sizeof(int));
}

void freeSelectedMarker(selected_marker_t* marker)
{
        free(marker->marker_fw);
        free(marker->marker_count_fw);
        free(marker->cluster_fw);

        free(marker->marker_cp);
        free(marker->marker_count_cp);
        free(marker->cluster_cp);
}

void printSelectedMarker(selected_marker_t* marker)
{
        int64_t i;
        for(i = marker->marker_fw_head; i < marker->marker_fw_tail; i++)
        {
                printf("%" PRIu64 ", ", marker->marker_fw[i]);
        }
        printf("\n");
        for(i = marker->marker_count_fw_head; i < marker->marker_count_fw_tail; i++)
        {
                printf("%" PRIu64 ", ", marker->marker_count_fw[i]);
        }
        printf("\n\n");

        for(i = marker->marker_cp_head; i < marker->marker_cp_tail; i++)
        {
                printf("%" PRIu64 ", ", marker->marker_cp[i]);
        }
        printf("\n");
        for(i = marker->marker_count_cp_head; i < marker->marker_count_cp_tail; i++)
        {
                printf("%" PRIu64 ", ", marker->marker_count_cp[i]);
        }
        printf("\n\n");
}

void printSelectedMarkerCount(selected_marker_t* marker)
{
	int64_t i;
	int64_t count_fw = 0, count_cp = 0;
	for(i = marker->marker_count_fw_head; i < marker->marker_count_fw_tail; i++)
        {
		count_fw += marker->marker_count_fw[i];
        }
	for(i = marker->marker_count_cp_head; i < marker->marker_count_cp_tail; i++)
        {
		count_cp += marker->marker_count_cp[i];
        }
	printf("%" PRIu64 "\t%" PRIu64 "\n", count_fw, count_cp);
}
