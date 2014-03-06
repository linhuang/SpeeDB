/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _SELECT_MARKER_H_
#define _SELECT_MARKER_H_

#include "datastruct.h"

selected_marker_t* newSelectedMarker(dbs_t* dbs);
void clearSelectedMarker(selected_marker_t* marker, dbs_t* dbs);
void freeSelectedMarker(selected_marker_t* marker);
void printSelectedMarker(selected_marker_t* marker);
void printSelectedMarkerCount(selected_marker_t* marker);

#endif
