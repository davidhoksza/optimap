/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef INDEXING_H
#define INDEXING_H

#include "types.h"

struct IndexRecord {
	int start_position;
	int end_position;
	int length;

	IndexRecord() : start_position(-1), end_position(-1), length(0) {};
	IndexRecord(int sp, int ep, int l) : start_position(sp), end_position(sp), length(l) {};
};

void analyze(std::vector<Fragment> optMap, std::vector<RMRead> refMap);
std::map<int, std::vector<IndexRecord> > init_index(RefMaps refMaps);
std::vector<IndexRecord> index_get_candidates(std::map<int, std::vector<IndexRecord> > &index, const Fragment &optMapFragment, const int threshold = INDEX_NEIGHBORHOOD_THRESHOLD);

#endif // INDEXING_H