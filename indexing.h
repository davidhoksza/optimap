/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef INDEXING_H
#define INDEXING_H

#include "types.h"
//#include "gzstream/gzstream.h"


struct CandidateRegion {
	std::string	chromosome;
	int			ixFrom;
	int			ixTo;
	int			score;
};

void analyze(std::vector<ExpMap> expMap, std::vector<RMRead> refMap);
void init_index(RefMaps &refMaps, int height);
void init_index(RefMap &refMap, int height);
//void serialize_index(RefMap &refMap, std::string fn);
//void deserialize_index(RefMap &refMap, std::string fn);
std::vector<CandidateRegion> index_get_candidates(const RefMaps &refMaps, const ExpMap &expMap, const Params &params);

#endif // INDEXING_H