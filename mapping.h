/*
* Copyright (C) 2017 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef MAPPING_H
#define MAPPING_H

#include "common.h"

void clean_dp_matrix(DpMatrixCell ** matrix, int height, int width);

inline SCORE_TYPE transform_prob(SCORE_TYPE p);

SCORE_TYPE score_segment(int expLength, int refLength, int cntExpFrags, int cntRefFrags, int ixExp, int ixRef, std::vector<int> &experiment, std::vector<RMRead> &reference);

void dp_fill_matrix(DpMatrixCell ** matrix, std::vector<int> &experiment, std::vector<RMRead> &reference, std::vector<SCORE_TYPE> &minScoresSoFar);

Mappings dp_backtrack(DpMatrixCell **matrix, int height, int width);

Mappings do_mapping(std::vector<int> &expMap, std::vector<RMRead> &refMap, std::vector<SCORE_TYPE> &minScoresSoFar);

void map_segment(int from, int to, std::vector<ExpMap> &expMap, RefMaps &refMaps, Mappings* resultSet);

Mappings* AlignOpticalMaps(std::vector<ExpMap> &expMap, RefMaps &refMaps);

#endif // MAPPING_H