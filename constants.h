/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "types.h"

const short int MAX_OPT_MAP_WINDOW = 10;
const short int INDEX_NEIGHBORHOOD_THRESHOLD = 100;
static const SCORE_TYPE SUB_MAX = std::numeric_limits<SCORE_TYPE>::max() / 1000;
const int CNT_PROB_BINS = 1000000;
const int MAX_FRAGMENT_LENGTH = 200000;

#endif // CONSTANTS_H