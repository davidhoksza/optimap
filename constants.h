/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>

const short int MAX_OPT_MAP_WINDOW = 6;
const short int INDEX_NEIGHBORHOOD_THRESHOLD = 100;
static const float SUB_MAX = std::numeric_limits<float>::max() / 1000;
const int CNT_PROB_BINS = 100000;
const int MAX_FRAGMENT_LENGTH = 5000000;
const std::string EXPERIMENT_FORMAT_TYPES = "opgen,bionano";
const std::string ERROR_MODELS = "valuev,li,valuev-lr";

#endif // CONSTANTS_H
