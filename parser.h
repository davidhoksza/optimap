/*
* Copyright (C) 2017 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef PARSE_H
#define PARSE_H

#include "common.h"

std::istream* open_map_file(std::string fileName);
RefMaps parse_ref_map(std::string fileName = params.rmFileName/*, float &avgRefLength*/);
std::vector<ExpMap> parse_exp_map(std::string fileName = params.omFileName, int topN = std::numeric_limits<int>::max());
void SmoothExpFragments(std::vector<ExpMap> &expMaps);
void SmoothRefFragments(RefMaps &refMaps);
void Parse(std::vector<ExpMap> &expMaps, RefMaps &refMaps);


#endif // PARSE_H