/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef MAPPING_H
#define MAPPING_H

#include "types.h"
#include "indexing.h"

Mappings* AlignOpticalMaps(std::vector<Fragment> &optMap, RefMaps &refMaps, std::map<int, std::vector<IndexRecord> > &index);

#endif // MAPPING_H