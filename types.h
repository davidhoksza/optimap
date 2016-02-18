/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <map>
#include <string>
#include <limits>

struct DpMatrixCell {
	int value;
	int sourceRow, sourceColumn;

	DpMatrixCell() : value(0), sourceRow(-1), sourceColumn(-1) {};
};

struct Fragment {
	int							length;
	std::vector<int>			reads;
	std::vector<std::string>	debugInfo;
	std::string					name;

	Fragment() : length(0) {};
};

struct RMRead {
	int		start;
	int		length;
	std::string	chromosome;
};

typedef std::map <std::string, std::vector<RMRead>> RefMaps; //one map per chromosome

//one optical map fragment can optimally map to multiple positions in the reference map
//each mapping consists of a vector of mapped positions, i.e. a pair
//the mapped positions do not need to be pair of indeces, but in case when the restrictions enzyme
//misses a position, multiple indeces in the OM fragment can map to single position 
//in the reference map (and possibly vice versa)

typedef std::pair <int, int> OmRmMatch;
typedef std::vector<OmRmMatch> OmRmPath;

struct Mapping {		
	OmRmPath	alignment;
	int			score;		//score from the dp
	double		quality;	//phred-scaled probablity (-10log_10(P(wrong mapping)) that the mapping can be mapped to wrong position in the reference
	std::string chromosome; 

	Mapping(int _score, OmRmPath _alignment) : alignment(_alignment), score(_score) {};
	void ComputeQuality() {
		quality = score;
	}
};
typedef std::vector<Mapping> Mappings;

struct Params {
	std::string	omFileName;
	std::string	rmFileName;
	std::string	outFileName;
	std::string	logFileName;
	std::string	chromosome;
	int		ixOmStart;
	int		ixOmEnd;
	int		cntThreads;
	int		topK;

	int		mapOmMissedPenalty;
	int		mapRmMissedPenalty;
	int		mapDpWindowSize;

	//Params() : logFileName(NULL), cntThreads(-1) {};
};

#endif // TYPES_H