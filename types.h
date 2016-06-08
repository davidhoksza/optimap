/*
* Copyright (C) 2016 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <string>
#include <limits>
#include <map>

#include "constants.h"

struct DpMatrixCell {
	float value;
	int sourceRow, sourceColumn;

	DpMatrixCell() : value(SUB_MAX), sourceRow(-1), sourceColumn(-1) {};
	void flush() { value = SUB_MAX; sourceRow = -1; sourceColumn = -1; };
};

struct ExpMap {
	int							length = 0;
	std::vector<int>			reads;
	std::vector<float>			qx11; //bionano qualities for each label (restriction site)
	std::vector<float>			qx12; //bionano noise-to-signal ratios for each label (restriction site)
	std::vector<std::string>	debugInfo;
	std::string					name;

	void Clear() {
		length = 0;
		reads.clear();
		debugInfo.clear();
		qx11.clear();
		qx12.clear();
		name = "";
	}
};

struct RMRead {
	int		start;
	int		length;
	std::string	chromosome;
};


struct SumTree {
	SumTree *left = NULL;
	SumTree *right = NULL;

	SumTree *rightSibling = NULL;

	bool	mostLeft = false; //used when deallocating where it is used to deallocate a node only once

	int		height = 0;
	int		sum = -1; //sum of the values in the interval the tree represents
	int		ixFrom; //start index in the map the tree represents
	int		ixTo; //end index in the map the tree represents

	SumTree() {};
	SumTree(SumTree *_left, SumTree *_right, int _height, int _sum, int _ixFrom, int _ixTo) :
		left(_left), right(_right), height(_height), sum(_sum), ixFrom(_ixFrom), ixTo(_ixTo) {};

	~SumTree()
	{
		if (mostLeft && left) delete left;
		if (right) delete right;
	}

	SumTree *GetMostLeftLeaf()
	{
		return left ? left->GetMostLeftLeaf() : this;
	}

	SumTree *GetMostRightLeaf()
	{
		return right ? right->GetMostRightLeaf() : this;
	}

	SumTree *GetMostLeftAtHeight(int _height)
	{
		if (height == 0) return NULL;
		else return _height == height ? this : left->GetMostLeftAtHeight(_height);
			
	}

	int GetSum()
	{
		SumTree *stl = GetMostLeftLeaf();
		SumTree *str = GetMostRightLeaf();

		int _sum = 0;

		SumTree *st = stl;
		while (st != str)
		{
			_sum += st->sum;
			st = st->rightSibling;
		}
		_sum += st->sum;

		return _sum;
	}	

};

typedef std::vector<SumTree *> SumForrest;

struct RefMap {
	std::vector<RMRead> fragments;
	SumForrest	sumForrest;

	~RefMap()
	{
		for (int ix = 0; ix < sumForrest.size(); ix++) delete sumForrest[ix];
	}

};

typedef std::map <std::string, RefMap> RefMaps; //one map per chromosome



//one optical map fragment can optimally map to multiple positions in the reference map
//each mapping consists of a vector of mapped positions, i.e. a pair
//the mapped positions do not need to be pair of indeces, but in case when the restrictions enzyme
//misses a position, multiple indeces in the OM fragment can map to single position 
//in the reference map (and possibly vice versa)

typedef std::pair <int, int> OmRmMatch;
typedef std::vector<OmRmMatch> OmRmPath;

typedef float SCORE_TYPE;

struct Mapping {		
	OmRmPath	alignment;
	std::vector<SCORE_TYPE> scores;
	SCORE_TYPE	score;		//score from the dp
	float		quality;	//phred-scaled probablity (-10log_10(P(wrong mapping)) that the mapping can be mapped to wrong position in the reference
	std::string chromosome; 
	bool		reversed = false;

	Mapping(float _score, std::vector<SCORE_TYPE> _scores, OmRmPath _alignment) : alignment(_alignment), score(_score), scores(_scores) {};
	void ComputeQuality() {
		quality = score;
	}
};

typedef std::vector<Mapping> Mappings;

struct Params {
	std::string	omFileName;
	std::string	omFormat;
	std::string	rmFileName;
	std::string	outFileName;
	std::string	logFileName;
	std::string	chromosome;
	int		ixOmStart;
	int		ixOmEnd;
	int		cntThreads;
	int		topK;
	//int		mapOmMissedPenalty;
	//int		mapRmMissedPenalty;
	int		maxDpWindowSize;
	float	sizingErrorStddev;
	int		smallFragmentThreshold;
	float	missRestrictionProb;
	float	noMissRestrictionProb; // computed as 1 - ((1 - pow(params.missRestrictionProb, params.maxDpWindowSize - 1)) / (1 - params.missRestrictionProb) - 1)
	float	falseCutProb;
	int		smoothingThreshold;

};

#endif // TYPES_H