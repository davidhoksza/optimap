/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#include "types.h"
#include "constants.h"
#include "indexing.h"

#include <deque>
#include <iostream>
#include <algorithm>

using namespace std;

map<int, vector<IndexRecord> > init_index(RefMaps refMaps)
{
	map<int, vector<IndexRecord> > index;
	deque<IndexRecord> runs; //the queue contains sum of the last MAX_OPT_MAP_WINDOW reference reads

	//TODO: multiple chromosomes
	/*cout << "Identification of possible fragment lenghts in reference map..." << endl;	
	for (int ix = 0; ix < refMap.size(); ++ix)
	{
		int actualLength = refMap[ix].length;
		for (int ixQ = 0; ixQ < runs.size(); ixQ++)
		{
			runs[ixQ].end_position = ix;
			runs[ixQ].length += actualLength;
		}
		runs.push_back(IndexRecord(ix, ix, actualLength));
		for (int ixQ = 0; ixQ < runs.size(); ixQ++)
		{
			if (index.find(runs[ixQ].length) == index.end()) index[runs[ixQ].length] = vector<IndexRecord>();
			index[runs[ixQ].length].push_back(runs[ixQ]);
		}
		if (runs.size() > MAX_OPT_MAP_WINDOW) runs.pop_front();
	}*/

	return index;
}

vector<IndexRecord> index_get_candidates(map<int, vector<IndexRecord> > &index, const ExpMap &expMap, const int threshold)
{
	vector<IndexRecord> candidates;

	map<int, vector<IndexRecord> >::iterator l = index.lower_bound(expMap.length - threshold);
	map<int, vector<IndexRecord> >::iterator u = index.upper_bound(expMap.length + threshold);

	while (l != u)
	{
		candidates.insert(candidates.end(), l->second.begin(), l->second.end());
		l++;
	}

	struct {
		bool operator()(IndexRecord a, IndexRecord b)
		{
			return a.start_position < b.start_position;
		}
	} ixPosLess;
	sort(candidates.begin(), candidates.end(), ixPosLess);
	//sort(candidates.begin(), candidates.end(), [](IndexRecord & a, IndexRecord & b) -> bool	{return a.start_position < b.start_position;});

	return candidates;
}

void analyze(vector<ExpMap> expMap, RefMaps refMaps)
{
	map<int, vector<IndexRecord> > refLengths = init_index(refMaps);

	cout << "Frequencies of the fragments in the reference map:" << endl;
	for (int ixOM = 0; ixOM < expMap.size(); ixOM++)
	{
		int cnt = 0;
		map<int, vector<IndexRecord> >::iterator l = refLengths.lower_bound(expMap[ixOM].length - INDEX_NEIGHBORHOOD_THRESHOLD);
		map<int, vector<IndexRecord> >::iterator u = refLengths.upper_bound(expMap[ixOM].length + INDEX_NEIGHBORHOOD_THRESHOLD);
		while (l != u)
		{
			cnt += l->second.size();
			l++;
		}

		cout << ixOM << ": " << cnt << endl;
	}
}