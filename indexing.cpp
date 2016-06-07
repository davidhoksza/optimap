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

void init_index(RefMaps &refMaps, int height)
{
	for (RefMaps::iterator refMap = refMaps.begin(); refMap != refMaps.end(); ++refMap) init_index(refMap->second, height);
}

void link_siblings(SumForrest &sf)
{
	for (int ixForrest = 0; ixForrest < sf.size() - 1; ++ixForrest) sf[ixForrest]->rightSibling = sf[ixForrest + 1];

}

void init_index(RefMap &refMap, int height)
{
	for (int i = 0; i < refMap.fragments.size(); i++) 
	{
		refMap.sumForrest.push_back(new SumTree(NULL, NULL, 0, refMap.fragments[i].length, i, i));
	}
	link_siblings(refMap.sumForrest);

	for (int ixHeight = 1; ixHeight < height; ++ixHeight)
	{
		cout << ixHeight << endl;
		for (int ixForrest = refMap.sumForrest.size() - 1; ixForrest >= 1; ixForrest--)
		{
			SumTree *st = new SumTree();
			st->left = refMap.sumForrest[ixForrest - 1];
			st->right = refMap.sumForrest[ixForrest];			
			st->sum = st->get_sum();
			st->ixFrom = st->left->ixFrom;
			st->ixTo = st->right->ixTo;

			refMap.sumForrest.insert(refMap.sumForrest.begin() + ixForrest, st);
			refMap.sumForrest.erase(refMap.sumForrest.begin() + ixForrest + 1);
		}
		refMap.sumForrest.erase(refMap.sumForrest.begin());
		refMap.sumForrest[0]->mostLeft = true;
		if (refMap.sumForrest.size() == 1) break;
		link_siblings(refMap.sumForrest);
	}	
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
		bool operator()(IndexRecord a, IndexRecord b) {
			return a.start_position < b.start_position;
		}
	} ixPosLess;
	sort(candidates.begin(), candidates.end(), ixPosLess);
	//sort(candidates.begin(), candidates.end(), [](IndexRecord & a, IndexRecord & b) -> bool	{return a.start_position < b.start_position;});

	return candidates;
}

void analyze(vector<ExpMap> expMap, RefMaps refMaps)
{
	/*map<int, vector<IndexRecord> > refLengths = init_index(refMaps);

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
	}*/
}