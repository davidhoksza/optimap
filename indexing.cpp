/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#include "types.h"
#include "constants.h"
//#include "common.h"
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
		//cout << ixHeight << endl;
		for (int ixForrest = refMap.sumForrest.size() - 1; ixForrest >= 1; ixForrest--)
		{
			SumTree *st = new SumTree();
			st->left = refMap.sumForrest[ixForrest - 1];
			st->right = refMap.sumForrest[ixForrest];	
			st->height = st->left->height + 1;
			st->sum = st->GetSum();
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

//void serialize_index(RefMap &refMap, std::string fn)
//{	
//	ofstream ofs(fn);
//	string errorMsg = "Index map file " + fn + " could not be created";
//	if (ofs.is_open()) error_exit(errorMsg);
//
//
//}

void merge_candidate_regions(vector<CandidateRegion> &crs)
{
	for (int ixCR1 = crs.size() - 1; ixCR1 > 0; --ixCR1)
	{
		for (int ixCR2 = ixCR1 - 1; ixCR2 >= 0; --ixCR2)
		{
			if ((crs[ixCR1].ixFrom >= crs[ixCR2].ixFrom && crs[ixCR1].ixFrom <= crs[ixCR2].ixTo) ||
				(crs[ixCR1].ixTo >= crs[ixCR2].ixFrom && crs[ixCR1].ixTo <= crs[ixCR2].ixTo))
			{
				crs[ixCR2].ixFrom = min(crs[ixCR2].ixFrom, crs[ixCR1].ixFrom);
				crs[ixCR2].ixTo = max(crs[ixCR2].ixTo, crs[ixCR1].ixTo);
				crs[ixCR2].score = max(crs[ixCR2].score, crs[ixCR1].score);

				crs.erase(crs.begin() + ixCR1);
				break;
			}
		}
	}
}

vector<CandidateRegion> _index_get_candidates(const RefMap &refMap, const ExpMap &expMap, const Params &params)
{
	vector<CandidateRegion> candidates;

	int optimalHeigt = expMap.reads.size() * (1 + params.missRestrictionProb) - 1; //-1 because node at heigh 1 has 2 childern
	SumTree *st = refMap.sumForrest[0]->GetMostLeftAtHeight(optimalHeigt);

	while (st)
	{
		int diff = abs(st->sum - expMap.length);
		if (candidates.size() < params.topK || diff < candidates[0].score)
		{
			CandidateRegion cr;
			cr.chromosome = refMap.fragments[0].chromosome;
			cr.score = diff;
			cr.ixFrom = st->ixFrom;
			cr.ixTo = st->ixTo;

			if (candidates.size() > params.topK) candidates.erase(candidates.begin());

			if (candidates.size() == 0) candidates.push_back(cr);
			else
			{
				vector<CandidateRegion>::iterator it = candidates.begin();
				while (it != candidates.end() && it->score < diff) it++;
				candidates.insert(it, cr);
			}
		}

		st = st->rightSibling;
	}

	merge_candidate_regions(candidates);

	return candidates;
}

vector<CandidateRegion> index_get_candidates(const RefMaps &refMaps, const ExpMap &expMap, const Params &params)
{
	vector<CandidateRegion> candidates;
	for (RefMaps::const_iterator refMap = refMaps.begin(); refMap != refMaps.end(); ++refMap)
	{
		vector<CandidateRegion> auxCandidates = _index_get_candidates(refMap->second, expMap, params);
		candidates.insert(candidates.end(), auxCandidates.begin(), auxCandidates.end());
	}
	
	struct { bool operator()(CandidateRegion a, CandidateRegion b){ return a.score < b.score; } } crLess;
	sort(candidates.begin(), candidates.end(), crLess);
	if (params.topK < candidates.size()) candidates.erase(candidates.begin() + params.topK, candidates.end());

	for (int i = 0; i < candidates.size(); i++) cout << candidates[i].score << " ";
	cout << endl;

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