/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#include "mapping.h"

#include <sstream>

using namespace std;

ostringstream ss;

int dp_fill_matrix(DpMatrixCell ** matrix, vector<int> &fragment, std::vector<RMRead> &reference, IndexRecord ir = IndexRecord())
{
	//cout << ir.start_position << " " << ir.end_position << endl;
	if (ir.start_position == -1) ir.start_position = 1;
	else ir.start_position++; //the start position is 0-based but here we have to account for the first zero coumn in the DP matrix
	if (ir.end_position == -1) ir.end_position = reference.size();
	else ir.end_position++; //the start position is 0-based but here we have to account for the first zero coumn in the DP matrix

	//first, let's initiliaze the first column with submax values which ensures 
	//that the resulting mapping will capture the whole fragemnt
	//we don't use max values because that might cause overflow
	//since we add to these values
	for (int ixRow = 1; ixRow <= fragment.size(); matrix[ixRow++][ir.start_position - 1].value = numeric_limits<int>::max() / 1000);

	int minMappingValue = numeric_limits<int>::max();

	for (int ixRow = 1; ixRow < fragment.size() + 1; ++ixRow)
	{
		bool isLastRow = ixRow == fragment.size() ? true : false;

		//we need to compute the first column in the last row where it makes sense to search for mins
		//for example the first column does not make sense since that would mean the whole OM fragment was
		//aligned with one read in the reference map
		//int ixColResultFrom = ceil(fragment.size() / (float) DP_WINDOW_SIZE);

		for (int ixCol = ir.start_position; ixCol <= ir.end_position; ++ixCol)
		{
			//matrix[ixRow][ixCol].value = matrix[ixRow - 1][ixCol - 1].value + abs(fragment[ixRow - 1] - reference[ixCol - 1]);
			DpMatrixCell minCell;
			minCell.value = numeric_limits<int>::max();

			int rowValue = 0;
			for (int ixWindowRow = 1; ixWindowRow <= params.mapDpWindowSize; ++ixWindowRow)
			{
				if (ixRow - ixWindowRow < 0) break; //should I touch position out of the array
				rowValue += fragment[ixRow - ixWindowRow];
				int colValue = 0;
				for (int ixWindowCol = 1; ixWindowCol <= params.mapDpWindowSize; ++ixWindowCol)
				{
					if (ixCol - ixWindowCol < ir.start_position - 1) break; //should I touch position out of the candidate window
					colValue += reference[ixCol - ixWindowCol].length; //since the maps and dp table are shifted by 1, this returns in the first iteration the inspected position ([ixRow,ixCol])
					int score = matrix[ixRow - ixWindowRow][ixCol - ixWindowCol].value + abs(rowValue - colValue);

					//penalty computation
					score += (ixWindowRow - 1) * params.mapOmMissedPenalty + (ixWindowCol - 1)* params.mapRmMissedPenalty;

					//cout << ixRow - ixWindowRow << ";" << ixCol - ixWindowCol << endl;
					scoresComputed++;
					if (score < minCell.value)
					{
						minCell.value = score;
						minCell.sourceColumn = ixCol - ixWindowCol;
						minCell.sourceRow = ixRow - ixWindowRow;
					}
				}
			}
			matrix[ixRow][ixCol] = minCell;
			if (isLastRow /*&& ixCol >= ixColResultFrom*/ && minCell.value < minMappingValue) minMappingValue = minCell.value;
		}
	}

	return minMappingValue;
}


Mappings dp_backtrack(DpMatrixCell **matrix, int height, int width, IndexRecord ir = IndexRecord())
{
	Mappings mappings;

	vector<pair<int, int>> minPositions; //top K min values and positions in increasing order
	if (ir.start_position == -1) ir.start_position = 1;
	else ir.start_position++;
	if (ir.end_position == -1) ir.end_position = width - 1;
	else ir.end_position++;

	for (int ixM = ir.start_position; ixM <= ir.end_position; ++ixM)
	{
		int alignmentScore = matrix[height - 1][ixM].value;
		if (minPositions.size() == 0) minPositions.push_back(make_pair(alignmentScore, ixM));
		else if (minPositions.size() < params.topK || (minPositions.end() - 1)->first > alignmentScore)
		{
			//find position in the sorted vector
			int pos = 0;
			while (pos < minPositions.size() && minPositions[pos].first < alignmentScore) pos++;
			//insert at that position
			minPositions.insert(minPositions.begin() + pos, make_pair(alignmentScore, ixM));
			//if there are more then K mappings, remove the one with lowest score
			if (minPositions.size() > params.topK) minPositions.erase(minPositions.end() - 1);
		}
	}

	//finally, we find the index with min value in the last row and extract the diagonal
	//we start at index height since we are interested only in diagonals of length at least "height"
	for (auto match : minPositions)
	{
		OmRmPath diagonal;
		//for (int ixDiagonal = height - 2; ixDiagonal >= 0; --ixDiagonal) diagonal.push_back(make_pair<int, int>(height - 2 - ixDiagonal, pos - 1 - ixDiagonal));

		int ixRow = height - 1, ixCol = match.second;
		while (ixRow >= 1 && ixCol >= ir.start_position)
		{
			diagonal.push_back(make_pair(ixRow, ixCol));
			int ixRowAux = matrix[ixRow][ixCol].sourceRow;
			int ixColAux = matrix[ixRow][ixCol].sourceColumn;
			ixRow = ixRowAux;
			ixCol = ixColAux;
		}
		//Finally, we add the last position which is not part of the alignment, but helps to identify
		//the boundaries. This is used when the alignment does not start at the beginning of the reference map (indexing)
		diagonal.push_back(make_pair(ixRow, ixCol));
		//Since the backtracking started from the end of the alignment, we reverse the result
		reverse(diagonal.begin(), diagonal.end());
		mappings.push_back(Mapping(match.first, diagonal));
	}

	return mappings;
}

Mappings do_mapping(vector<int> &fragment, std::vector<RMRead> &refMap, vector<IndexRecord> &candidates)
{
	DpMatrixCell **dpMatrix = new DpMatrixCell*[fragment.size() + 1];
	for (int ixDPM = 0; ixDPM < fragment.size() + 1; ++ixDPM) dpMatrix[ixDPM] = new DpMatrixCell[refMap.size() + 1];

	Mappings matchSequences;
	if (candidates.size() == 0)
	{
		dp_fill_matrix(dpMatrix, fragment, refMap);
		matchSequences = dp_backtrack(dpMatrix, fragment.size() + 1, refMap.size() + 1);
	}
	else
	{
		//first, we remove candidates which are too short, i.e. the number segemnts < fragment_length/DP_WINDOW_SIZE
		//second, the candidates need to be merged in case they overlap
		//for (int ixCnd = candidates.size() - 1; ixCnd >= 0; ixCnd--)
		//{
		//	int cndLength = candidates[ixCnd].end_position - candidates[ixCnd].start_position + 1;
		//	if (cndLength < fragment.size() / (float)DP_WINDOW_SIZE) candidates.erase(candidates.begin() + ixCnd);
		//}
		//ss << "#candidates after removing short candidates: " << candidates.size() << endl; logger.Log(Logger::LOGFILE, ss);

		//vector<IndexRecord> candidatesMerged;
		//candidatesMerged.push_back(candidates[0]);		
		//for (int ixCnd = 1; ixCnd < candidates.size(); ixCnd++)
		//{
		//	//we merge even if the stretches adjoin
		//	if (candidates[ixCnd].start_position <= candidates[ixCnd - 1].end_position + 1) (candidatesMerged.end() - 1)->end_position = candidates[ixCnd].end_position;
		//	else candidatesMerged.push_back(candidates[ixCnd]);
		//}
		//int lengthSum = 0;
		//for (int ixCnd = 0; ixCnd < candidatesMerged.size(); ixCnd++) lengthSum += candidatesMerged[ixCnd].end_position - candidatesMerged[ixCnd].start_position;
		//candidateSectionLengths.push_back(lengthSum);
		//ss << "#candidates after meging: " << candidatesMerged.size() << endl; logger.Log(Logger::LOGFILE, ss);		

		//vector<pair<int, int> > bestMapping; //indexes and values of best candidate (there can be multiple mappings with minimal value)
		//bestMapping.push_back(pair<int, int>(0, dp_fill_matrix(dpMatrix, fragment, refMap, candidatesMerged[0])));		
		//for (int ixCnd = 1; ixCnd < candidatesMerged.size(); ixCnd++)
		//{
		//	int aux = dp_fill_matrix(dpMatrix, fragment, refMap, candidatesMerged[ixCnd]);
		//	if (aux < bestMapping[0].second) {
		//		bestMapping.clear();
		//		bestMapping.push_back(pair<int, int>(ixCnd, aux));
		//	}
		//	else if (aux == bestMapping[0].second) bestMapping.push_back(pair<int, int>(ixCnd, aux));
		//}
		////cout << "mappings: " << bestMapping.size() << endl;		
		//matchSequences.value = bestMapping[0].second;
		//for (int ixM = 0; ixM < bestMapping.size(); ixM++)
		//{
		//	Mappings aux = dp_backtrack(dpMatrix, fragment.size() + 1, refMap.map.size() + 1, candidatesMerged[bestMapping[ixM].first], bestMapping[0].second);
		//	matchSequences.matches.insert(matchSequences.matches.end(), aux.matches.begin(), aux.matches.end());
		//}
		//cout << "mapped" << endl;
	}

	for (int i = 0; i < fragment.size() + 1; ++i) delete[] dpMatrix[i];
	delete[] dpMatrix;

	return matchSequences;
}

void verify_candidates(Fragment &optMap, vector<int> &refMap, vector<IndexRecord> candidates)
{
	int omSize = optMap.length;
	for (int iC = 0; iC < candidates.size(); iC++)
	{
		ss << "Verifying candidate (start: " << candidates[iC].start_position << ", end: " << candidates[iC].end_position
			<< ", length:" << candidates[iC].length << ")" << endl; logger.Log(Logger::STDOUT, ss);

		int rmSize = 0;
		for (int i = candidates[iC].start_position; i <= candidates[iC].end_position; i++) rmSize += refMap[i];

		ss << "om: " << omSize << ", "; logger.Log(Logger::STDOUT, ss);
		ss << "rf: " << rmSize << ",  "; logger.Log(Logger::STDOUT, ss);
		ss << "diff: " << abs(omSize - rmSize) << endl; logger.Log(Logger::STDOUT, ss);
	}

}

void map_segment(int from, int to, vector<Fragment> &optMap, RefMaps &refMaps, Mappings* resultSet, map<int, vector<IndexRecord> > &index)
{
	for (int ixOM = from; ixOM <= to; ++ixOM)
	{
		int threshold = INDEX_NEIGHBORHOOD_THRESHOLD;
		vector<IndexRecord> candidates;
		while (!index.empty() && candidates.empty())
		{
			threshold++;
			candidates = index_get_candidates(index, optMap[ixOM], threshold);
		}
		mutexOM.lock();
		ss << "Thresold: " << threshold << ". Candidates: " << candidates.size() << endl; logger.Log(Logger::LOGFILE, ss);
		mutexOM.unlock();

		Mappings mappings;
		for each (auto refMap in refMaps)
		{
			//resultSet[ixOM] = do_mapping(optMap[ixOM].reads, refMap.second, candidates);
			Mappings aux_mapping = do_mapping(optMap[ixOM].reads, refMap.second, candidates);
			for (auto it = aux_mapping.begin(); it != aux_mapping.end(); ++it){
				it->chromosome = refMap.first;
				it->ComputeQuality();
			}
			mappings.insert(mappings.end(), aux_mapping.begin(), aux_mapping.end());
		}

		//keep top params.topK mappings
		sort(mappings.begin(), mappings.end(), [](Mapping & a, Mapping & b) -> bool	{return a.score < b.score; });
		mappings.erase(mappings.begin() + params.topK, mappings.end());
		resultSet[ixOM] = mappings;

		//verify_candidates(optMap[ixOM], refMap, candidates);		
		scoresCalculations.push_back(scoresComputed); scoresComputed = 0;

		mutexOM.lock();
		omProcessed++;
		ss << omProcessed << " "; logger.Log(Logger::STDOUT, ss);
		mutexOM.unlock();
	}
}

Mappings* AlignOpticalMaps(vector<Fragment> &optMap, RefMaps &refMaps, map<int, vector<IndexRecord> > &index)
{
	ss << "======= ALIGNING - START =======" << endl; logger.Log(Logger::STDOUT, ss);
	clock_t begin_time = clock();
	//initialize a vector of results where the threads will store the mappings of individual reads
	//after all the opt maps will be mapped, the vector will be serialized
	Mappings *omMappings = new Mappings[optMap.size()];

	int cntThreads = thread::hardware_concurrency();
	if (params.cntThreads > 0) cntThreads = params.cntThreads;
	ss << "Using " << cntThreads << " threads" << endl; logger.Log(Logger::STDOUT, ss);
	int batchSize = ceil(optMap.size() / (double)cntThreads);
	assert(batchSize > 0);
	vector<thread> threads;

	int ixFrom = 0, ixTo = optMap.size() - 1;
	if (params.ixOmStart > 0) ixFrom = params.ixOmStart;
	if (params.ixOmEnd >= 0) ixTo = params.ixOmEnd;

	while (ixFrom <= ixTo)
	{
		int ixToAux = ixFrom + batchSize;
		if (ixToAux > ixTo) ixToAux = ixTo;
		threads.push_back(thread(map_segment, ixFrom, ixToAux, std::ref(optMap), std::ref(refMaps), omMappings, std::ref(index)));
		ixFrom = ixToAux + 1;
	}
	for (int i = 0; i < threads.size(); i++) threads[i].join();

	ss << endl << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl; logger.Log(Logger::STDOUT, ss);
	ss << "======= ALIGNING - END =======" << endl; logger.Log(Logger::STDOUT, ss);

	return omMappings;
}