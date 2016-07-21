/*
* Copyright (C) 2016 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#define _SCL_SECURE_NO_WARNINGS

#include "common.h"
#include "constants.h"
//#include "indexing.h"
#include "types.h"
#include "logger.h"
#include "string_functions.h"
#include "stats.h"

#include "tclap/CmdLine.h"

#include "gzstream/gzstream.h"

#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>  

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <thread>
#include <mutex>
#include <sstream>
#include <algorithm>

using namespace std;

mutex mutexOM;
long unsigned int omProcessed = 0;

ostringstream ss;
vector<int> candidateSectionLengths;
vector<int> scoresCalculations;
unsigned long int scoresComputed = 0;

Params params;

istream* open_map_file(string fileName)
{
	string errorMsg = "Reference map file " + fileName + " could not be opened";
	istream *ifs;
	if (strings::ends_with(fileName, ".gz"))
	{
		ifs = new igzstream(fileName.c_str());
		if (!((igzstream*)ifs)->is_open()) error_exit(errorMsg);
	}
	else {
		ifs = new ifstream(fileName);
		if (!((ifstream*)ifs)->is_open()) error_exit(errorMsg);
	}
	if (!ifs->good()) error_exit(errorMsg);

	return ifs;
}

RefMaps parse_ref_map(string fileName)
{
	RefMaps refMaps;

	istream *ifs = open_map_file(fileName);

	string line;
	for (string line; getline(*ifs, line);) {
		vector<string> strs = strings::split(line, "\t");

		RMRead auxRMRead;
		auxRMRead.chromosome = strings::trim(strs[0]);
		if (params.chromosome != "" && strings::upper(auxRMRead.chromosome) != strings::upper(params.chromosome)) continue;
		auxRMRead.start = stof(strs[1]);
		auxRMRead.length = stof(strs[3]) * 1000;
		if (refMaps.count(auxRMRead.chromosome) == 0) refMaps[auxRMRead.chromosome] = vector<RMRead>();
		refMaps[auxRMRead.chromosome].push_back(auxRMRead);
	}

	delete ifs;
	return refMaps;
}

vector<ExpMap> parse_exp_map(string fileName, int topN = numeric_limits<int>::max())
{
	vector<ExpMap> expMap;

	istream *ifs = open_map_file(fileName);

	string line;
	char buffer[20];
	vector<string> cleavageSites;

	string fragName;
	ExpMap f;
	for (std::string line; getline(*ifs, line);)
	{
		if (params.omFormat == "opgen")
		{
			if (line.find("debug") != string::npos)
			{
				stringstream ss(line); // Insert the string into a stream
				string strBuffer;
				ss >> strBuffer;
				while (ss >> strBuffer) cleavageSites.push_back(strBuffer);
			}

			if (line.find("KpnI") != string::npos)
			{			
				f.Clear();
				f.name = fragName;

				int pos = -1;

				for (string::const_iterator it = line.begin(), end = line.end(); it != end; ++it)
				{
					pos++;
					if (*it == '\t')
					{
						buffer[pos] = 0;
						if (buffer[0] != 'K' && pos > 0) {
							int aux = 1000 * atof(buffer);
							f.reads.push_back(aux);
							f.length += aux;
						}
						pos = -1;
					}
					else buffer[pos] = *it;
				}
				if (pos >= 0)
				{
					buffer[pos + 1] = 0;
					int aux = 1000 * atof(buffer);
					f.reads.push_back(aux);
					f.length += aux;
				}

				if (cleavageSites.size() > 0)
				{
					assert(cleavageSites.size() - 2 == f.reads.size()); // debug info includes chromosome and first and last position -> +2
					f.debugInfo = cleavageSites;
					cleavageSites.clear();
				}
				expMap.push_back(f);
				if (expMap.size() >= topN) break;
			}
			else fragName = strings::trim(line);
		}
		if (params.omFormat == "bionano")
		{
			vector<string> sLine = strings::split(line, "\t");

			if (sLine[0] == "0")
			{
				f.Clear();
				ss << sLine[1] << "_" << sLine[2];
				f.name = ss.str(); ss.str(string());				
			}
			if (sLine[0] == "1")
			{
				f.reads.push_back(strings::toFloat(sLine[0]));
				int s = sLine.size();
				for (int ix = 2; ix < s; ix++)
				{
					f.reads.push_back(strings::toInt(sLine[ix]) - strings::toInt(sLine[ix - 1]));
				}
				f.length = strings::toFloat(sLine[s - 1]);
			}
			if (sLine[0] == "QX11")
			{
				for (int ix = 1; ix < sLine.size(); f.qx11.push_back(strings::toFloat(sLine[ix++])));
			}
			if (sLine[0] == "QX12")
			{
				for (int ix = 1; ix < sLine.size(); f.qx12.push_back(strings::toFloat(sLine[ix++])));
				expMap.push_back(f);
				if (expMap.size() >= topN) break;
			}
		}
	}

	delete ifs;
	return expMap;
}

void clean_dp_matrix(DpMatrixCell ** matrix, int height, int width)
{
	for (int ixRow = 0; ixRow < height; ixRow++)
		for (int ixCol = 0; ixCol < width; ixCol++)
			matrix[ixRow][ixCol].flush();
}

inline SCORE_TYPE transform_prob(SCORE_TYPE p)
{
	if (p == 0) return SUB_MAX;
	SCORE_TYPE aux = -log(p);
	return (aux > SUB_MAX) ? SUB_MAX : aux;
}

SCORE_TYPE score_segment(int expLength, int refLength, int cntExpFrags, int cntRefFrags, int ixExp, int ixRef, vector<int> &experiment, std::vector<RMRead> &reference)
{
	SCORE_TYPE score = 0;

	if (params.errorModel == "valuev")
	{

		//float stddev = params.sizingErrorStddev * (exp_length > params.smallFragmentThreshold ? exp_length : params.smallFragmentThreshold);
		float stddev = sqrt(float(refLength)) * (refLength > params.smallFragmentThreshold ? 5 : 6.4);

		float x;
		if (stddev == 0)
		{
			if (expLength == refLength) x = 1;
			else x = 0;
		}
		else x = stats::pdf_gaussian((expLength - refLength) / stddev, 0, 1);
		//else x = stats::pdf_gaussian_full((exp_length - ref_length) / stddev, 0, 1);
		//else x = stats::pdf_gaussian_full((exp_length - ref_length) / sqrt((double)ref_length), 0, 5);

		score += stats::transform_prob(x);

		//penalty computation
		//score += (ixWindowRow - 1) * params.mapOmMissedPenalty + (ixWindowCol - 1)* params.mapRmMissedPenalty;
		SCORE_TYPE auxP;
		if (params.maxDpWindowSize == 1) auxP = 1;
		else
		{
			if (cntRefFrags > 1) auxP = pow(params.missRestrictionProb, cntRefFrags - 1);
			else auxP = params.noMissRestrictionProb;
		}
		score += stats::transform_prob(auxP);

		//x = stats::pdf_poisson_full(ixWindowRow - 1, colValue * params.falseCutProb);
		// The Poisson PDF is precomputed with rowValue * params.falseCutProb. But since falseCutProb is constant,
		// we can provide only the rowValue and use it as index to the array with precomputed values.
		x = stats::pdf_poisson(cntExpFrags - 1, expLength);
		score += stats::transform_prob(x);
	}
	else if (params.errorModel == "li")
	{
		// Sizing error

		float location, scale;
		if (refLength < 2400)
		{
			location = 0.858181;
			scale = 0.180196;
		}
		else if (refLength < 3600)
		{
			location = 0.980760;
			scale = 0.071176;
		}
		else if (refLength < 4800)
		{
			location = 1.003354;
			scale = 0.052800;
		}
		else 
		{
			location = 1.00482;
			scale = 0.042428;
		}

		score += stats::transform_prob(stats::pdf_laplace_full(expLength/(float)refLength, location, scale));

		// Missing cuts + aligned cut
		float digestion_rate;
		SCORE_TYPE auxP = 0;
		for (int ix = 1; ix <= cntRefFrags; ix++)
		{
			if (ixRef + ix == reference.size())
			{
				digestion_rate = 1;
				
			}
			else
			{
				int d1 = reference[ixRef + ix - 1].length;
				int d2 = reference[ixRef + ix].length;
				float dAvg = (d1 + d2) / (2.0 * 1200);
				digestion_rate = 0.0003089 * dAvg * dAvg * dAvg - 0.01069 * dAvg * dAvg + 0.1253 * dAvg + 0.3693;				
			}
			if (ix < cntRefFrags) auxP *= 1 - digestion_rate;
			else auxP *= digestion_rate;				
		}
		score += stats::transform_prob(auxP);

		//False cuts
		//probability of seeing given number of false cuts
		auxP = 0.18 * stats::pdf_poisson_full(cntExpFrags - 1, 0) + 0.6 * stats::pdf_poisson_full(cntExpFrags - 1, 1) + 0.22 * stats::pdf_poisson_full(cntExpFrags - 1, 3);
		//now this needs to be modified based on where the given cuts are
		//TODO
	}

	return score;
}

void dp_fill_matrix(DpMatrixCell ** matrix, vector<int> &experiment, std::vector<RMRead> &reference, vector<SCORE_TYPE> &minScoresSoFar)
{
	//first, let's initiliaze the first column with submax values which ensures 
	//that the resulting mapping will capture the whole fragemnt
	//we don't use max values because that might cause overflow
	//since we add to these values
	for (int ixRow = 1; ixRow <= experiment.size(); matrix[ixRow++][0].value = SUB_MAX);

	//and first row by 0 so that the alignment can't start anywhere in the reference
	for (int ixCol = 0; ixCol <= reference.size(); matrix[0][ixCol++].value = stats::transform_prob(1));

	for (int ixRow = 1; ixRow < experiment.size() + 1; ++ixRow)
	{
		bool isLastRow = ixRow == experiment.size() ? true : false;

		for (int ixCol = 1; ixCol <= reference.size(); ++ixCol)
		{
			DpMatrixCell minCell;
			minCell.value = SUB_MAX;
			minCell.sourceColumn = ixCol - 1; //This value is used when no good match is found
			minCell.sourceRow = ixRow - 1; //This value is used when no good match is found
			
			//First and last fragments are scored 0 AND experiment fragment is shorter 
			//(otherwise we would omit a clear cut, i.e. the alignment would be on wrong place due to this first/last fragment)
			if ((ixRow == 1 || isLastRow) && experiment[ixRow - 1] - reference[ixCol - 1].length < 0)
			{
				isLastRow ? minCell.value = matrix[ixRow - 1][ixCol - 1].value + stats::transform_prob(1) : minCell.value = stats::transform_prob(1);
			}
			else
			{
				//check whether the column and row are not too different				
				//float stddev = params.sizingErrorStddev * (experiment[ixRow] > params.smallFragmentThreshold ? experiment[ixRow] : params.smallFragmentThreshold);
				//float ds = diff / stddev;
				//if (ds < -stats::max__reasonable_stddev || ds > stats::max__reasonable_stddev) continue;

				int rowValue = 0;
				for (int ixWindowRow = 1; ixWindowRow <= params.maxDpWindowSize; ++ixWindowRow)
				{
					int ixExp = ixRow - ixWindowRow;
					if (ixRow - ixWindowRow < 0) break; //should I touch position out of the array
					
					//ends of fragments can be aligned with zero score since
					//the molecules forming fragments were not created with a restriction enzyme
					rowValue += experiment[ixExp];
					int colValue = 0;
					for (int ixWindowCol = 1; ixWindowCol <= params.maxDpWindowSize; ++ixWindowCol)
					{
						if (ixCol - ixWindowCol < 0) break; //should I touch position out of the candidate window

						int ixRef = ixCol - ixWindowCol;
						colValue += reference[ixRef].length; //since the maps and dp table are shifted by 1, this returns in the first iteration the inspected position ([ixRow,ixCol])
						//float score = matrix[ixRow - ixWindowRow][ixCol - ixWindowCol].value + pow(rowValue - colValue, 2)/(colValue*1.05);	

						SCORE_TYPE score = matrix[ixRow - ixWindowRow][ixRef].value;
						if (score >= minScoresSoFar[0]) continue;
						//we use ixCol-1 and not ixCol, since in score_segment we will index experiment and refernece maps, which are not +1 indexed as the DP matrix
						score += score_segment(rowValue, colValue, ixWindowRow, ixWindowCol, ixCol-1, ixRow-1, experiment, reference);
						if (score >= minScoresSoFar[0]) continue;
						
						if (score < minCell.value)
						{
							minCell.value = score;
							minCell.sourceColumn = ixRef;
							minCell.sourceRow = ixExp;
						}
					}
				}
			}
			matrix[ixRow][ixCol] = minCell;
			if (isLastRow /*&& ixCol >= ixColResultFrom*/ && minCell.value < minScoresSoFar[0])
			{
				minScoresSoFar.erase(minScoresSoFar.begin());
				minScoresSoFar.push_back(minCell.value);
				sort(minScoresSoFar.begin(), minScoresSoFar.end(), std::greater<SCORE_TYPE>());
				//for (auto &s : minScoresSoFar) cout << s << " "; cout << endl;
			}
		}
	}
}

Mappings dp_backtrack(DpMatrixCell **matrix, int height, int width)
{
	Mappings mappings;

	vector<pair<SCORE_TYPE, int>> minPositions; //top K min values and positions in increasing order

	for (int ixM = 1; ixM <= width - 1; ++ixM)
	{
		SCORE_TYPE alignmentScore = matrix[height - 1][ixM].value;
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
	for (int ix = 0; ix < minPositions.size(); ix++)
	{
		pair<SCORE_TYPE, int> match = minPositions[ix];
		OmRmPath diagonal;
		vector<SCORE_TYPE> scores;

		int ixRow = height - 1, ixCol = match.second;
		while (ixRow >= 1 && ixCol >= 1)
		{
			diagonal.push_back(make_pair(ixRow, ixCol));
			scores.push_back(matrix[ixRow][ixCol].value);
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
		reverse(scores.begin(), scores.end());
		mappings.push_back(Mapping(match.first, scores, diagonal));
	}

	return mappings;
}

Mappings do_mapping(vector<int> &expMap, std::vector<RMRead> &refMap, vector<SCORE_TYPE> &minScoresSoFar)
{
	DpMatrixCell **dpMatrix = new DpMatrixCell*[expMap.size() + 1];
	for (int ixDPM = 0; ixDPM < expMap.size() + 1; ++ixDPM) dpMatrix[ixDPM] = new DpMatrixCell[refMap.size() + 1];

	Mappings matchSequences, matchSequencesRev;

	dp_fill_matrix(dpMatrix, expMap, refMap, minScoresSoFar);
	matchSequences = dp_backtrack(dpMatrix, expMap.size() + 1, refMap.size() + 1);

	//repeat the process with reversed experimental map
	clean_dp_matrix(dpMatrix, expMap.size() + 1, refMap.size() + 1);
	std::reverse(expMap.begin(), expMap.end());
	dp_fill_matrix(dpMatrix, expMap, refMap, minScoresSoFar);
	matchSequencesRev = dp_backtrack(dpMatrix, expMap.size() + 1, refMap.size() + 1);
	for (int ixMSR = 0; ixMSR < matchSequencesRev.size(); ixMSR++) matchSequencesRev[ixMSR].reversed = true;

	//add the reverse matches into the normal matches vector
	matchSequences.insert(matchSequences.end(), matchSequencesRev.begin(), matchSequencesRev.end());

	//restrict the list of sequences to topK		
	struct {
		bool operator()(Mapping a, Mapping b)
		{
			return a.score < b.score;
		}
	} mapLess;
	sort(matchSequences.begin(), matchSequences.end(), mapLess);
	if (params.topK < matchSequences.size()) matchSequences.erase(matchSequences.begin() + params.topK, matchSequences.end());

	//reverse the fragment back
	std::reverse(expMap.begin(), expMap.end());


	for (int i = 0; i < expMap.size() + 1; ++i) delete[] dpMatrix[i];
	delete[] dpMatrix;

	return matchSequences;
}


void map_segment(int from, int to, vector<ExpMap> &expMap, RefMaps &refMaps, Mappings* resultSet)
{
	for (int ixRM = from; ixRM <= to; ++ixRM)
	{
		int threshold = INDEX_NEIGHBORHOOD_THRESHOLD;

		Mappings mappings;
		vector<SCORE_TYPE> minScoresSoFar;
		for (int ix = 0; ix < params.topK; ++ix) minScoresSoFar.push_back(SUB_MAX);
		for (RefMaps::iterator refMap = refMaps.begin(); refMap != refMaps.end(); ++refMap)
		{
			Mappings aux_mapping = do_mapping(expMap[ixRM].reads, refMap->second, minScoresSoFar);
			for (Mappings::iterator itAux = aux_mapping.begin(); itAux != aux_mapping.end(); ++itAux)
			{
				itAux->chromosome = refMap->first;
				itAux->ComputeQuality();
			}
			mappings.insert(mappings.end(), aux_mapping.begin(), aux_mapping.end());
		}

		//keep top params.topK mappings
		struct {
			bool operator()(Mapping a, Mapping b)
			{
				return a.score < b.score;
			}
		} mapLess;
		sort(mappings.begin(), mappings.end(), mapLess);
		//sort(mappings.begin(), mappings.end(), [](Mapping & a, Mapping & b) -> bool	{return a.score < b.score; });
		mappings.erase(mappings.begin() + params.topK, mappings.end());
		resultSet[ixRM] = mappings;

		scoresCalculations.push_back(scoresComputed); scoresComputed = 0;

		mutexOM.lock();
		omProcessed++;
		ss << omProcessed << " "; logger.Log(Logger::LOGFILE, ss);
		mutexOM.unlock();
	}
}

void InitLogging()
{
	if (!params.logFileName.empty())	logger.InitChannel(Logger::LOGFILE, params.logFileName);
	if (params.outFileName != "") logger.InitChannel(Logger::RESFILE, params.outFileName);
}

void SmoothRefFragments(RefMaps &refMaps)
{
	for (RefMaps::iterator it = refMaps.begin(); it != refMaps.end(); it++)
	{
		for (int ixFrag = it->second.size() - 1; ixFrag >= 0; ixFrag--)
		{
			//we will proceed from the end and join every short fragment to its preceeding fragment
			if (it->second[ixFrag].length < params.smoothingThreshold && ixFrag > 0)
			{
				it->second[ixFrag - 1].length += it->second[ixFrag].length;
				it->second.erase(it->second.begin() + ixFrag);
			}
		}
		//first fragment can be too short so it needs to be joined with its successor 
		if (it->second[0].length < params.smoothingThreshold && it->second.size() > 1)
		{
			it->second[1].length += it->second[0].length;
			it->second.erase(it->second.begin());
		}
	}
}

void SmoothExpFragments(vector<ExpMap> &expMaps)
{
	for (vector<ExpMap>::iterator it = expMaps.begin(); it != expMaps.end(); it++)
	{
		for (int ixFrag = it->reads.size() - 1; ixFrag >= 0; ixFrag--)
		{
			//we will proceed from the end and join every short fragment to its preceeding fragment
			if (it->reads[ixFrag] < params.smoothingThreshold && ixFrag > 0)
			{
				it->reads[ixFrag - 1] += it->reads[ixFrag];
				it->reads.erase(it->reads.begin() + ixFrag);
			}
		}
		//first fragment can be too short so it needs to be joined with its successor 
		if (it->reads[0] < params.smoothingThreshold && it->reads.size() > 1)
		{
			it->reads[1] += it->reads[0];
			it->reads.erase(it->reads.begin());
		}
	}
}


void Parse(vector<ExpMap> &expMaps, RefMaps &refMaps)
{
	ss << "======= PARSE - START =======" << endl;
	logger.Log(Logger::LOGFILE, ss);
	clock_t begin_time = clock();
	//vector<ExpMap> expMap = parse_exp_map("../CASTEiJ_Alldata.maps", 1000);
	//vector<ExpMap> expMap = parse_exp_map("../ref.map.split", 100);
	expMaps = parse_exp_map(params.omFileName, 50);	
	refMaps = parse_ref_map(params.rmFileName);	
	SmoothExpFragments(expMaps);
	SmoothRefFragments(refMaps);
	ss << "ref. chromosomes: " << refMaps.size() << "\n"; logger.Log(Logger::LOGFILE, ss);
	int sum = 0;
	for (RefMaps::iterator it = refMaps.begin(); it != refMaps.end(); it++) sum += it->second.size();
	ss << "ref. maps total size: " << sum << "\n"; logger.Log(Logger::LOGFILE, ss);
	ss << "exp. maps length: " << expMaps.size() << "\n"; logger.Log(Logger::LOGFILE, ss);
	ss << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n"; logger.Log(Logger::LOGFILE, ss);
	ss << "======= PARSE - END =======" << endl; logger.Log(Logger::LOGFILE, ss);
}


Mappings* AlignOpticalMaps(vector<ExpMap> &expMap, RefMaps &refMaps)
{
	ss << "======= ALIGNING - START =======" << endl; logger.Log(Logger::LOGFILE, ss);
	clock_t begin_time = clock();
	//initialize a vector of results where the threads will store the mappings of individual reads
	//after all the exp maps will be mapped, the vector will be serialized
	Mappings *omMappings = new Mappings[expMap.size()];

	int cntThreads = thread::hardware_concurrency();
	if (params.cntThreads > 0) cntThreads = params.cntThreads;
	ss << "Using " << cntThreads << " threads" << endl; logger.Log(Logger::LOGFILE, ss);
	int batchSize = ceil(expMap.size() / (double)cntThreads);
	assert(batchSize > 0);
	vector<thread> threads;

	int ixFrom = 0, ixTo = expMap.size() - 1;
	if (params.ixOmStart > 0) ixFrom = params.ixOmStart;
	if (params.ixOmEnd >= 0) ixTo = min(params.ixOmEnd, (int)expMap.size() - 1);

	while (ixFrom <= ixTo)
	{
		int ixToAux = ixFrom + batchSize;
		if (ixToAux > ixTo) ixToAux = ixTo;
		threads.push_back(thread(map_segment, ixFrom, ixToAux, std::ref(expMap), std::ref(refMaps), omMappings));
		ixFrom = ixToAux + 1;
	}
	for (int i = 0; i < threads.size(); i++) threads[i].join();

	ss << endl << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl; logger.Log(Logger::LOGFILE, ss);
	ss << "======= ALIGNING - END =======" << endl; logger.Log(Logger::LOGFILE, ss);

	return omMappings;
}

void SerializeMappings(Mappings *omMappings, vector<ExpMap> &expMap, RefMaps &refMaps)
{
	ss << endl << "Outputting results..." << endl; logger.Log(Logger::LOGFILE, ss);
	ss << "ix;om_length;rm_length;length_diff;candidate_sections_length;score_calucations" << endl; logger.Log(Logger::STATSFILE, ss);
	ss << "#QX11 qaulity_score1;quality_score2;... (available in case of Bionano experimental maps)" << endl; logger.Log(Logger::RESFILE, ss);
	ss << "#QX12 signal_to_noise_ratio1;signal_to_noise_ratio2;... (available in case of Bionano experimental maps)" << endl; logger.Log(Logger::RESFILE, ss);
	ss << "#LEN_DIFF total_refmap_length - total_expmap_length" << endl; logger.Log(Logger::RESFILE, ss);
	ss << "#ALN aligned_ref_frags_len-aligned_exp_frags_len,#aligned_ref_frags:#aligned_exp_frags,aligned_ref_frags_len ..." << endl; logger.Log(Logger::RESFILE, ss);	
	ss << "#ALN_DETAIL aligned_ref_frags1:aligned_exp_frags1 aligned_ref_frags2:aligned_exp_frags2 ... (frags separated by comma)" << endl; logger.Log(Logger::RESFILE, ss);

	int cntIncorrectlyMapped = 0;;
	for (int ixOM = 0; ixOM < expMap.size(); ixOM++)
	{
		if (ixOM < params.ixOmStart || (ixOM > params.ixOmEnd && params.ixOmEnd != -1)) continue;
		ss << "EXP_OPTMAP_IX: " << ixOM << endl; logger.Log(Logger::RESFILE, ss);
		ss << "NAME: " << expMap[ixOM].name << endl; logger.Log(Logger::RESFILE, ss);	
		for (int ix = 0; ix < 2; ix++)
		{
			vector<float> q = (ix == 0 ? expMap[ixOM].qx11 : expMap[ixOM].qx12);
			ix == 0 ? ss << "QX11: " : ss << "QX12: ";
			for (int ixQ = 0; ixQ < q.size(); ixQ++)
			{
				if (ixQ > 0) ss << ";";
				ss << q[ixQ];
			}
			ss << endl;
			logger.Log(Logger::RESFILE, ss);
		}

		Mappings mappings = omMappings[ixOM];
		for (int ixMappings = 0; ixMappings < mappings.size(); ixMappings++)
		{
			string chr = mappings[ixMappings].chromosome;
			ss << "Mapping " << ixMappings << ": "; logger.Log(Logger::LOGFILE, ss);

			int omLength = 0, rmLength = 0;
			for (int ixAux = mappings[ixMappings].alignment.begin()->first + 1; ixAux <= (mappings[ixMappings].alignment.end() - 1)->first; ixAux++) omLength += expMap[ixOM].reads[ixAux - 1];
			for (int ixAux = mappings[ixMappings].alignment.begin()->second + 1; ixAux <= (mappings[ixMappings].alignment.end() - 1)->second; ixAux++) rmLength += refMaps[chr][ixAux - 1].length;

			string posOmFirst = "0";
			string posOmChrom = "";
			string posRmFirst = std::to_string(refMaps[chr][mappings[ixMappings].alignment[0].second].start);	//the DP matrix (and hence indeces) in the alignment is +1 shifted (init row and col)
			//with respect to the real sequences. But beginning of the alignment starts
			//-1 position before the "real" mapping starts -> no -1 shift needed
			string posRmChrom = refMaps[chr][mappings[ixMappings].alignment[0].second].chromosome;
			string posOmLast = std::to_string(expMap[0].length);
			auto al = mappings[ixMappings].alignment;
			string posRmLast = std::to_string(refMaps[chr][(mappings[ixMappings].alignment.end() - 1)->second - 1].start + refMaps[chr][(mappings[ixMappings].alignment.end() - 1)->second - 1].length);
			if (expMap[ixOM].debugInfo.size() > 0)
			{
				posOmChrom = expMap[ixOM].debugInfo[0];
				posOmFirst = expMap[ixOM].debugInfo[1];
				posOmLast = *(expMap[ixOM].debugInfo.end() - 1);
				int ixColon = posOmLast.find(':');
				if (ixColon != string::npos) posOmLast = posOmLast.substr(ixColon + 1);
			}

			ss << "REF_POS: " << posRmChrom << ":" << posRmFirst << "-" << posRmLast << endl; logger.Log(Logger::RESFILE, ss);
			ss << "QUALITY: " << mappings[ixMappings].quality << endl; logger.Log(Logger::RESFILE, ss);
			ss << "DP_SCORE: " << mappings[ixMappings].score << endl; logger.Log(Logger::RESFILE, ss);
			ss << "LEN_DIFF: " << rmLength - omLength << endl; logger.Log(Logger::RESFILE, ss);
			ss << "REVERSED: ";
			mappings[ixMappings].reversed ? ss << "1" : ss << "0";
			ss << endl; logger.Log(Logger::RESFILE, ss);

			std::ostringstream ssAln, ssAlnDetail;
			ssAln << "ALN: ";
			ssAlnDetail << "ALN_DETAIL: ";

			for (int ixAlignment = 1; ixAlignment < mappings[ixMappings].alignment.size(); ixAlignment++)
			{
				if (ixAlignment > 1){
					ssAln << " ";
					ssAlnDetail << " ";
				}
				int ixOMAux, ixRMAux;
				ostringstream ssRmPos, ssRmLengths, ssOmLengths;

				//get the previously matched pair + 1
				ixOMAux = mappings[ixMappings].alignment[ixAlignment - 1].first + 1;
				ixRMAux = mappings[ixMappings].alignment[ixAlignment - 1].second + 1;

				int sumOM = 0, sumRM = 0; //length between this and last matched position
				vector<int> auxRmIxs, auxOmIxs; //indeces of the aligned regions in RM and OM

				int cntOM = 0;
				while (ixOMAux <= mappings[ixMappings].alignment[ixAlignment].first)
				{
					int ixOmAuxReal = ixOMAux - 1;
					if (mappings[ixMappings].reversed) ixOmAuxReal = expMap[ixOM].reads.size() - ixOMAux;
					int length = expMap[ixOM].reads[ixOmAuxReal];
					sumOM += length; //ends of mapping are scored 0
					if (cntOM > 0) ssOmLengths << ",";
					ssOmLengths << length;
					ixOMAux++;
					cntOM++;
				}
				int cntRM = 0;
				while (ixRMAux <= mappings[ixMappings].alignment[ixAlignment].second)
				{
					int length = refMaps[chr][ixRMAux - 1].length;
					sumRM += length; //ends of mapping are scored 0
					ssRmPos << " " << refMaps[chr][ixRMAux - 1].chromosome << "_" << refMaps[chr][ixRMAux - 1].start;
					if (cntRM > 0) ssRmLengths << ",";
					ssRmLengths << length;
					ixRMAux++;
					cntRM++;
				}

				ss << ixOMAux + 1 << " - " << ixRMAux + 1 << " (" << sumOM << " - " << sumRM << ") "; logger.Log(Logger::LOGFILE, ss);
				ssAln << sumRM - sumOM << "," << cntRM << ":" << cntOM << "," << sumRM;// << sumOM << " ";
				ssAlnDetail << strings::trim(ssRmLengths.str()) << ":" << strings::trim(ssOmLengths.str());// << ":" << mappings[ixMappings].scores[ixAlignment - 1];
			}
			ss << endl; logger.Log(Logger::LOGFILE, ss);
			ssAln << endl; logger.Log(Logger::RESFILE, ssAln);
			ssAlnDetail << endl; logger.Log(Logger::RESFILE, ssAlnDetail);
		}
		ss << endl; logger.Log(Logger::RESFILE, ss);
		ss << "-----------------" << endl; logger.Log(Logger::LOGFILE, ss);
	}
	ss << "Incorreclty mapped fragments: " << cntIncorrectlyMapped << endl; logger.Log(Logger::LOGFILE, ss);

}

void ParseCmdLine(int argc, char** argv)
{
	try
	{
		TCLAP::CmdLine cmd("Optical mapping", ' ', "0.8");
		TCLAP::ValueArg<std::string> omFileNameArg("o", "expmap", "Experimental optical maps file (either plain text or gzipped)", true, "", "filename");
		TCLAP::ValueArg<std::string> rmFileNameArg("r", "refmap", "Reference map file (either plain text or gzipped)", true, "", "filename");
		ss.str(std::string());  ss << "Format of experiment file. Supported formats: [" << EXPERIMENT_FORMAT_TYPES << "]";
		TCLAP::ValueArg<std::string> formatArg("", "expformat", ss.str(), false, "opgen", "string"); ss.str(std::string());
		TCLAP::ValueArg<std::string> outFileNameArg("m", "outfile", "Output mapping file (if not present, the standard output will be used)", false, "", "filename");
		TCLAP::ValueArg<std::string> logFileNameArg("l", "logfile", "Log file", false, "", "filename");
		TCLAP::ValueArg<int> ixStartArg("b", "begin", "Index (zero-based) of the first fragment to map in the OM", false, 0, "int");
		TCLAP::ValueArg<int> ixEndArg("e", "end", "Index (zero-based) of the last fragment to map in the OM", false, -1, "int");
		TCLAP::ValueArg<int> cntThreadsArg("t", "threads", "Number of threads to use", false, 1, "int");
		TCLAP::ValueArg<int> topK("k", "topk", "Returns top K best mappings for each experimental map", false, 1, "int");
		TCLAP::ValueArg<string> chromosome("c", "chromosome", "Target chromosome (empty string = no restriction)", false, "", "string");
		TCLAP::ValueArg<int> dpwindowsize("", "miss-cnt", "Maximum number of missed or false restriction sites per aligned segment (maximal allowed value is 3).", false, 3, "int");
		TCLAP::ValueArg<float> smoothingThreshold("", "smooth-threshold", "Fragments shorther than this threshold will be merged with the neighbouring fragment", false, 1000, "int");

		ss.str(std::string());  ss << "Error model. Currently supported models: [" << ERROR_MODELS << "]";
		TCLAP::ValueArg<std::string> errorModelArg("", "errmodel", ss.str(), false, "valuev", "string"); ss.str(std::string());
		
		//TCLAP::ValueArg<int> omMissed("", "omissed", "Penalty for missing restriction site in an experimental optical map", false, 2000, "int");
		//TCLAP::ValueArg<int> rmMissed("", "rmmissed", "Penalty for missing restriction site in an refernce map", false, 2000, "int");
		
		TCLAP::ValueArg<float> sizingErrorStddev("", "read-error-stddev", "Fragment read error stddev. Size estimation error for a fragment  of length R is moddeled as N(0, est-error-stddev*R*R)", false, 0.02, "float");
		TCLAP::ValueArg<int> smallFragmentThreshold("", "small-fragment-threshold", "Sizing error stddev. Stddev for small fragments is constant ~ N(mean, est-error-stddev)", false, 4000, "int");
		TCLAP::ValueArg<float> digEff("", "cut-eff", "Cut (digestion) efficiency. Probabily of missing N restriction sites is (1 - cut-eff)^N", false, 0.8, "float");
		TCLAP::ValueArg<float> falseCutProb("", "false-cut-p", "Probability of false cut per base. Probability of N false cuts is modelled by Poisson distribution with mean = false-cut-p*segment_length", false, 0.00000001, "float");
		

		cmd.add(omFileNameArg);
		cmd.add(formatArg);
		cmd.add(rmFileNameArg);
		cmd.add(outFileNameArg);
		cmd.add(logFileNameArg);
		cmd.add(ixStartArg);
		cmd.add(ixEndArg);
		cmd.add(cntThreadsArg);
		cmd.add(topK);
		cmd.add(chromosome);
		//cmd.add(omMissed);
		//cmd.add(rmMissed);
		cmd.add(dpwindowsize);
		cmd.add(sizingErrorStddev);
		cmd.add(smallFragmentThreshold);
		cmd.add(digEff);
		cmd.add(falseCutProb);
		cmd.add(smoothingThreshold);
		cmd.add(errorModelArg);

		cmd.parse(argc, argv);

		params.omFileName = omFileNameArg.getValue();
		params.rmFileName = rmFileNameArg.getValue();
		params.omFormat = strings::lower(formatArg.getValue());
		params.outFileName = outFileNameArg.getValue();
		params.logFileName = logFileNameArg.getValue();
		params.ixOmStart = ixStartArg.getValue();
		params.ixOmEnd = ixEndArg.getValue();
		params.cntThreads = cntThreadsArg.getValue();
		params.topK = topK.getValue();
		params.chromosome = chromosome.getValue();
		params.errorModel = strings::lower(errorModelArg.getValue());
		//params.mapOmMissedPenalty = omMissed.getValue();
		//params.mapRmMissedPenalty = rmMissed.getValue();
		params.maxDpWindowSize = dpwindowsize.getValue() + 1;
		params.sizingErrorStddev = sizingErrorStddev.getValue();
		params.smallFragmentThreshold = smallFragmentThreshold.getValue();
		params.missRestrictionProb = 1 - digEff.getValue();
		params.noMissRestrictionProb = 1 - ((1 - pow(params.missRestrictionProb, params.maxDpWindowSize - 1)) / (1 - params.missRestrictionProb) - 1);
		params.falseCutProb = falseCutProb.getValue();
		params.smoothingThreshold = smoothingThreshold.getValue();

		if (params.maxDpWindowSize > MAX_OPT_MAP_WINDOW)
		{
			ss << "The maximum number of the miss-cnt parameter is " << MAX_OPT_MAP_WINDOW << "." << std::endl;
			error_exit(ss.str());
		}

		vector<string> allowedFormats = strings::split(EXPERIMENT_FORMAT_TYPES, ",");
		if (std::find(allowedFormats.begin(), allowedFormats.end(), params.omFormat) == allowedFormats.end())
		{
			ss << "The allowed experiment format values are " << EXPERIMENT_FORMAT_TYPES << "." << std::endl;
			error_exit(ss.str());
		}

		vector<string> allowedErrorModels = strings::split(ERROR_MODELS, ",");
		if (std::find(allowedErrorModels.begin(), allowedErrorModels.end(), params.errorModel) == allowedErrorModels.end())
		{
			ss << "The allowed error models are " << ERROR_MODELS << "." << std::endl;
			error_exit(ss.str());
		}
	}
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		ss << e.error() << " for arg " << e.argId() << std::endl;
		error_exit(ss.str());
	}
}

int main(int argc, char** argv)
{
	vector<ExpMap> expMap;
	RefMaps refMaps;	

	ParseCmdLine(argc, argv);
	stats::init_stats(params.falseCutProb);
	InitLogging();	
	Parse(expMap, refMaps);
	Mappings *omMatches = AlignOpticalMaps(expMap, refMaps);
	SerializeMappings(omMatches, expMap, refMaps);
	delete[] omMatches;

	return EXIT_SUCCESS;
}
