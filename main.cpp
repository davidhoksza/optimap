/*
* Copyright (C) 2016 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#define _SCL_SECURE_NO_WARNINGS

#include "common.h"
#include "constants.h"
#include "indexing.h"
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

vector<Fragment> parse_opt_map(string fileName, int topN = numeric_limits<int>::max())
{
	vector<Fragment> optMap;

	istream *ifs = open_map_file(fileName);

	string line;
	char buffer[20];
	vector<string> cleavageSites;

	string fragName;
	for (std::string line; getline(*ifs, line);) 	
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
			Fragment f;
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

			//boost::split(strs, boost::trim_copy(line), boost::is_any_of("\t"));
			//std::transform(strs.begin() + 2, strs.end(), std::back_inserter(numbers), [](const std::string& str) { return std::stof(str); });		
			if (cleavageSites.size() > 0)
			{
				assert(cleavageSites.size() - 2 == f.reads.size()); // debug info includes chromosome and first and last position -> +2
				f.debugInfo = cleavageSites;
				cleavageSites.clear();
			}
			optMap.push_back(f);
			if (optMap.size() >= topN) break;
		}
		else fragName = strings::trim(line);
	}

	delete ifs;
	return optMap;
}

void clean_dp_matrix(DpMatrixCell ** matrix, int height, int width )
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

void dp_fill_matrix(DpMatrixCell ** matrix, vector<int> &experiment, std::vector<RMRead> &reference, IndexRecord ir = IndexRecord())
{
	if (ir.start_position == -1) ir.start_position = 1;
	else ir.start_position++; //the start position is 0-based but here we have to account for the first zero coumn in the DP matrix
	if (ir.end_position == -1) ir.end_position = reference.size();
	else ir.end_position++; //the start position is 0-based but here we have to account for the first zero coumn in the DP matrix

	//first, let's initiliaze the first column with submax values which ensures 
	//that the resulting mapping will capture the whole fragemnt
	//we don't use max values because that might cause overflow
	//since we add to these values
	for (int ixRow = 1; ixRow <= experiment.size(); matrix[ixRow++][ir.start_position - 1].value = SUB_MAX);

	SCORE_TYPE minMappingValue = numeric_limits<SCORE_TYPE>::max();

	for (int ixRow = 1; ixRow < experiment.size() + 1; ++ixRow)
	{
		bool isLastRow = ixRow == experiment.size() ? true : false;

		//we need to find the first column in the last row where it makes sense to search for mins
		//for example the first column does not make sense since that would mean the whole OM fragment was
		//aligned with one read in the reference map
		//int ixColResultFrom = ceil(fragment.size() / (float) DP_WINDOW_SIZE);

		for (int ixCol = ir.start_position; ixCol <= ir.end_position; ++ixCol)
		{
			//matrix[ixRow][ixCol].value = matrix[ixRow - 1][ixCol - 1].value + abs(fragment[ixRow - 1] - reference[ixCol - 1]);
			DpMatrixCell minCell;
			minCell.value = numeric_limits<SCORE_TYPE>::max();
			minCell.sourceColumn = ixCol - 1; //This value is used when no good match is found
			minCell.sourceRow = ixRow - 1; //This value is used when no good match is found

			//First and last fragments are scored 0
			if (ixRow == 1 || isLastRow)
			{
				isLastRow ? minCell.value = matrix[ixRow - 1][ixCol - 1].value + stats::transform_prob(1) : minCell.value = stats::transform_prob(1);				
				matrix[ixRow][ixCol] = minCell;
				continue;
			}
			
			
			int rowValue = 0;
			for (int ixWindowRow = 1; ixWindowRow <= params.mapDpWindowSize; ++ixWindowRow)
			{
				if (ixRow - ixWindowRow < 0) break; //should I touch position out of the array
				int ixExp = ixRow - ixWindowRow;
				//ends of fragments can be aligned with zero score since
				//the molecules forming fragments were not created with a restriction enzyme
				//if (ixFrag > 0 && ixFrag < fragment.size()-1) 
				rowValue += experiment[ixExp];
				int colValue = 0;
				for (int ixWindowCol = 1; ixWindowCol <= params.mapDpWindowSize; ++ixWindowCol)
				{
					if (ixCol - ixWindowCol < ir.start_position - 1) break; //should I touch position out of the candidate window
					colValue += reference[ixCol - ixWindowCol].length; //since the maps and dp table are shifted by 1, this returns in the first iteration the inspected position ([ixRow,ixCol])
					//float score = matrix[ixRow - ixWindowRow][ixCol - ixWindowCol].value + pow(rowValue - colValue, 2)/(colValue*1.05);

					//if (ixFrag > 0 && ixFrag < fragment.size() - 1) rowValue = colValue; //if the aligned segment includes first fragment from the experimental map,
																						//we have to treat exp and reference segments as having equal lengths since we do not know
																						//the real length of the first experimental fragment. 
																						//(molecules were not created by the same restriction enzyme as the interior cleavage sites)
										
					SCORE_TYPE score = matrix[ixRow - ixWindowRow][ixCol - ixWindowCol].value;
					if (score >= SUB_MAX) continue;
					float stddev = params.sizingErrorStddev * (colValue > params.smallFragmentThreshold ? colValue : params.smallFragmentThreshold);

					float x = stats::pdf_gaussian((rowValue - colValue) / stddev, 0, 1);
					score += stats::transform_prob(x);
					if (score >= SUB_MAX) continue;

					//penalty computation
					//score += (ixWindowRow - 1) * params.mapOmMissedPenalty + (ixWindowCol - 1)* params.mapRmMissedPenalty;
					SCORE_TYPE auxP;
					if (params.mapDpWindowSize == 1) auxP = 1;
					else 
					{
						if (ixWindowCol > 1) auxP = pow(params.missRestrictionProb, ixWindowCol - 1);
						else auxP = params.noMissRestrictionProb;
					}
					score += stats::transform_prob(auxP);


					x = stats::pdf_poisson_full(ixWindowRow - 1, rowValue * params.falseCutProb);
					// The Poisson PDF is precomputed with rowValue * params.falseCutProb. But since falseCutProb is constant,
					// we can provide only the rowValue and use it as index to the array with precomputed values.
					//x = stats::pdf_poisson(ixWindowRow - 1, rowValue);
					score += stats::transform_prob(x);

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
}


Mappings dp_backtrack(DpMatrixCell **matrix, int height, int width, IndexRecord ir = IndexRecord())
{
	Mappings mappings;

	vector<pair<SCORE_TYPE, int>> minPositions; //top K min values and positions in increasing order
	if (ir.start_position == -1) ir.start_position = 1;
	else ir.start_position++;
	if (ir.end_position == -1) ir.end_position = width - 1;
	else ir.end_position++;

	for (int ixM = ir.start_position; ixM <= ir.end_position; ++ixM)
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

	Mappings matchSequences, matchSequencesRev;
	if (candidates.size() == 0)
	{
		dp_fill_matrix(dpMatrix, fragment, refMap);
		matchSequences = dp_backtrack(dpMatrix, fragment.size() + 1, refMap.size() + 1);

		//repeat the process with reversed fragment
		clean_dp_matrix(dpMatrix, fragment.size() + 1, refMap.size() + 1);
		std::reverse(fragment.begin(), fragment.end());
		dp_fill_matrix(dpMatrix, fragment, refMap);
		matchSequencesRev = dp_backtrack(dpMatrix, fragment.size() + 1, refMap.size() + 1);
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
		std::reverse(fragment.begin(), fragment.end());
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
		for (RefMaps::iterator refMap = refMaps.begin(); refMap != refMaps.end(); ++refMap)
		{	
			//resultSet[ixOM] = do_mapping(optMap[ixOM].reads, refMap.second, candidates);
			Mappings aux_mapping = do_mapping(optMap[ixOM].reads, refMap->second, candidates);
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
		resultSet[ixOM] = mappings;

		//verify_candidates(optMap[ixOM], refMap, candidates);		
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
	if (params.outFileName!= "") logger.InitChannel(Logger::RESFILE, params.outFileName); 
	//logger.InitChannel(Logger::STATSFILE, "../stats.csv");
}

void ParseCmdLine(int argc, char** argv)
{
	try
	{
		TCLAP::CmdLine cmd("Optical mapping", ' ', "0.8");
		TCLAP::ValueArg<std::string> omFileNameArg("o", "optmap", "Optical maps file (either plain text or gzipped)", true, "", "filename");
		TCLAP::ValueArg<std::string> rmFileNameArg("r", "refmap", "Reference map file (either plain text or gzipped)", true, "", "filename");
		TCLAP::ValueArg<std::string> outFileNameArg("m", "outfile", "Output mapping file (if not present, the standard output will be used)", false, "", "filename");
		TCLAP::ValueArg<std::string> logFileNameArg("l", "logfile", "Log file", false, "", "filename");
		TCLAP::ValueArg<int> ixStartArg("b", "begin", "Index (zero-based) of the first fragment to map in the OM", false, 0, "int");
		TCLAP::ValueArg<int> ixEndArg("e", "end", "Index (zero-based) of the last fragment to map in the OM", false, -1, "int");
		TCLAP::ValueArg<int> cntThreadsArg("t", "threads", "Number of threads to use", false, 1, "int");
		TCLAP::ValueArg<int> topK("k", "topk", "Finds top K best mappings for each optical map", false, 1, "int");
		TCLAP::ValueArg<string> chromosome("c", "chromosome", "Target chromosome (empty string = no restriction)", false, "", "string");
		//TCLAP::ValueArg<int> omMissed("", "omissed", "Penalty for missing restriction site in an experimental optical map", false, 2000, "int");
		//TCLAP::ValueArg<int> rmMissed("", "rmmissed", "Penalty for missing restriction site in an refernce map", false, 2000, "int");
		TCLAP::ValueArg<int> dpwindowsize("", "miss-cnt", "Maximum number of missed or false restriction sites per aligned segment (maximal allowed value is 10).", false, 3, "int");		
		TCLAP::ValueArg<float> sizingErrorStddev("", "read-error-stddev", "Fragment read error stddev. Size estimation error for a fragment \
																	of length R is moddeled as N(0, est-error-stddev*R*R)", true, 0.02, "float");
		TCLAP::ValueArg<int> smallFragmentThreshold("", "small-fragment-threshold", "Sizing error stddev. \
																					Stddev for small fragments is constant ~ N(mean, est-error-stddev)", false, 4000, "int");
		TCLAP::ValueArg<float> digEff("", "cut-eff", "Cut (digestion) efficiency. Probabily of missing N restriction sites is (1 - cut-eff)^N", false, 0.8, "float");
		TCLAP::ValueArg<float> falseCutProb("", "false-cut-p", "Probability of false cut per base. Probability of N false cuts is modelled by \
															   Poisson distribution with mean = false-cut-p*segment_length", false, 0.00000001, "float");

		cmd.add(omFileNameArg);
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

		cmd.parse(argc, argv);

		params.omFileName = omFileNameArg.getValue();
		params.rmFileName = rmFileNameArg.getValue();
		params.outFileName = outFileNameArg.getValue();
		params.logFileName = logFileNameArg.getValue();
		params.ixOmStart = ixStartArg.getValue();
		params.ixOmEnd = ixEndArg.getValue();
		params.cntThreads = cntThreadsArg.getValue();
		params.topK = topK.getValue();
		params.chromosome = chromosome.getValue();
		//params.mapOmMissedPenalty = omMissed.getValue();
		//params.mapRmMissedPenalty = rmMissed.getValue();
		params.mapDpWindowSize = dpwindowsize.getValue();
		params.sizingErrorStddev = sizingErrorStddev.getValue();
		params.smallFragmentThreshold = smallFragmentThreshold.getValue();
		params.missRestrictionProb = 1 - digEff.getValue();
		params.noMissRestrictionProb =  1 - ((1 - pow(params.missRestrictionProb, params.mapDpWindowSize - 1)) / (1 - params.missRestrictionProb) - 1);
		params.falseCutProb = falseCutProb.getValue();

		if (params.mapDpWindowSize > MAX_OPT_MAP_WINDOW) error_exit("The maximum number of the miss-cnt parameter is 10.");
	}
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		ss << e.error() << " for arg " << e.argId() << std::endl;
		error_exit(ss.str());
	}
}

void Parse(vector<Fragment> &optMap, RefMaps &refMaps)
{
	ss << "======= PARSE - START =======" << endl;
	logger.Log(Logger::LOGFILE, ss);
	clock_t begin_time = clock();
	//vector<Fragment> optMap = parse_opt_map("../CASTEiJ_Alldata.maps", 1000);
	//vector<Fragment> optMap = parse_opt_map("../ref.map.split", 100);
	optMap = parse_opt_map(params.omFileName);
	refMaps = parse_ref_map(params.rmFileName); //vector<RMRead> refMap = parse_ref_map("../ref.map"/*, 100000*/);
	ss << "ref. chromosomes: " << refMaps.size() << "\n"; logger.Log(Logger::LOGFILE, ss);
	int sum = 0;
	for (RefMaps::iterator it = refMaps.begin(); it != refMaps.end(); it++) sum += it->second.size();
	ss << "ref. maps total size: " << sum << "\n"; logger.Log(Logger::LOGFILE, ss);
	ss << "opt. map length: " << optMap.size() << "\n"; logger.Log(Logger::LOGFILE, ss);
	ss << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n"; logger.Log(Logger::LOGFILE, ss);
	ss << "======= PARSE - END =======" << endl; logger.Log(Logger::LOGFILE, ss);
}

void InitializeIndex(map<int, vector<IndexRecord> > &index, RefMaps &refMaps)
{
	ss << "======= INDEX INITIALIZATION - START =======" << endl; logger.Log(Logger::LOGFILE, ss);
	clock_t begin_time = clock();
	index = init_index(refMaps);
	ss << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl; logger.Log(Logger::LOGFILE, ss);
	ss << "======= INDEX INITIALIZATION  - END =======" << endl; logger.Log(Logger::LOGFILE, ss);

}

Mappings* AlignOpticalMaps(vector<Fragment> &optMap, RefMaps &refMaps, map<int, vector<IndexRecord> > &index)
{
	ss << "======= ALIGNING - START =======" << endl; logger.Log(Logger::LOGFILE, ss);
	clock_t begin_time = clock();
	//initialize a vector of results where the threads will store the mappings of individual reads
	//after all the opt maps will be mapped, the vector will be serialized
	Mappings *omMappings = new Mappings[optMap.size()];

	int cntThreads = thread::hardware_concurrency();
	if (params.cntThreads > 0) cntThreads = params.cntThreads;
	ss << "Using " << cntThreads << " threads" << endl; logger.Log(Logger::LOGFILE, ss);
	int batchSize = ceil(optMap.size() / (double)cntThreads);
	assert(batchSize > 0);
	vector<thread> threads;

	int ixFrom = 0, ixTo = optMap.size() - 1;
	if (params.ixOmStart > 0) ixFrom = params.ixOmStart;
	if (params.ixOmEnd >= 0) ixTo = min(params.ixOmEnd, (int)optMap.size() - 1);

	while (ixFrom <= ixTo)
	{
		int ixToAux = ixFrom + batchSize;
		if (ixToAux > ixTo) ixToAux = ixTo;
		threads.push_back(thread(map_segment, ixFrom, ixToAux, std::ref(optMap), std::ref(refMaps), omMappings, std::ref(index)));
		ixFrom = ixToAux + 1;
	}
	for (int i = 0; i < threads.size(); i++) threads[i].join();

	ss << endl << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl; logger.Log(Logger::LOGFILE, ss);
	ss << "======= ALIGNING - END =======" << endl; logger.Log(Logger::LOGFILE, ss);

	return omMappings;
}

void SerializeMappings(Mappings *omMappings, vector<Fragment> &optMap, RefMaps &refMaps)
{
	ss << endl << "Outputting results..." << endl; logger.Log(Logger::LOGFILE, ss);
	ss << "ix;om_length;rm_length;length_diff;candidate_sections_length;score_calucations" << endl; logger.Log(Logger::STATSFILE, ss);
	//int ixCSL = 0;

	/*ss << "#om id;om indeces;rm positions;om lengths; rm lengths" << endl; logger.Log(Logger::RESFILE, ss);*/

	ss << "#LEN_DIFF total_refmap_length - total_expmap_length" << endl; logger.Log(Logger::RESFILE, ss);
	ss << "#ALN aligned_ref_frags_len-aligned_exp_frags_len,#aligned_ref_frags:#aligned_exp_frags,aligned_ref_frags_len ..." << endl; logger.Log(Logger::RESFILE, ss);
	ss << "#ALN_DETAIL aligned_ref_frags1:aligned_exp_frags1 aligned_ref_frags2:aligned_exp_frags2 ... (frags separated by comma)" << endl; logger.Log(Logger::RESFILE, ss);

	int cntIncorrectlyMapped = 0;;
	for (int ixOM = 0; ixOM < optMap.size(); ixOM++)
	{
		if (ixOM < params.ixOmStart || (ixOM > params.ixOmEnd && params.ixOmEnd != -1)) continue;
		ss << "EXP_OPTMAP_IX: " << ixOM << endl; logger.Log(Logger::RESFILE, ss);
		ss << "NAME: " << optMap[ixOM].name << endl; logger.Log(Logger::RESFILE, ss);
		Mappings mappings = omMappings[ixOM];
		for (int ixMappings = 0; ixMappings < mappings.size(); ixMappings++)
		{
			string chr = mappings[ixMappings].chromosome;
			ss << "Mapping " << ixMappings << ": "; logger.Log(Logger::LOGFILE, ss);

			int omLength = 0, rmLength = 0;
			for (int ixAux = mappings[ixMappings].alignment.begin()->first + 1; ixAux <= (mappings[ixMappings].alignment.end() - 1)->first; ixAux++) omLength += optMap[ixOM].reads[ixAux - 1];
			for (int ixAux = mappings[ixMappings].alignment.begin()->second + 1; ixAux <= (mappings[ixMappings].alignment.end() - 1)->second; ixAux++) rmLength += refMaps[chr][ixAux - 1].length;
			//int aux = (candidateSectionLengths.size() > ixCSL ? candidateSectionLengths[ixCSL++] : -1);
			//ss << ixOM << ";" << omLength << ";" << rmLength << ";" << abs(omLength - rmLength) << ";" << aux << ";" << scoresCalculations[ixOM] << endl; logger.Log(Logger::STATSFILE, ss);

			//ss << "# " << ixMappings + 1 << ". mapping with score " << mappings[ixMappings].score << " and real length difference "  << abs(omLength - rmLength) << endl; logger.Log(Logger::RESFILE, ss);
			//int posOmFirst = optMap[ixOM].reads[mappings[ixMappings].alignment[0].first];
			string posOmFirst = "0";
			string posOmChrom = "";
			string posRmFirst = std::to_string(refMaps[chr][mappings[ixMappings].alignment[0].second].start);	//the DP matrix (and hence indeces) in the alignment is +1 shifted (init row and col)
			//with respect to the real sequences. But beginning of the alignment starts
			//-1 position before the "real" mapping starts -> no -1 shift needed
			string posRmChrom = refMaps[chr][mappings[ixMappings].alignment[0].second].chromosome;
			//int posOmLast = optMap[ixOM].reads[(mappings[ixMappings].alignment.end() - 1)->first]; 
			string posOmLast = std::to_string(optMap[0].length);
			auto al = mappings[ixMappings].alignment;
			string posRmLast = std::to_string(refMaps[chr][(mappings[ixMappings].alignment.end() - 1)->second - 1].start + refMaps[chr][(mappings[ixMappings].alignment.end() - 1)->second - 1].length);
			if (optMap[ixOM].debugInfo.size() > 0)
			{
				posOmChrom = optMap[ixOM].debugInfo[0];
				posOmFirst = optMap[ixOM].debugInfo[1];
				posOmLast = *(optMap[ixOM].debugInfo.end() - 1);
				int ixColon = posOmLast.find(':');
				if (ixColon != string::npos) posOmLast = posOmLast.substr(ixColon + 1);
			}

			/*ss << "# " << ixMappings + 1 << "\tscore: " << mappings[ixMappings].score << "\tlength-diff: " << abs(omLength - rmLength)
			<< "\tmap: " << posOmChrom << ":" << posOmFirst << "-" << posOmLast << "->" << posRmChrom << ":" << posRmFirst << "-" << posRmLast << endl; logger.Log(Logger::RESFILE, ss);*/

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
				ostringstream ssRmPos, ssRmLengths, ssOmIxs, ssOmLengths;

				//get the previously matched pair + 1
				ixOMAux = mappings[ixMappings].alignment[ixAlignment - 1].first + 1;
				ixRMAux = mappings[ixMappings].alignment[ixAlignment - 1].second + 1;

				int sumOM = 0, sumRM = 0; //length between this and last matched position
				vector<int> auxRmIxs, auxOmIxs; //indeces of the aligned regions in RM and OM

				int cntOM = 0;
				while (ixOMAux <= mappings[ixMappings].alignment[ixAlignment].first)
				{
					int length = optMap[ixOM].reads[ixOMAux - 1];
					//if (cntOM > 0 && ixOMAux < mappings[ixMappings].alignment[ixAlignment].first) 
					sumOM += length; //ends of mapping are scored 0
					ssOmIxs << " " << ixOMAux - 1;
					if (cntOM > 0) ssOmLengths << ",";
					ssOmLengths << length;
					ixOMAux++;
					cntOM++;
				}
				int cntRM = 0;
				while (ixRMAux <= mappings[ixMappings].alignment[ixAlignment].second)
				{
					int length = refMaps[chr][ixRMAux - 1].length;
					//if (cntRM > 0 && ixRMAux < mappings[ixMappings].alignment[ixAlignment].second) 
					sumRM += length; //ends of mapping are scored 0
					ssRmPos << " " << refMaps[chr][ixRMAux - 1].chromosome << "_" << refMaps[chr][ixRMAux - 1].start;
					if (cntRM > 0) ssRmLengths << ",";
					ssRmLengths << length;
					ixRMAux++;
					cntRM++;
				}
				/*ss << "\t" << ixOM << ";" << strings::trim(ssOmIxs.str()) << ";" << strings::trim(ssRmPos.str()) << ";"
				<< strings::trim(ssOmLengths.str()) << ";" << strings::trim(ssRmLengths.str()) << endl; logger.Log(Logger::RESFILE, ss);*/

				ss << ixOMAux + 1 << " - " << ixRMAux + 1 << " (" << sumOM << " - " << sumRM << ") "; logger.Log(Logger::LOGFILE, ss);
				ssAln << sumRM - sumOM << "," << cntRM << ":" << cntOM << "," << sumRM;// << sumOM << " ";
				ssAlnDetail << strings::trim(ssRmLengths.str()) << ":" << strings::trim(ssOmLengths.str());
			}
			ss << endl; logger.Log(Logger::LOGFILE, ss);
			ssAln << endl; logger.Log(Logger::RESFILE, ssAln);
			ssAlnDetail << endl; logger.Log(Logger::RESFILE, ssAlnDetail);

			//For debugging purposes we want to check whether the OM comes from the same place in RM (debugInfo in OM = position in RM)	
			//Beggining of the OM segement is stored in debugInfo[2] and end in the last element of debugInfo, matching segements in refmap is stored in matches.matches
			//int omStart, omEnd, rmStart, rmEnd;
			//omStart = stoi(optMap[ixOM].debugInfo[1]);
			//omEnd = stoi(optMap[ixOM].debugInfo[optMap[ixOM].debugInfo.size() - 1].substr(2));
			//rmStart = refMap.map[matches.matches[ixMatches][matches.matches[ixMatches].size() - 1].second].start;
			//rmEnd = refMap.map[matches.matches[ixMatches][0].second - 1].start + refMap.map[matches.matches[ixMatches][0].second - 1].length + 6; //6 for the cleavage
			//if ((omStart != rmStart && omStart != rmStart - 1 && omStart != rmStart + 1) || (omEnd != rmEnd && omEnd != rmEnd - 1 && omEnd != rmEnd + 1))
			//{
			//	cntIncorrectlyMapped++;
			//	ss << "OM info (reads): ";
			//	if (optMap[ixOM].debugInfo.size() > 0)
			//	{
			//		for (int ixDI = 0; ixDI < optMap[ixOM].debugInfo.size(); ixDI++) ss << optMap[ixOM].debugInfo[ixDI] << " ";
			//		ss << endl; logger.Log(Logger::LOGFILE, ss);
			//	}
			//	ss << "RM info (first-last position): " << rmStart << " - " << rmEnd;
			//	ss << endl << endl; logger.Log(Logger::LOGFILE, ss);
			//}
		}		
		ss << endl; logger.Log(Logger::RESFILE, ss);
		ss << "-----------------" << endl; logger.Log(Logger::LOGFILE, ss);
	}
	ss << "Incorreclty mapped fragments: " << cntIncorrectlyMapped << endl; logger.Log(Logger::LOGFILE, ss);

}

int main(int argc, char** argv)
{
	vector<Fragment> optMap;
	RefMaps refMaps;
	map<int, vector<IndexRecord> > index;

	ParseCmdLine(argc, argv);
	stats::init_stats(params.falseCutProb);
	InitLogging();
	Parse(optMap, refMaps);
	//InitializeIndex(index, refMap);	
	
	Mappings *omMatches = AlignOpticalMaps(optMap, refMaps, index);
	SerializeMappings(omMatches, optMap, refMaps);
	delete[] omMatches;

	return EXIT_SUCCESS;
}
