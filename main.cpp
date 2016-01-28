/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#define _SCL_SECURE_NO_WARNINGS

#include "constants.h"
#include "indexing.h"
#include "types.h"
#include "logger.h"
#include "string_functions.h"

#include "tclap/CmdLine.h"

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

RefMaps parse_ref_map(string fileName)
{
	RefMaps refMaps;
	

	ifstream ifs(fileName);

	if (!ifs.is_open())
	{
		cout << "ERROR: Reference map file " + fileName + " could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	for (string line; getline(ifs, line);) {
		vector<string> strs = strings::split(line, "\t");

		RMRead auxRMRead;
		auxRMRead.chromosome = strings::trim(strs[0]);		
		if (params.chromosome != "" && strings::upper(auxRMRead.chromosome) != strings::upper(params.chromosome)) continue;
		auxRMRead.start = stof(strs[1]);
		auxRMRead.length = stof(strs[3]) * 1000;
		if (refMaps.count(auxRMRead.chromosome) == 0) refMaps[auxRMRead.chromosome] = vector<RMRead>();
		refMaps[auxRMRead.chromosome].push_back(auxRMRead);
	}

	return refMaps;
}

vector<Fragment> parse_opt_map(string fileName, int topN = numeric_limits<int>::max())
{
	vector<Fragment> optMap;

	ifstream ifs(fileName);

	if (!ifs.is_open())
	{
		cout << "ERROR: Optical map file " + fileName + " could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	char buffer[20];
	vector<string> cleavageSites;
	for (std::string line; getline(ifs, line);) 
	{		
		if (line.find_first_of("debug") != string::npos)
		{
			stringstream ss(line); // Insert the string into a stream
			string strBuffer;
			ss >> strBuffer;
			while (ss >> strBuffer) cleavageSites.push_back(strBuffer);
		}
		if (line.find_first_of("KpnI") != string::npos)
		{
			Fragment f;

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
				buffer[pos+1] = 0;
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
	}

	return optMap;
}

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
	for (int ixRow = 1; ixRow <= fragment.size(); matrix[ixRow++][ir.start_position-1].value = numeric_limits<int>::max() / 1000);

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
					if (ixCol - ixWindowCol < ir.start_position-1) break; //should I touch position out of the candidate window
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
	if (ir.end_position == -1) ir.end_position = width-1;
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
		for (auto refMap : refMaps)
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

void InitLogging()
{
	if (!params.logFileName.empty())	logger.InitChannel(Logger::LOGFILE, params.logFileName);
	logger.InitChannel(Logger::RESFILE, params.outFileName); //logger.InitChannel(Logger::RESFILE, "../mapping.out");
	logger.InitChannel(Logger::STATSFILE, "../stats.csv");
}

void ParseCmdLine(int argc, char** argv)
{
	try
	{
		TCLAP::CmdLine cmd("Optical mapping", ' ', "0.8");
		TCLAP::ValueArg<std::string> omFileNameArg("o", "optmap", "Optical maps file", true, "", "filename");
		TCLAP::ValueArg<std::string> rmFileNameArg("r", "refmap", "Reference map file", true, "", "filename");
		TCLAP::ValueArg<std::string> outFileNameArg("m", "outfile", "Output mapping file", true, "", "filename");
		TCLAP::ValueArg<std::string> logFileNameArg("l", "logfile", "Log file", false, "", "filename");
		TCLAP::ValueArg<int> ixStartArg("b", "begin", "Index (zero-based) of the first fragment to map in the OM", false, 0, "int");
		TCLAP::ValueArg<int> ixEndArg("e", "end", "Index (zero-based) of the last fragment to map in the OM", false, -1, "int");
		TCLAP::ValueArg<int> cntThreadsArg("t", "threads", "Number of threads to use", false, -1, "int");
		TCLAP::ValueArg<int> topK("k", "topk", "Finds top K best mappings for each optical map", false, 1, "int");
		TCLAP::ValueArg<string> chromosome("c", "chromosome", "Target chromosome (empty string = no restriction)", false, "", "string");
		TCLAP::ValueArg<int> omMissed("", "omissed", "Penalty for missing restriction site in an experimental optical map", false, 1000, "int");
		TCLAP::ValueArg<int> rmMissed("", "rmmissed", "Penalty for missing restriction site in an refernce map", false, 1000, "int");
		TCLAP::ValueArg<int> dpwindowsize("", "dpwindowsize", "Size of the dynamic programming window, i.e. how many restriction sites can be missed", false, 2, "int");
		
		cmd.add(omFileNameArg);
		cmd.add(rmFileNameArg);
		cmd.add(outFileNameArg);
		cmd.add(logFileNameArg);
		cmd.add(ixStartArg);
		cmd.add(ixEndArg);
		cmd.add(cntThreadsArg);
		cmd.add(topK);
		cmd.add(chromosome);
		cmd.add(omMissed);
		cmd.add(rmMissed);
		cmd.add(dpwindowsize);
		
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
		params.mapOmMissedPenalty = omMissed.getValue();
		params.mapRmMissedPenalty = rmMissed.getValue();
		params.mapDpWindowSize = dpwindowsize.getValue();

	}
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;

		exit(EXIT_FAILURE);
	}
}

void Parse(vector<Fragment> &optMap, RefMaps &refMaps)
{
	ss << "======= PARSE - START =======" << endl;
	logger.Log(Logger::STDOUT, ss);
	clock_t begin_time = clock();
	//vector<Fragment> optMap = parse_opt_map("../CASTEiJ_Alldata.maps", 1000);
	//vector<Fragment> optMap = parse_opt_map("../ref.map.split", 100);
	optMap = parse_opt_map(params.omFileName);
	refMaps = parse_ref_map(params.rmFileName); //vector<RMRead> refMap = parse_ref_map("../ref.map"/*, 100000*/);
	ss << "ref. chromosomes: " << refMaps.size() << "\n"; logger.Log(Logger::STDOUT, ss);
	int sum = 0;
	for (auto m : refMaps) sum += m.second.size();
	ss << "ref. maps total size: " << sum << "\n"; logger.Log(Logger::STDOUT, ss);
	ss << "opt. map length: " << optMap.size() << "\n"; logger.Log(Logger::STDOUT, ss);
	ss << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n"; logger.Log(Logger::STDOUT, ss);
	ss << "======= PARSE - END =======" << endl; logger.Log(Logger::STDOUT, ss);
}

void InitializeIndex(map<int, vector<IndexRecord> > &index, RefMaps &refMaps)
{
	ss << "======= INDEX INITIALIZATION - START =======" << endl; logger.Log(Logger::STDOUT, ss);
	clock_t begin_time = clock();
	index = init_index(refMaps);
	ss << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl; logger.Log(Logger::STDOUT, ss);
	ss << "======= INDEX INITIALIZATION  - END =======" << endl; logger.Log(Logger::STDOUT, ss);

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

void SerializeMappings(Mappings *omMappings, vector<Fragment> &optMap, RefMaps &refMaps)
{
	ss << endl << "Outputting results..." << endl; logger.Log(Logger::STDOUT, ss);
	ss << "ix;om_length;rm_length;length_diff;candidate_sections_length;score_calucations" << endl; logger.Log(Logger::STATSFILE, ss);
	int ixCSL = 0;

	/*ss << "#om id;om indeces;rm positions;om lengths; rm lengths" << endl; logger.Log(Logger::RESFILE, ss);*/

	int cntIncorrectlyMapped = 0;;
	for (int ixOM = 0; ixOM < optMap.size(); ixOM++)
	{
		if (ixOM < params.ixOmStart || (ixOM > params.ixOmEnd && params.ixOmEnd != -1)) continue;
		ss << "EXP_OPTMAP_IX: " << ixOM << endl; logger.Log(Logger::RESFILE, ss);
		Mappings mappings = omMappings[ixOM];
		for (int ixMappings = 0; ixMappings < mappings.size(); ixMappings++)
		{
			string chr = mappings[ixMappings].chromosome;
			ss << "Mapping " << ixMappings << ": "; logger.Log(Logger::LOGFILE, ss);

			int omLength = 0, rmLength = 0;
			for (int ixAux = mappings[ixMappings].alignment.begin()->first + 1; ixAux <= (mappings[ixMappings].alignment.end() - 1)->first; ixAux++) omLength += optMap[ixOM].reads[ixAux - 1];
			for (int ixAux = mappings[ixMappings].alignment.begin()->second + 1; ixAux <= (mappings[ixMappings].alignment.end() - 1)->second; ixAux++) rmLength += refMaps[chr][ixAux - 1].length;
			int aux = (candidateSectionLengths.size() > ixCSL ? candidateSectionLengths[ixCSL++] : -1);
			ss << ixOM << ";" << omLength << ";" << rmLength << ";" << abs(omLength - rmLength) << ";" << aux << ";" << scoresCalculations[ixOM] << endl; logger.Log(Logger::STATSFILE, ss);

			//ss << "# " << ixMappings + 1 << ". mapping with score " << mappings[ixMappings].score << " and real length difference "  << abs(omLength - rmLength) << endl; logger.Log(Logger::RESFILE, ss);
			//int posOmFirst = optMap[ixOM].reads[mappings[ixMappings].alignment[0].first];
			string posOmFirst = "0";
			string posOmChrom = "";			
			string posRmFirst = std::to_string(refMaps[chr][mappings[ixMappings].alignment[0].second].start);
			string posRmChrom = refMaps[chr][mappings[ixMappings].alignment[0].second].chromosome;
			//int posOmLast = optMap[ixOM].reads[(mappings[ixMappings].alignment.end() - 1)->first]; 
			string posOmLast = std::to_string(optMap[0].length);
			string posRmLast = std::to_string(refMaps[chr][(mappings[ixMappings].alignment.end() - 1)->second].start);
			if (optMap[ixOM].debugInfo.size() > 0)
			{
				posOmChrom = optMap[ixOM].debugInfo[0];
				posOmFirst =  optMap[ixOM].debugInfo[1];
				posOmLast = *(optMap[ixOM].debugInfo.end()-1);
				int ixColon = posOmLast.find(':');
				if (ixColon != string::npos) posOmLast = posOmLast.substr(ixColon+1);
			}

			/*ss << "# " << ixMappings + 1 << "\tscore: " << mappings[ixMappings].score << "\tlength-diff: " << abs(omLength - rmLength) 
				<< "\tmap: " << posOmChrom << ":" << posOmFirst << "-" << posOmLast << "->" << posRmChrom << ":" << posRmFirst << "-" << posRmLast << endl; logger.Log(Logger::RESFILE, ss);*/
						
			ss << "REF_POS: " << posRmChrom << ":" << posRmFirst << "-" << posRmLast << endl; logger.Log(Logger::RESFILE, ss);
			ss << "QUALITY: " << mappings[ixMappings].quality << endl; logger.Log(Logger::RESFILE, ss);
			ss << "DP_SCORE: " << mappings[ixMappings].score << endl; logger.Log(Logger::RESFILE, ss);
			std::ostringstream ssAln; 
			ssAln << "ALN: ";

			for (int ixAlignment = 1; ixAlignment < mappings[ixMappings].alignment.size(); ixAlignment++)
			{
				if (ixAlignment > 1) ssAln << ";";
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
					sumOM += length;
					ssOmIxs << " " << ixOMAux - 1;
					ssOmLengths << " " << length;
					ixOMAux++;
					cntOM++;
				}
				int cntRM = 0;
				while (ixRMAux <= mappings[ixMappings].alignment[ixAlignment].second)
				{
					int length = refMaps[chr][ixRMAux - 1].length;
					sumRM += length;
					ssRmPos << " " << refMaps[chr][ixRMAux].chromosome << "_" << refMaps[chr][ixRMAux - 1].start;
					ssRmLengths << " " << length;
					ixRMAux++;
					cntRM++;
				}
				/*ss << "\t" << ixOM << ";" << strings::trim(ssOmIxs.str()) << ";" << strings::trim(ssRmPos.str()) << ";"
					<< strings::trim(ssOmLengths.str()) << ";" << strings::trim(ssRmLengths.str()) << endl; logger.Log(Logger::RESFILE, ss);*/

				ss << ixOMAux + 1 << " - " << ixRMAux + 1 << " (" << sumOM << " - " << sumRM << ") "; logger.Log(Logger::LOGFILE, ss);
				ssAln << sumRM - sumOM << "," << cntRM << ":" << cntOM << "," << sumRM;// << sumOM << " ";
			}
			ss << endl; logger.Log(Logger::LOGFILE, ss);
			ssAln << endl;
			logger.Log(Logger::RESFILE, ssAln);

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
	InitLogging();	
	Parse(optMap, refMaps);
	//InitializeIndex(index, refMap);	
	Mappings *omMatches = AlignOpticalMaps(optMap, refMaps, index);
	SerializeMappings(omMatches, optMap, refMaps);

	return EXIT_SUCCESS;
}