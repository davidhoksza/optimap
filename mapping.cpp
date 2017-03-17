/*
* Copyright (C) 2017 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#include <assert.h>  

#include <iomanip>
#include <algorithm>
#include <thread>
#include <mutex>


#include "mapping.h"
#include "string_functions.h"
#include "stats.h"

using namespace std;



mutex mutexOM;
long unsigned int omProcessed = 0;

static ostringstream ss;
vector<int> candidateSectionLengths;
vector<int> scoresCalculations;
unsigned long int scoresComputed = 0;

//double *scoresCache;
//int		stepLength = 5;
//int		maxSeqLength = 200000;


float avgRefLength;

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
		float stddev = sqrt(float(refLength)) * (refLength > params.smallFragmentThreshold ? 5 : 7);

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
	else if (params.errorModel == "valuev-lr")
	{
		/********

		The score is computed as likelihood ratio as described in Valuev et. al. (DOI: 10.1089/cmb.2006.13.442).
		Score of each aligned region consists of two likelihood ratios:
		1. LR for matching the region of size x com prised of m of optical map fragments
		to the region of size y comprised of n reference map fragments isone representing sizing error
		2. LR for the matching region comprised of m optical map fragments and the reference region
		comprised of n reference map fragments given the reference region of size y

		S(q_i − q_g, r_j −r_h, i − g, j − h) = −log(LR(q_i − q_g; r_j −r_h,i − g, j − h)) − log(LR(i − g; r_j − r_h, j − h))

		********/
		//float sigma = refLength > params.smallFragmentThreshold ? 5 : 7;  //standard deviation of error
		//float sigma2 = refLength > params.smallFragmentThreshold ? 25 : 49;  //standard deviation of error
		//float twosigma2 = refLength > params.smallFragmentThreshold ? 50 : 98;  //standard deviation of error
		//float sqrt_2pi_sigma = stats::sqrt_2pi *sigma;		
		static float sigma_short = 5;
		static float sigma_long = 7;
		static float twosigma2_short = 2 * sigma_short *sigma_short;
		static float twosigma2_long = 2 * sigma_long *sigma_long;
		static float sqrt_2pi_sigma_short = stats::sqrt_2pi *sigma_short;
		static float sqrt_2pi_sigma_long = stats::sqrt_2pi *sigma_long;
		//float sigma = refLength > params.smallFragmentThreshold ? sqrt(0.3) : sqrt(5);  //standard deviation of error
		static float zeta = 0.0000065; //breakage rate
		static float dgst_prob = 1 - params.missRestrictionProb;
		static float lambda = avgRefLength; //mean of reference fragments (,which have exponential density)		
		static float tau = 1 / (zeta + dgst_prob / lambda);
		static float theta_short = 1 / (1 / sigma_short *sqrt(2 / tau + 1 / (sigma_short *sigma_short)) - 1 / (sigma_short*sigma_short)); //mean of fragment sizes of experimental maps (,which have exponential density)		
		static float theta_long = 1 / (1 / sigma_long *sqrt(2 / tau + 1 / (sigma_long *sigma_long)) - 1 / (sigma_long*sigma_long)); //mean of fragment sizes of experimental maps (,which have exponential density)		
		static float f_M_m = 1.0 / params.maxDpWindowSize;

		//float lr_size, lr_cnt;
		float f_size_H0, f_size_HA, f_cnt_H0, f_cnt_HA;


		if (refLength > 4000)
		{
			/*lr_size = (stats::sqrt_2pi * sqrt(refLength) * sigma * pow(expLength, cntExpFrags - 1)) / (stats::factorial(cntExpFrags - 1) * pow(theta, cntExpFrags)) *
			exp(((expLength - refLength) * (expLength - refLength)) / (twosigma2 * refLength) - expLength / theta);*/
			f_size_H0 = pow((float)expLength, cntExpFrags - 1)*exp(-(float)expLength / theta_long) / (stats::factorial(cntExpFrags - 1)*pow(theta_long, cntExpFrags));
			f_size_HA = exp(-(((float)expLength - refLength)*(expLength - refLength)) / (twosigma2_long*refLength)) / (sqrt_2pi_sigma_long * sqrt((float)refLength));
		}
		else
		{
			/*lr_size = (stats::sqrt_2pi * sigma * pow(expLength, cntExpFrags - 1)) / (stats::factorial(cntExpFrags - 1) * pow(theta, cntExpFrags)) *
			exp(((expLength - refLength) * (expLength - refLength)) / (twosigma2)-expLength / theta);*/
			f_size_H0 = pow((float)expLength, cntExpFrags - 1)*exp(-(float)expLength / theta_short) / (stats::factorial(cntExpFrags - 1)*pow(theta_short, cntExpFrags));
			f_size_HA = exp(-(((float)expLength - refLength)*(expLength - refLength)) / (twosigma2_short)) / (sqrt_2pi_sigma_short);
		}

		//cout << f_size_H0 << ";" << f_size_HA << endl;

		//lr_cnt = (exp(zeta*refLength)*stats::factorial(cntExpFrags - 1)*f_M_m) / (pow(1 - dgst_prob, cntRefFrags - 1) * pow(zeta*refLength, cntExpFrags - 1));
		f_cnt_H0 = f_M_m;
		f_cnt_HA = (pow(1 - dgst_prob, cntRefFrags - 1) * exp(-zeta*refLength) * pow(zeta*refLength, cntExpFrags - 1)) / stats::factorial(cntExpFrags - 1);


		/*if (isinf(lr_size) || isinf(lr_cnt)) score = stats::transform_prob(0);
		else score = -log(lr_size) -log(lr_cnt);*/
		if (isinf(f_size_HA) || f_size_HA == 0 || isinf(f_size_H0) || f_size_H0 == 0 || isinf(f_cnt_HA) || f_cnt_HA == 0 || isinf(f_cnt_H0) || f_cnt_H0 == 0) score = stats::transform_prob(0);
		else
		{
			//score = stats::transform_prob(f_size_HA) - stats::transform_prob(f_size_H0) + stats::transform_prob(f_cnt_HA) - stats::transform_prob(f_cnt_H0);
			score = -log(f_size_HA / f_size_H0) - log(f_cnt_HA / f_cnt_H0);
			//score = - (log(f_size_HA) - log(f_size_H0)) - (log(f_cnt_HA) -log(f_cnt_H0));
			//cout << score << ";" << f_size_HA / f_size_H0 << endl;
		}

		//cout << score << endl;

		//exit(0);
	}
	else if (params.errorModel == "li")
	{
		// Sizing error

		//float location, scale;
		int laplace_type;
		if (refLength < 2400)
		{
			//location = 0.858181;
			//scale = 0.180196;
			laplace_type = 0;
		}
		else if (refLength < 3600)
		{
			//location = 0.980760;
			//scale = 0.071176;
			laplace_type = 1;
		}
		else if (refLength < 4800)
		{
			//location = 1.003354;
			//scale = 0.052800;
			laplace_type = 2;
		}
		else
		{
			//location = 1.00482;
			//scale = 0.042428;
			laplace_type = 3;
		}

		//score += stats::transform_prob(stats::pdf_laplace_full(expLength/(float)refLength, location, scale));
		//float aux = stats::pdf_laplace_full(expLength / (float)refLength, location, scale);
		//float aux = stats::pdf_laplace(expLength / (float)refLength, laplace_type);
		//cout << aux << ";";
		score += stats::transform_prob(stats::pdf_laplace(expLength / (float)refLength, laplace_type));

		// Missing cuts + aligned cut
		float digestion_rate;
		SCORE_TYPE auxP = 1;
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
				float dAvg = (d1 + d2) / 2.0; //unit length is taken to be 1.2kb
				if (dAvg > 18000) digestion_rate = 0.9;
				else {
					dAvg /= 1200;
					digestion_rate = 0.0003089 * dAvg * dAvg * dAvg - 0.01069 * dAvg * dAvg + 0.1253 * dAvg + 0.3693;
				}

			}
			if (ix < cntRefFrags) auxP *= 1 - digestion_rate;
			else auxP *= digestion_rate;
		}
		score += stats::transform_prob(auxP);

		//False cuts
		//probability of seeing given number of false cuts
		//float lambdaFactor = expLength / 200000.0; //the model takes 200kb as a unit		
		//auxP = 0.18 * stats::pdf_poisson_full(cntExpFrags - 1, 0) + 0.6 * stats::pdf_poisson_full(cntExpFrags - 1, 1 * lambdaFactor) + 0.22 * stats::pdf_poisson_full(cntExpFrags - 1, 3 * lambdaFactor);
		auxP = 0.18 * stats::pdf_poisson_200kb(0, cntExpFrags - 1, expLength) + 0.6 * stats::pdf_poisson_200kb(1, cntExpFrags - 1, expLength) + 0.22 * stats::pdf_poisson_200kb(3, cntExpFrags - 1, expLength);
		//now this needs to be modified based on where the given cuts are
		//TODO
		//location of false cut is modeled as hybrid of three distributions
		// U[0.1, 0.9], 0.1 ≤ lfp ≤ 0.9, w.p. 0.8852
		// N(0.1, 0.044186), lfp < 0.1, w.p. 0.0574
		// N(0.9, 0.044186), lfp > 0.9, w.p. 0.0574
		score += stats::transform_prob(auxP);
	}

	return score;
}

void dp_fill_matrix(DpMatrixCell ** matrix, vector<int> &experiment, std::vector<RMRead> &reference, vector<SCORE_TYPE> &minScoresSoFar)
{
	//first, let's initiliaze the first column with submax values which ensures 
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
						score += score_segment(rowValue, colValue, ixWindowRow, ixWindowCol, ixRow - 1, ixCol - 1, experiment, reference);
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
			if (params.errorModel != "valuev-lr" && isLastRow /*&& ixCol >= ixColResultFrom*/ && minCell.value < minScoresSoFar[0])
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
		//int threshold = INDEX_NEIGHBORHOOD_THRESHOLD;

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

Mappings* AlignOpticalMaps(vector<ExpMap> &expMap, RefMaps &refMaps)
{
	ss << "======= PROCOMPUTING STATS - START =======" << endl; logger.Log(Logger::LOGFILE, ss);
	stats::init_stats(params.falseCutProb);
	ss << "======= PROCOMPUTING STATS - END =======" << endl; logger.Log(Logger::LOGFILE, ss);

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

	ss << endl << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl; logger.Log(Logger::STDOUT, ss);
	ss << "======= ALIGNING - END =======" << endl; logger.Log(Logger::LOGFILE, ss);

	return omMappings;
}