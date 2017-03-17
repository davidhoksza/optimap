/*
* Copyright (C) 2016 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#define _SCL_SECURE_NO_WARNINGS

//#include "indexing.h"
#include "common.h"
#include "string_functions.h"
#include "parser.h"
#include "mapping.h"

#include "tclap/CmdLine.h"

using namespace std;

static ostringstream ss;

Params params;
Logger logger;

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

//void PrecomputeScores(vector<ExpMap> &expMap, RefMaps &refMaps)
//{
//	if (params.errorModel == "valuev-lr")
//	{
//		for (int ixRefCnt; ixRefCnt <= params.maxDpWindowSize; ixRefCnt++)
//		{
//			for (int ixExpCnt; ixExpCnt <= params.maxDpWindowSize; ixExpCnt++)
//			{
//				for (int ixRefLen = 0; ixRefLen <= (maxSeqLength * params.maxDpWindowSize) / stepLength; ixRefLen++)
//				{
//					for (int ixExpLen = 0; ixExpLen <= (maxSeqLength * params.maxDpWindowSize) / stepLength; ixExpLen++)
//					{
//						SCORE_TYPE score = score_segment(ixExpLen, ixRefLen, ixExpCnt, ixRefCnt, -1, -1, expMap[0].reads, refMaps.begin()->second);
//					}
//				}
//			}
//		}
//	}
//}

void InitLogging()
{
	if (!params.logFileName.empty())	logger.InitChannel(Logger::LOGFILE, params.logFileName);
	if (params.outFileName != "") logger.InitChannel(Logger::RESFILE, params.outFileName);
}

int main(int argc, char** argv)
{
	vector<ExpMap> expMap;
	RefMaps refMaps;

	ParseCmdLine(argc, argv);	
	InitLogging();
	Parse(expMap, refMaps);
	//PrecomputeScores(expMap, refMaps);
	//clock_t begin_time = clock();	
	Mappings *omMatches = AlignOpticalMaps(expMap, refMaps);
	//cout << endl << "Time(s): " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl; 
	SerializeMappings(omMatches, expMap, refMaps);
	delete[] omMatches;

	return EXIT_SUCCESS;
}
