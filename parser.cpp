/*
* Copyright (C) 2017 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#include <assert.h>
#include <time.h>

#include "parser.h"
#include "string_functions.h"
#include "gzstream/gzstream.h"

using namespace std;

static ostringstream ss; //static because the same definition is used in different files (by default the implicit behaivour is extern linkage)

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

RefMaps parse_ref_map(string fileName/*, float &avgRefLength*/)
{
	RefMaps refMaps;

	istream *ifs = open_map_file(fileName);

	string line;
	//avgRefLength = 0;
	int cnt = 0;
	float auxLength;
	for (string line; getline(*ifs, line);) {
		cnt++;
		vector<string> strs = strings::split(line, "\t");

		RMRead auxRMRead;
		auxRMRead.chromosome = strings::trim(strs[0]);
		if (params.chromosome != "" && strings::upper(auxRMRead.chromosome) != strings::upper(params.chromosome)) continue;
		auxRMRead.start = stof(strs[1]);
		auxLength = stof(strs[3]);
		if (auxLength < 0) auxLength = 0; //if length two restriction sites overlap
		auxRMRead.length = auxLength * 1000;
		//avgRefLength += auxLength;
		if (refMaps.count(auxRMRead.chromosome) == 0) refMaps[auxRMRead.chromosome] = vector<RMRead>();
		refMaps[auxRMRead.chromosome].push_back(auxRMRead);
	}
	//avgRefLength /= cnt;
	//avgRefLength *= 1000;

	delete ifs;
	return refMaps;
}

vector<ExpMap> parse_exp_map(string fileName, int topN)
{
	ostringstream ss;
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

void Parse(vector<ExpMap> &expMaps, RefMaps &refMaps)
{
	ss << "======= PARSE - START =======" << endl;
	logger.Log(Logger::LOGFILE, ss);
	clock_t begin_time = clock();
	//vector<ExpMap> expMap = parse_exp_map("../CASTEiJ_Alldata.maps", 1000);
	//vector<ExpMap> expMap = parse_exp_map("../ref.map.split", 100);
	expMaps = parse_exp_map();
	refMaps = parse_ref_map();
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