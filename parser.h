#ifndef PARSE_H
#define PARSE_H

#include <fstream>

#include "common.h"
#include "string_functions.h"
#include "gzstream/gzstream.h"
#include "types.h"

using namespace std;


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

RefMaps parse_ref_map(string fileName, string ch, float &avgRefLength)
{
	RefMaps refMaps;

	istream *ifs = open_map_file(fileName);

	string line;
	avgRefLength = 0;
	int cnt = 0;
	float auxLength;
	for (string line; getline(*ifs, line);) {
		cnt++;
		vector<string> strs = strings::split(line, "\t");

		RMRead auxRMRead;
		auxRMRead.chromosome = strings::trim(strs[0]);
		if (ch != "" && strings::upper(auxRMRead.chromosome) != strings::upper(ch)) continue;
		auxRMRead.start = stof(strs[1]);
		auxLength = stof(strs[3]);
		if (auxLength < 0) auxLength = 0; //if length two restriction sites overlap
		auxRMRead.length = auxLength * 1000;
		avgRefLength += auxLength;
		if (refMaps.count(auxRMRead.chromosome) == 0) refMaps[auxRMRead.chromosome] = vector<RMRead>();
		refMaps[auxRMRead.chromosome].push_back(auxRMRead);
	}
	avgRefLength /= cnt;
	avgRefLength *= 1000;

	delete ifs;
	return refMaps;
}

vector<ExpMap> parse_exp_map(string fileName, string format, int topN = numeric_limits<int>::max())
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
		if (format == "opgen")
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
		if (format == "bionano")
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


#endif // PARSE_H