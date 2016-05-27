/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef STRING_FUNCTIONS_H
#define STRING_FUNCTIONS_H

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;

namespace strings{
	inline vector<string> split(string s, string delimiter)
	{
		vector<string> strings;
		size_t pos = 0;
		std::string token;
		while ((pos = s.find(delimiter)) != std::string::npos) {
			token = s.substr(0, pos);
			strings.push_back(token);
			s.erase(0, pos + delimiter.length());
		}
		strings.push_back(s);
		return strings;
	}

	inline std::string trim(const std::string &s)
	{
		std::string::const_iterator it = s.begin();
		while (it != s.end() && isspace(*it))
			it++;

		std::string::const_reverse_iterator rit = s.rbegin();
		while (rit.base() != it && isspace(*rit))
			rit++;

		return std::string(it, rit.base());
	}

	inline std::string upper(const std::string s)
	{
		string sUpper = s;
		std::transform(sUpper.begin(), sUpper.end(), sUpper.begin(), ::toupper);
		return sUpper;
	}

	inline std::string lower(const std::string s)
	{
		string sLower = s;
		std::transform(sLower.begin(), sLower.end(), sLower.begin(), ::tolower);
		return sLower;
	}

	inline bool ends_with(std::string const & str, std::string const & suffix)
	{
		if (suffix.size() > str.size()) return false;
		return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
	}

	inline bool starts_with(std::string const & str, std::string const & prefix)
	{
		if (prefix.size() > str.size()) return false;
		return std::equal(prefix.begin(), prefix.end(), str.begin());;
	}

	inline bool isFloat(string str) 
	{
		std::istringstream iss(str);
		float f;
		iss >> noskipws >> f; // noskipws considers leading whitespace invalid
		// Check the entire string was consumed and if either failbit or badbit is set
		return iss.eof() && !iss.fail();
	}

	inline int toInt(string str)
	{
		try {
			return std::stoi(str);
		}
		catch (...)
		{
			return 0;
		}
		// return isFloat(str) ? std::stoi(str) : 0;
	}

	inline float toFloat(string str)
	{		
		try {
			return std::stof(str);
		}
		catch (...)
		{
			return 0;
		}
		// return isFloat(str) ? std::stoi(str) : 0;
	}
}



#endif // STRING_FUNCTIONS_H