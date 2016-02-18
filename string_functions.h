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

	inline bool ends_with(std::string const & value, std::string const & ending)
	{
		if (ending.size() > value.size()) return false;
		return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
	}
}



#endif // STRING_FUNCTIONS_H