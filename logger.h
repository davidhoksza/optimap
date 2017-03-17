/*
* Copyright (C) 2017 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>

class Logger
{
public:
	enum Channel
	{
		STDOUT = 1 << 0,
		LOGFILE = 1 << 1,
		RESFILE = 1 << 2,
		STATSFILE = 1 << 3
	};

	Logger();
	~Logger();
	void InitChannel(Channel channel, std::string filename);
	void Log(int channels, std::string message);
	void Log(int channels, std::ostringstream &messageStream, bool clearStream = true);

private:
	std::ofstream ofsLog, ofsStats, ofsResults;

};

#endif // LOGGER_H