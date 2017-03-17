/*
* Copyright (C) 2017 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#include <iostream>

#include <string>
#include <sstream>
#include <mutex>

#include "logger.h"

std::mutex mutexLog;

Logger::Logger()
{
}

Logger::~Logger()
{
	if (ofsLog.is_open()) ofsLog.close();
	if (ofsStats.is_open()) ofsStats.close();
	if (ofsResults.is_open()) ofsResults.close();
}

void Logger::InitChannel(Channel channel, std::string filename)
{
	switch (channel)
	{
	case Channel::LOGFILE:
		ofsLog.open(filename);
		break;
	case Channel::STATSFILE:
		ofsStats.open(filename);
		break;
	case Channel::RESFILE:
		if (filename != "")	ofsResults.open(filename);
		break;
	default:
		break;
	}
}

inline void Logger::Log(int channels, std::string message)
{
	if (channels & Channel::STDOUT) std::cout << message;
	if ((channels & Channel::LOGFILE) && ofsLog.is_open()) ofsLog << message;
	if ((channels & Channel::STATSFILE) && ofsStats.is_open()) ofsStats << message;
	if (channels & Channel::RESFILE) ofsResults.is_open() ? ofsResults << message : std::cout << message;
}

void Logger::Log(int channels, std::ostringstream &messageStream, bool clearStream)
{
	mutexLog.lock();
	Log(channels, messageStream.str());
	if (clearStream) messageStream.str(std::string());
	mutexLog.unlock();
}
