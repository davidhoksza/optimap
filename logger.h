/*
* Copyright (C) 2015 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <mutex>

std::mutex mutexLog;

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
	void Log(int channels, std::ostringstream &messageStream);

private:
	std::ofstream ofsLog, ofsStats, ofsResults;

};

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
		ofsResults.open(filename);
		break;
	default:
		break;
	}
}

void Logger::Log(int channels, std::string message)
{
	if (channels & Channel::STDOUT) std::cout << message;
	if (channels & Channel::LOGFILE) ofsLog << message;
	if (channels & Channel::STATSFILE) ofsStats << message;
	if (channels & Channel::RESFILE) ofsResults << message;
}

void Logger::Log(int channels, std::ostringstream &messageStream)
{
	mutexLog.lock();
	Log(channels, messageStream.str());
	messageStream.str(std::string());
	mutexLog.unlock();
}

Logger logger;

#endif // LOGGER_H