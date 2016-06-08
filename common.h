#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <iostream>

void error_exit(std::string message)
{
	std::cerr << "ERROR: " << message << std::endl;
	exit(EXIT_FAILURE);
}


#endif // COMMON_H