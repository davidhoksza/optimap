/*
* Copyright (C) 2017 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#include "common.h"

void error_exit(std::string message)
{
	std::cerr << "ERROR: " << message << std::endl;
	exit(EXIT_FAILURE);
}