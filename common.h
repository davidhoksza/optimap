#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <iostream>
#include <logger.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h> 


#include "types.h"
#include "logger.h"
#include "constants.h"

void error_exit(std::string message);

extern Params params;
extern Logger logger;


#endif // TYPES_H