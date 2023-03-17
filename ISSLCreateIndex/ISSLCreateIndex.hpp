// ISSLCreateIndex.hpp
/*

Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o isslCreateIndex isslCreateIndex.cpp -O3 -std=c++11 -fopenmp -mpopcnt

*/
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <regex>
#include <filesystem>
#include <fstream>
#include <mutex>
#define FMT_HEADER_ONLY
#include "../include/fmt/format.h"
#include "../include/libpopcnt/libpopcnt.h"