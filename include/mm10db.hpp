// mm10db.hpp
#pragma once
#include <string>
#include <map>
#include <array>
#include <iterator>
#include <fstream>
#include <list>
#include <regex>
#include <sstream>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Helpers.hpp>

class mm10db
{
public:
	mm10db(ConfigManager& cm);

	void run(std::map<std::string, std::map<std::string, std::string>>& candidateGuides);

	bool static leadingT(std::string candidateGuide);

	float static AT_percentage(std::string candidateGuide);

	bool static polyT(std::string candidateGuide);

	std::string static transToDNA(std::string RNA);

private:
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
	int threadCount;
	std::string RNAFoldOutFile;
	std::string RNAFoldInFile;
	std::string RNAFoldBin;
	float lowEnergyThreshold;
	float highEnergyThreshold;
	int RNAFoldPageLength;
	char printingBuffer[1024];

};