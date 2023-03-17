#include "../include/outputModule.hpp"

using std::ofstream;

outputModule::outputModule(const cracklingConfig& config)
{
	this->outputFile = config.output.filename;
	this->firstWrite = true;
}

void outputModule::run() {};

void outputModule::run(std::vector<guideResults>& candidateGuides)
{
	ofstream outFile;
	outFile.open(outputFile, std::ios::out | std::ios::binary | std::ios::app);
	if (firstWrite)
		outFile << "seq, header, start, strand, end, isUnique, passedG20, passedAvoidLeadingT, passedATPercent, passedTTTT, passedSecondaryStructure, acceptedByMm10db, acceptedBySgRnaScorer2, consensusCount, passedBowtie2, passedOffTargetScore, AT, ssL1, ssStructure, ssEnergy, sgrnascorer2score, bowtie2Chr, bowtie2Start, bowtie2End, mitOfftargetscore, cfdOfftargetscore\n";
		firstWrite = false;

	for (guideResults& candidate : candidateGuides)
	{
		outFile << candidate.seq << "," << candidate.header << ",";
		candidate.start == ULLONG_MAX && candidate.end == ULLONG_MAX ? 
			outFile << CODE_AMBIGUOUS << "," << CODE_AMBIGUOUS << ",": 
			outFile << candidate.start << "," << candidate.end << ",";
		outFile
			<< candidate.strand << ","
			<< candidate.isUnique << ","
			<< candidate.passedG20 << ","
			<< candidate.passedAvoidLeadingT << ","
			<< candidate.passedATPercent << ","
			<< candidate.passedTTTT << ","
			<< candidate.passedSecondaryStructure << ","
			<< candidate.acceptedByMm10db << ","
			<< candidate.acceptedBySgRnaScorer2 << ","
			<< std::to_string(candidate.consensusCount) << ","
			<< candidate.passedBowtie2 << ","
			<< candidate.passedOffTargetScore << ",";
		candidate.AT == -DBL_MAX ?
			outFile << CODE_UNTESTED << ",":
			outFile << candidate.AT << ",";	
		outFile
			<< candidate.ssL1 << ","
			<< candidate.ssStructure << ",";
		candidate.ssEnergy == -DBL_MAX ?
			outFile << CODE_UNTESTED << ",":
			outFile << candidate.ssEnergy << ",";
		candidate.sgrnascorer2score == -DBL_MAX ?
			outFile << CODE_UNTESTED << ",":
			outFile << candidate.sgrnascorer2score << ",";
		outFile << candidate.bowtie2Chr << ",";
		candidate.bowtie2Start == ULLONG_MAX && candidate.bowtie2End == ULLONG_MAX ?
			outFile << CODE_UNTESTED << "," << CODE_UNTESTED << ",":
			outFile << candidate.bowtie2Start << "," << candidate.bowtie2End << ",";
		candidate.mitOfftargetscore == -DBL_MAX ?
			outFile << CODE_UNTESTED << ",":
			outFile << candidate.mitOfftargetscore << ",";
		candidate.cfdOfftargetscore == -DBL_MAX ?
			outFile << CODE_UNTESTED << ",":
			outFile << candidate.cfdOfftargetscore << ",";
		outFile << "\n";
	}
	outFile.close();
}