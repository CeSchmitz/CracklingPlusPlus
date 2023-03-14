// CracklingPlusPlus.cpp : Defines the entry point for the application.
#include "../include/Logger.hpp"
#include "../include/configParserModule.hpp"
#include "../include/cas9InputModule.hpp"
#include "../include/chopchopModule.hpp"
#include "../include/mm10dbModule.hpp"
#include "../include/sgrnascorer2Module.hpp"
#include "../include/bowtie2Module.hpp"
#include "../include/ISSLScoringModule.hpp"
#include "../include/outputModule.hpp"
#include <chrono>

int main(int argc, char** argv)
{
	// Check input arguments
	if (argc != 2) {
		std::cout << fmt::format("Usage: {} [Crackling Config File]\n", argv[0]);
		exit(1);
	}

	// Load config
	configParserModule cmNew; 
	cracklingConfig config = cmNew.run(argv[1]);

	// Create logger objects and apply number formatting
	Logger coutLogger(std::cout, config.output.log.string());
	Logger cerrLogger(std::cerr, config.output.errLog.string());
	std::cout.imbue(comma_locale);
	std::cerr.imbue(comma_locale);

	// Create pipeline modules
	cas9InputModule			cas9IM(config);
	chopchopModule			chopchop(config);
	mm10dbModule			mm10db(config);
	sgrnascorer2Module		sgrnascorer2(config);
	bowtie2Module			bowtie2(config);
	ISSLScoringModule		ISSLScoring(config);
	outputModule			output(config);

	// Record start time
	auto startTime = std::chrono::steady_clock::now();

	cas9IM.run();

	uint64_t batchNum = 0;
	// Process guides in batches
	for (std::vector<guideResults>* currentBatch; currentBatch = cas9IM.next();)
	{
		// Record batch start time
		auto batchStartTime = std::chrono::steady_clock::now();

		std::cout << "Processing batch " << ++batchNum << std::endl;

		//// Consensus scoring
		//chopchop.run(*currentBatch);
		//mm10db.run(*currentBatch);
		//sgrnascorer2.run(*currentBatch);

		//// Complete consensus evaluation
		//std::cout << "Evaluating efficiency via consensus approach." << std::endl;
		//uint64_t failedCount = 0;
		//uint64_t testedCount = 0;
		//for (guideResults& candidate : *currentBatch)
		//{
		//	if (candidate.consensusCount < config.consensus.n) { failedCount++; }
		//	testedCount++;
		//}
		//std::cout << fmt::format("\t{:L} of {:L} failed here.", failedCount, testedCount) << std::endl;

		// Specificity scoring
		//bowtie2.run(*currentBatch);
		ISSLScoring.run(*currentBatch);

		// Print results to file
		output.run(*currentBatch);

		// Report time taken (batch)
		std::chrono::nanoseconds nanoSec = std::chrono::steady_clock::now() - batchStartTime;
		std::chrono::duration<uint64_t, std::ratio<86400>> days = std::chrono::duration_cast<std::chrono::duration<uint64_t, std::ratio<86400>>>(nanoSec);
		std::chrono::hours hours = std::chrono::duration_cast<std::chrono::hours>(nanoSec - days);
		std::chrono::minutes minutes = std::chrono::duration_cast<std::chrono::minutes>(nanoSec - days - hours);
		std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(nanoSec - days - hours - minutes);
		std::chrono::milliseconds milliSec = std::chrono::duration_cast<std::chrono::milliseconds>(nanoSec - days - hours - minutes - sec);
		std::chrono::microseconds microSec = std::chrono::duration_cast<std::chrono::microseconds>(nanoSec - days - hours - minutes - sec - milliSec);
		std::cout << fmt::format("This batch ran in {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days.count(), hours.count(), minutes.count(), sec.count(), std::chrono::duration_cast<std::chrono::seconds>(nanoSec).count()) <<  std::endl;

	}
	// Report time taken (total)
	std::chrono::nanoseconds nanoSec = std::chrono::steady_clock::now() - startTime;
	std::chrono::duration<uint64_t, std::ratio<86400>> days = std::chrono::duration_cast<std::chrono::duration<uint64_t, std::ratio<86400>>>(nanoSec);
	std::chrono::hours hours = std::chrono::duration_cast<std::chrono::hours>(nanoSec - days);
	std::chrono::minutes minutes = std::chrono::duration_cast<std::chrono::minutes>(nanoSec - days - hours);
	std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(nanoSec - days - hours - minutes);
	std::chrono::milliseconds milliSec = std::chrono::duration_cast<std::chrono::milliseconds>(nanoSec - days - hours - minutes - sec);
	std::chrono::microseconds microSec = std::chrono::duration_cast<std::chrono::microseconds>(nanoSec - days - hours - minutes - sec - milliSec);
	std::cout << fmt::format("Total run time {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days.count(), hours.count(), minutes.count(), sec.count(), std::chrono::duration_cast<std::chrono::seconds>(nanoSec).count()) << std::endl;

	coutLogger.close();
	cerrLogger.close();
	cas9IM.cleanup();
}
