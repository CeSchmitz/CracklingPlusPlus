#include "../include/mm10dbModule.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ios;
using std::ofstream;
using std::ifstream;
using boost::algorithm::starts_with;
using boost::algorithm::ends_with;
using boost::algorithm::trim_right;
using boost::algorithm::split;
using boost::regex;
using boost::smatch;
using boost::regex_search;

const vector<char> AT_vector = { 'A', 'T' };
const string guide = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU";
const regex pattern_RNAstructure(".{28}\\({4}\\.{4}\\){4}\\.{3}\\){4}.{21}\\({4}\\.{4}\\){4}\\({7}\\.{3}\\){7}\\.{3}\\s\\((.+)\\)");
const regex pattern_RNAenergy("\\s\\((.+)\\)");

mm10dbModule::mm10dbModule(cracklingConfig config) : consensusModule(config)
{
	this->toolIsSelected = config.consensus.mm10db;
	this->rnafoldConfig = config.rnafold;
}

void mm10dbModule::run(std::vector<guideResults>& candidateGuides)
{
	if (!toolIsSelected)
	{
		cout << "mm10db has been configured not to run. Skipping mm10db" << endl;
		return;
	}

	cout << "mm10db - remove all targets with a leading T(+) or trailing A(-)." << endl;
	uint64_t failedCount = 0;
	uint64_t testedCount = 0;
	for (guideResults& candidate : candidateGuides)
	{
		// Run time filtering
		if (!processGuide(candidate)) { continue; }

		if (leadingT(candidate.seq))
		{
			candidate.passedAvoidLeadingT = CODE_REJECTED;
			failedCount++;
		}
		else
		{
			candidate.passedAvoidLeadingT = CODE_ACCEPTED;
		}
		testedCount++;
	}
	cout << fmt::format("\t{} of {} failed here.", failedCount, testedCount) << endl;

	cout << "mm10db - remove based on AT percent." << endl;
	failedCount = 0;
	testedCount = 0;
	for (guideResults& candidate : candidateGuides)
	{
		// Run time filtering
		if (!processGuide(candidate)) { continue; }

		double AT = AT_percent(candidate.seq.substr(0,20));
		if (AT < 20.0 || AT > 65.0)
		{
			candidate.passedATPercent = CODE_REJECTED;
			failedCount++;
		}
		else
		{
			candidate.passedATPercent = CODE_ACCEPTED;
		}
		candidate.AT = AT;
		testedCount++;
	}
	cout << fmt::format("\t{} of {} failed here.", failedCount, testedCount) << endl;

	cout << "mm10db - remove all targets that contain TTTT." << endl;
	failedCount = 0;
	testedCount = 0;
	for (guideResults& candidate : candidateGuides)
	{
		// Run time filtering
		if (!processGuide(candidate)) { continue; }

		if (polyT(candidate.seq))
		{
			candidate.passedTTTT = CODE_REJECTED;
			failedCount++;
		}
		else
		{
			candidate.passedTTTT = CODE_ACCEPTED;
		}
		testedCount++;
	}
	cout << fmt::format("\t{} of {} failed here.", failedCount, testedCount) << endl;

	cout << "mm10db - check secondary structure." << endl;

	testedCount = 0;
	failedCount = 0;
	uint64_t errorCount = 0;
	uint64_t guidesInPage = 0;
	uint64_t pgIdx = 1;
	auto paginatorIterator = candidateGuides.begin();
	auto pageStart = candidateGuides.begin();
	auto pageEnd = candidateGuides.begin();

	// Outer loop deals with changing iterator start and end points (Pagination)
	while (pageEnd != candidateGuides.end())
	{
		if (rnafoldConfig.pageLen > 0)
		{
			// Advance the pageEnd pointer
			std::advance(pageEnd, std::min((uint64_t)std::distance(pageEnd, candidateGuides.end()), rnafoldConfig.pageLen));
			// Record page start
			pageStart = paginatorIterator;
			// Print page information
			cout << fmt::format("\tProcessing page {} ({} per page).", pgIdx, rnafoldConfig.pageLen) << endl;
		}
		else {
			// Process all guides at once
			pageEnd = candidateGuides.end();
		}
		cout << "\t\tConstructing the RNAfold input file." << endl;

		// Open input file 
		ofstream out;
		out.open(rnafoldConfig.inFile, ios::binary | ios::out);

		guidesInPage = 0;
		while (paginatorIterator != pageEnd)
		{
			// Run time filtering
			if (!processGuide(*paginatorIterator)) {
				// Advance page end
				if (pageEnd != candidateGuides.end())
				{
					pageEnd++;
				}
				paginatorIterator++;
				continue;
			}
			out << "G" << paginatorIterator->seq.substr(1, 19) << guide << "\n";
			guidesInPage++;
			paginatorIterator++;
		}

		out.close();

		cout << fmt::format("\t\t{} guides in this page.", guidesInPage) << endl;

		// Call RNAFold
		runner(fmt::format("{} --noPS -j{} -i {} > {}", rnafoldConfig.binary.string(), rnafoldConfig.threads, rnafoldConfig.inFile.string(), rnafoldConfig.outFile.string()).c_str());

		cout << "\t\tStarting to process the RNAfold results." << endl;

		// Reset paginatorIterator to page start
		paginatorIterator = pageStart;

		// Open output file
		ifstream in;
		in.open(rnafoldConfig.outFile, ios::binary | ios::in);

		while (paginatorIterator != pageEnd)
		{
			// Run time filtering
			if (!processGuide(*paginatorIterator)) {
				paginatorIterator++;
				continue;
			}

			// Results are across two lines L1 and L2
			string L1;
			string L2;
			string target;
			std::getline(in, L1);
			std::getline(in, L2);
			trim_right(L1);
			trim_right(L2);
			target = L1.substr(0, 20);
			vector<string> L2_split;
			split(L2_split, L2, boost::is_any_of(" "));

			paginatorIterator->ssL1 = L1;
			paginatorIterator->ssStructure = L2_split[0];
			paginatorIterator->ssEnergy = std::stod(L2_split[1]);

			if ((transToDNA(target) != paginatorIterator->seq.substr(0, 20)) &&
				((transToDNA("C" + target.substr(1))) != paginatorIterator->seq.substr(0, 20)) &&
				((transToDNA("A" + target.substr(1))) != paginatorIterator->seq.substr(0, 20)))
			{
				paginatorIterator->passedSecondaryStructure = CODE_ERROR;
				errorCount++;
				paginatorIterator++;
				continue;
			}

			smatch match_structure;
			smatch match_energy;
			if (regex_search(L2, match_structure, pattern_RNAstructure))
			{
				double energy = std::stod(match_structure[1].str());
				if (energy < rnafoldConfig.lowEngeryThreshold)
				{
					paginatorIterator->passedSecondaryStructure = CODE_REJECTED;
					failedCount++;
				}
				else
				{
					paginatorIterator->passedSecondaryStructure = CODE_ACCEPTED;
				}
			}
			else if (regex_search(L2, match_energy, pattern_RNAenergy))
			{
				double energy = std::stod(match_energy[1].str());
				if (energy <= rnafoldConfig.highEngeryThreshold)
				{
					paginatorIterator->passedSecondaryStructure = CODE_REJECTED;
					failedCount++;
				}
				else
				{
					paginatorIterator->passedSecondaryStructure = CODE_ACCEPTED;
				}
			}
			testedCount++;
			paginatorIterator++;
		}

		// Advance paginatorIterator to page end for next loop
		paginatorIterator = pageEnd;
		pgIdx++;
	}

	cout << fmt::format("\t{} of {} failed here.", failedCount, testedCount) << endl;
	if (errorCount > 0)
	{
		cout << fmt::format("\t{} of {} errored here.", errorCount, testedCount) << endl;
	}

	cout << "Calculating mm10db final result." << endl;
	uint64_t acceptedCount = 0;
	failedCount = 0;
	for (guideResults& candidate : candidateGuides)
	{
		if ((candidate.passedAvoidLeadingT == CODE_ACCEPTED) &&
			(candidate.passedATPercent == CODE_ACCEPTED) &&
			(candidate.passedTTTT == CODE_ACCEPTED) &&
			(candidate.passedSecondaryStructure == CODE_ACCEPTED))
		{
			candidate.acceptedByMm10db = CODE_ACCEPTED;
			acceptedCount++;
		}
		else
		{
			candidate.acceptedByMm10db = CODE_REJECTED;
			failedCount++;
		}
	}
	cout << fmt::format("\t{} accepted.\n\t{} rejected", acceptedCount, failedCount) << endl;
}


bool mm10dbModule::leadingT(std::string_view guide)
{
	return	(starts_with(guide, "T") || ends_with(guide, "A"));
}

double mm10dbModule::AT_percent(std::string_view guide)
{
	double total = 0.0;
	double length = (double)guide.size();

	for (char c : guide) {
		// Check if the char is present in the 'AT' array
		if (std::find(AT_vector.begin(), AT_vector.end(), c) != AT_vector.end())
		{
			total += 1.0;
		}
	}

	return (100.0 * total / length);
}

bool mm10dbModule::polyT(std::string_view guide)
{
	for (int i = 0; i < guide.size() - 4; i++)
	{
		if (guide.substr(i, 4) == "TTTT")
		{
			return true;
		}
	}
	return false;
}

string mm10dbModule::transToDNA(string RNA)
{
	// Swap U with T
	std::replace(RNA.begin(), RNA.end(), 'U', 'T');
	return RNA;
}

bool mm10dbModule::processGuide(const guideResults& guide)
{
	// Process all guides at this level
	if (optimsationLevel == optimisationLevel::ultralow) { return true; }

	// For all levels above `ultralow`
	if (!guide.isUnique)
	{
		// Reject all guides that have been seen more than once
		return false;
	}

	// For optimisation levels `medium` and `high`
	if (optimsationLevel == optimisationLevel::medium || optimsationLevel == optimisationLevel::high)
	{
		bool failedMM10DBStage =
			guide.passedAvoidLeadingT == CODE_REJECTED ||
			guide.passedATPercent == CODE_REJECTED ||
			guide.passedTTTT == CODE_REJECTED ||
			guide.passedSecondaryStructure == CODE_REJECTED;
		// Stop testing a guide if it has failed any stage of mm10db
		if (failedMM10DBStage) { return false; }
	}

	// For optimisation level `high`
	if (optimsationLevel = optimisationLevel::high)
	{
		int countAlreadyAccepted =
			(int)guide.passedG20 == CODE_ACCEPTED +
			(int)guide.acceptedByMm10db == CODE_ACCEPTED +
			(int)guide.acceptedBySgRnaScorer == CODE_ACCEPTED;

		int countAlreadyAssessed =
			(int)guide.passedG20 != CODE_UNTESTED +
			(int)guide.acceptedByMm10db != CODE_UNTESTED +
			(int)guide.acceptedBySgRnaScorer != CODE_UNTESTED;

		// Reject if the consensus has already been passed
		if (countAlreadyAccepted >= consensusN) { return false; }

		// Reject if there is not enough tests remaining to pass consensus
		if (toolCount - countAlreadyAssessed < consensusN - countAlreadyAccepted) { return false; }
	}

	// None of the failure conditions have been meet, return true
	return true;
}