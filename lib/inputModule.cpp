#include "../include/inputModule.hpp"

using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::fstream;
using std::vector;
using phmap::flat_hash_set;
using std::filesystem::path;
using std::filesystem::file_size;
using std::filesystem::is_directory;
using std::filesystem::create_directory;
using std::filesystem::remove_all;
using std::filesystem::temp_directory_path;
using boost::algorithm::trim;
using boost::algorithm::to_upper;
using boost::algorithm::split;
using boost::sregex_iterator;
using boost::smatch;



inputModule::inputModule(const cracklingConfig& config)
{
	this->filesToProcess = config.input.filesToProcess;
	this->batchSize = config.input.batchLen;
	this->tempWorkingDir = path(temp_directory_path() / "Crackling");
}

void inputModule::run()
{
	cout << "Analysing files..." << endl;

	// Progress reporting
	uint64_t totalSizeBytes = 0;
	uint64_t completedSizeBytes = 0;
	double completedPercent = 0.0;
	for (const path& file : filesToProcess)
	{
		totalSizeBytes += file_size(file);
	}

	// Setup temp working dir
	if (is_directory(tempWorkingDir)) { remove_all(tempWorkingDir); }
	if (!create_directory(tempWorkingDir)) { throw tempFileSystemError(); }

	cout << fmt::format("Storing batch files in: {}", tempWorkingDir.string()) << endl;

	// Create first batch file
	path outFileName = tempWorkingDir / fmt::format("{}_batchFile.txt", batchFiles.size());
	batchFiles.push_back(outFileName);
	currentBatchFile.open(outFileName, std::ios::binary | std::ios::out);

	guidesInBatch = 0;
	numDuplicateGuides = 0;
	numIdentifiedGuides = 0;

	// Begin processing files
	for (const path& file : filesToProcess)
	{
		cout << fmt::format("Identifying possible target sites in : {}", file.string()) << endl;;
		ifstream inFile;
		string inputLine;

		inFile.open(file, std::ios::binary | std::ios::in);
		std::getline(inFile, inputLine);
		// file is Fasta formatted
		if (inputLine[0] == '>')
		{
			// Remove all line breaks between sequences segments
			fstream tempFasta;
			tempFasta.open(tempWorkingDir / "tempFastsa.fa", std::ios::binary | std::ios::out);
			trim(inputLine);
			tempFasta << inputLine << "\n";
			for (inputLine; std::getline(inFile, inputLine);)
			{
				trim(inputLine);
				if (inputLine[0] == '>')
					tempFasta << "\n" << inputLine << "\n";
				else
					tempFasta << inputLine;
			}
			inFile.close();
			tempFasta.close();

			// Process in pairs
			inFile.open(tempWorkingDir / "tempFastsa.fa", std::ios::binary | std::ios::in);
			for (inputLine; std::getline(inFile, inputLine);)
			{
				trim(inputLine);
				string header = inputLine.substr(1);
				boost::erase_all(header, ",");
				std::getline(inFile, inputLine);
				trim(inputLine);
				string seq = inputLine;
				if (recordedSequences.find(header) == recordedSequences.end())
				{
					recordedSequences.insert(header);
					processSeqeunce(seq, header);
				}
			}
		}
		// file is plain text, assume one sequence per line
		else
		{
			do
			{
				// Clean up input line
				trim(inputLine);
				to_upper(inputLine);

				processSeqeunce(inputLine, file.stem().string());

			} while (std::getline(inFile, inputLine));
		}
		inFile.close();

		// Report overall progress before processing next file
		completedSizeBytes += file_size(file);
		completedPercent = static_cast<double>(completedSizeBytes) / totalSizeBytes * 100.0;

		cout << fmt::format("\tProcessed {:.2f}% of input.", completedPercent) << endl;
	}

	// Finished processing files
	currentBatchFile.close();

	// Report results
	double duplicatePercent = (static_cast<double>(numDuplicateGuides) / static_cast<double>(numIdentifiedGuides)) * 100.0;
	cout << 
		fmt::format(comma_locale,
			"\tIdentified {:L} possible target sites.\n"
			"\tOf these, {:L} are not unique. These sites occur a total of {:L} times.\n"
			"\t{:L} of {:L} ({:.2f}%) of guides will be ignored for optimisation levels over ultralow.\n"
			"\t{:L} distinct guides were identified.",
			numIdentifiedGuides,
			duplicateGuides.size(),
			numDuplicateGuides,
			numDuplicateGuides,
			numIdentifiedGuides,
			duplicatePercent,
			candidateGuides.size()
		)
		<< endl;
	// Remove unneed sets
	candidateGuides.clear();
	recordedSequences.clear();
}

void inputModule::cleanup() {
	remove_all(tempWorkingDir);
}

void inputModule::processSeqeunce(const std::string& seqeunce, const std::string& header)
{
	for (sregex_iterator regexItr(seqeunce.begin(), seqeunce.end(), fwdExp); regexItr != sregex_iterator(); regexItr++)
	{
		numIdentifiedGuides++;
		string guide = (*regexItr)[1].str();
		uint64_t pos = regexItr->position();
		if (candidateGuides.find(guide) == candidateGuides.end())
		{
			candidateGuides.insert(guide);
			if (++guidesInBatch > batchSize)
			{
				currentBatchFile.close();
				path outFileName = tempWorkingDir / fmt::format("{}_batchFile.txt", std::to_string(batchFiles.size()));
				batchFiles.push_back(outFileName);
				currentBatchFile.open(outFileName, std::ios::binary | std::ios::out);
				guidesInBatch = 1;
			}
			currentBatchFile << guide << "," << header << "," << pos << "," << (pos + 23) << "," << "+" << "\n";
		}
		else
		{
			numDuplicateGuides++;
			duplicateGuides.insert(guide);
		}
	}

	for (sregex_iterator regexItr(seqeunce.begin(), seqeunce.end(), bwdExp); regexItr != sregex_iterator(); regexItr++)
	{
		numIdentifiedGuides++;
		string guide = rc((*regexItr)[1].str());
		uint64_t pos = regexItr->position();
		if (candidateGuides.find(guide) == candidateGuides.end())
		{
			candidateGuides.insert(guide);
			if (++guidesInBatch > batchSize)
			{
				currentBatchFile.close();
				path outFileName = tempWorkingDir / fmt::format("{}_batchFile.txt", std::to_string(batchFiles.size()));
				batchFiles.push_back(outFileName);
				currentBatchFile.open(outFileName, std::ios::binary | std::ios::out);
				guidesInBatch = 1;
			}
			currentBatchFile << guide << "," << header << "," << pos << "," << (pos + 23) << "," << "-" << "\n";
		}
		else
		{
			numDuplicateGuides++;
			duplicateGuides.insert(guide);
		}
	}
}

vector<guideResults>* inputModule::next()
{
	if (batchFiles.empty()) return NULL;
	 
	guideBatch.clear();
	guideBatch.reserve(batchSize);
	currentBatchFile.open(batchFiles[0], std::ios::binary | std::ios::in);
	for (string line; std::getline(currentBatchFile, line);)
	{
		trim(line);
		vector<string> splitLine;
		split(splitLine, line, boost::is_any_of(","));
		guideResults newGuide;
		newGuide.seq = splitLine[0];
		if (candidateGuides.find(newGuide.seq) == candidateGuides.end()) {
			newGuide.header = splitLine[1];
			newGuide.start = stoull(splitLine[2]);
			newGuide.end = stoull(splitLine[3]);
			newGuide.strand = splitLine[4][0];
			newGuide.isUnique = true;
		}
		guideBatch.push_back(newGuide);
	}
	currentBatchFile.close();
	batchFiles.erase(batchFiles.begin());
	return &guideBatch;
}
