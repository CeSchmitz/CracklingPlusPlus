#include "ExtractOfftargets.hpp"

using std::vector;
using std::string;
using std::fstream;
using std::ofstream;
using std::ifstream;
using boost::algorithm::trim;
using boost::algorithm::to_upper;
using boost::algorithm::split;
using boost::regex;
using boost::sregex_iterator;
using boost::smatch;
namespace fs = std::filesystem;

const regex fwdExp = regex("(?=([ACG][ACGT]{19}[ACGT][AG]G))");
const regex bwdExp = regex("(?=(C[CT][ACGT][ACGT]{19}[TGC]))");

int main(int argc, char** argv)
{
    // Check number of args
    if (argc < 3)
    {
        std::cerr << fmt::format("Usage: {} <output-file>  {<input-file-1> <input-file-2> ... <input-file-n> | <input-dir>}\n", argv[0]) << std::endl;
        exit(1);
    }

    auto startTime = std::chrono::steady_clock::now();
    fs::path tempWorkingDir = fs::path(fs::temp_directory_path() / "Crackling-extractOfftargets");
    fs::create_directory(tempWorkingDir);
    std::atomic_ullong fileCounter = 0;

    std::cout << "Spliting input(s)" << std::endl;
    // Split each sequence to it's own file
    #pragma omp parallel for schedule(static,1)
    for (int i = 2; i < argc; i++)
    {
        // Command line input
        fs::path input(argv[i]);
        vector<string> filesToProcess;
        // File 
        if (fs::is_regular_file(input))
        {
            std::cout << "file" << std::endl;
            filesToProcess.push_back(input.string());
        }
        // Directory
        else if (fs::is_directory(input))
        {
            std::cout << "dir" << std::endl;
            for (const fs::path& file : fs::directory_iterator(input))
            {
                filesToProcess.push_back(file.string());
            }
        }
        // Skip otherwise
        else
        {
            string errorMsg = fmt::format("Error processing: {}\nSkipping entry, Please check the path and try again.");
            std::cout << errorMsg << std::endl;
            continue;
        }

        std::cout << filesToProcess.size() << std::endl;

        for (const string& file : filesToProcess)
        {
            ofstream tempOutFile;
            ifstream inFile;
            string inputLine;
            inFile.open(file, std::ios::binary | std::ios::in);
            std::getline(inFile, inputLine);
            // file is Fasta formatted
            if (inputLine[0] == '>')
            {
                // Remove all line breaks between sequences segments
                tempOutFile.open(tempWorkingDir / fmt::format("{}.txt", fileCounter++), std::ios::binary | std::ios::out);
                for (inputLine; std::getline(inFile, inputLine);)
                {
                    trim(inputLine);
                    to_upper(inputLine);
                    if (inputLine[0] == '>')
                    {
                        tempOutFile.close();
                        tempOutFile.open(tempWorkingDir / fmt::format("{}.txt", fileCounter++), std::ios::binary | std::ios::out);
                    }
                    else
                    {
                        tempOutFile << inputLine;
                    }
                }
            }
            // file is plain text, assume one sequence per line
            else
            {
                do
                {
                    tempOutFile.open(tempWorkingDir / fmt::format("{}.txt", fileCounter++), std::ios::binary | std::ios::out);
                    trim(inputLine);
                    to_upper(inputLine);
                    tempOutFile << inputLine;
                    tempOutFile.close();
                } while (std::getline(inFile, inputLine));
            }
            tempOutFile.close();
            inFile.close();
        }
    }
    std::cout << "Done" << std::endl;

    std::cout << "Sorting intermediate files" << std::endl;
    // Mulithread process each extact and sort
    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < fileCounter; i++)
    {
        ofstream tempOutFile;
        ifstream inFile;
        string inputLine;
        vector<string> offTargets;
        inFile.open(tempWorkingDir / fmt::format("{}.txt", i), std::ios::binary | std::ios::in);
        for (inputLine; std::getline(inFile, inputLine);)
        {
            std::getline(inFile, inputLine);
            trim(inputLine);
            // Add forward matches
            for (sregex_iterator regexItr(inputLine.begin(), inputLine.end(), fwdExp); regexItr != sregex_iterator(); regexItr++)
            {
                offTargets.push_back((*regexItr)[1].str().substr(0,20));
            }
            // Add reverse matches
            for (sregex_iterator regexItr(inputLine.begin(), inputLine.end(), bwdExp); regexItr != sregex_iterator(); regexItr++)
            {
                offTargets.push_back(rc((*regexItr)[1].str()).substr(0, 20));
            }
        }
        inFile.close();
        std::sort(offTargets.begin(), offTargets.end());
        tempOutFile.open(tempWorkingDir / fmt::format("{}_sorted.txt", i), std::ios::binary | std::ios::out);
        for (const string& s : offTargets)
        {
            tempOutFile << s << "\n";
        }
        tempOutFile.close();

    }
    std::cout << "Done" << std::endl;

    std::cout << "Joining intermediate files" << std::endl;
    // Merge sorted files
    vector<ifstream> sortedFiles(fileCounter);
    vector<string> offTargets(fileCounter);
    for (int i = 0; i < fileCounter; i++)
    {
        sortedFiles[i].open((tempWorkingDir / fmt::format("{}_sorted.txt", i)).string());
        string offTarget;
        std::getline(sortedFiles[i], offTarget);
        if (!sortedFiles[i].eof()) {
            offTargets[i] = offTarget;
        }
    }

    ofstream finalOutput;
    finalOutput.open(argv[1], std::ios::binary | std::ios::out);
    while(sortedFiles.size() > 1)
    {
        // Find index of lowest off-target
        int lowest = 0;
        for (int i = 0; i < sortedFiles.size(); i++)
        {
            if (offTargets[i] < offTargets[lowest])
            {
                lowest = i;
            }
        }

        // Write to file
        finalOutput << offTargets[lowest] << "\n";

        // Update offtargets
        std::getline(sortedFiles[lowest], offTargets[lowest]);

        // If at EOF remove from list
        if (sortedFiles[lowest].eof())
        {
            sortedFiles.erase(sortedFiles.begin() + lowest);
            offTargets.erase(offTargets.begin() + lowest);
        }
    }

    while (std::getline(sortedFiles[0], offTargets[0])) 
    {
        finalOutput << offTargets[0] << "\n";
    }
    sortedFiles[0].close();
    finalOutput.close();

    std::cout << "Done" << std::endl;

    std::cout << "Cleaning intermediate files" << std::endl;
    //fs::remove(tempWorkingDir);
    fs::remove_all(tempWorkingDir);
    std::cout << "Done" << std::endl;

    // Report time taken (total)
    std::chrono::nanoseconds nanoSec = std::chrono::steady_clock::now() - startTime;
    std::chrono::duration<uint64_t, std::ratio<86400>> days = std::chrono::duration_cast<std::chrono::duration<uint64_t, std::ratio<86400>>>(nanoSec);
    std::chrono::hours hours = std::chrono::duration_cast<std::chrono::hours>(nanoSec - days);
    std::chrono::minutes minutes = std::chrono::duration_cast<std::chrono::minutes>(nanoSec - days - hours);
    std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(nanoSec - days - hours - minutes);
    std::chrono::milliseconds milliSec = std::chrono::duration_cast<std::chrono::milliseconds>(nanoSec - days - hours - minutes - sec);
    std::chrono::microseconds microSec = std::chrono::duration_cast<std::chrono::microseconds>(nanoSec - days - hours - minutes - sec - milliSec);
    std::cout << fmt::format("Total run time {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days.count(), hours.count(), minutes.count(), sec.count(), std::chrono::duration_cast<std::chrono::seconds>(nanoSec).count()) << std::endl;
}