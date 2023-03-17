#include "ISSLCreateIndex.hpp"

using std::vector;
using std::map;
using std::pair;
using std::string;
using std::ifstream;
using std::regex;
using std::regex_iterator;
using std::filesystem::path;
using std::filesystem::exists;
using std::filesystem::file_size;

const regex extractNumbers("[1234567890]+");
const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };
// 3bit
// const vector<uint8_t> nucleotideIndex{ 1,0,2,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,7 };
// const vector<char> signatureIndex{ '0','A','C','3','G','5','6','T' };
uint64_t seqLength;

/**
 * getLineEnding
 *
 * Determines what line ending the provided file is using.
 *
 * @param file The file to check
 * @return The linefeed charater(s) that the file is using.
 */
string getLineEnding(const path& file)
{
    string lineEnding;
    char currentChar;
    ifstream inFile;
    inFile.open(file);
    inFile.ignore(seqLength);
    for (;;)
    {
        inFile.read(&currentChar, 1);
        if (currentChar == '\r' || currentChar == '\n')
        {
            lineEnding.push_back(currentChar);
        }
        else
        {
            break;
        }  
    }
    inFile.close();
    return lineEnding;
}


/**
 * encode3Bit
 *
 * Converts from sequence of letters to three bit encoded signature
 *
 * @param ptr A pointer to the beginning of a char seqeuence
 * @return The three bit encoded signature of the sequence at `ptr`
 */
uint64_t encode3Bit(const char* ptr)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLength; j++) {
        signature |= static_cast<uint64_t>(nucleotideIndex[(ptr[j]) - 65]) << (j * 3);
    }
    return signature;
}

/**
 * decode3Bit
 *
 * Converts from three bit encoded signature to sequence of letters
 *
 * @param signature The three bit encoded signature
 * @return The letter sequence of `signature`
 */
string decode3Bit(uint64_t signature)
{
    string sequence = string(seqLength, ' ');
    for (uint64_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 3)) & 0x7];
    }
    return sequence;
}

/**
 * decode3Bit
 *
 * Converts from three bit encoded signature to sequence of letters
 *
 * @param signature The three bit encoded signature
 * @return The letter sequence of `signature`
 */
string decode3Bit(uint64_t signature, int len)
{
    string sequence = string(len, ' ');
    for (uint64_t j = 0; j < len; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 3)) & 0x7];
    }
    return sequence;
}

/**
 * computeMasksThreeBit
 *
 * This function generates all of the three bit combinations of mismatches.
 *
 * @param seqLength, The length of the seqeunce
 * @param mismatches, The number of mismatches allowed
 * @return A vector of two bit combinations
 */
vector<uint64_t> computeMasksThreeBit(int seqLength, int mismatches) {
    vector<uint64_t> masks;

    if (mismatches != seqLength) {
        if (mismatches > 0) {
            for (auto mask : computeMasksThreeBit(seqLength - 1, mismatches - 1)) {
                masks.push_back((1ULL << (seqLength - 1) * 3) + mask);
            }

            for (auto mask : computeMasksThreeBit(seqLength - 1, mismatches)) {
                masks.push_back(mask);
            }
        }
        else {
            masks.push_back(0ULL);
        }
    }
    else {
        uint64_t tempMask = 0;
        for (int i = 0; i < seqLength; i++) {
            tempMask |= (1ULL << i * 3);
        }
        masks.push_back(tempMask);
    }
    return masks;
}

/**
 * sequenceToSignature
 * 
 * Converts from sequence of letters to two bit encoded signature
 * 
 * @param ptr, A pointer to the beginning of a char seqeuence
 * @return The two bit encoded signature of the sequence at `ptr`
 */
uint64_t sequenceToSignature(const char* ptr)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLength; j++) {
        signature |= static_cast<uint64_t>(nucleotideIndex[*ptr]) << (j * 2);
        ptr++;
    }
    return signature;
}

/**
 * sequenceToSignature
 * 
 * Converts from two bit encoded signature to sequence of letters
 * 
 * @param signature, The two bit encoded signature
 * @return The letter sequence of `signature`
 */
string signatureToSequence(uint64_t signature)
{
    string sequence = string(seqLength, ' ');
    for (uint64_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}

/**
 * computeMasksTwoBit
 *
 * This function generates all of the two bit combinations of mismatches.
 *
 * @param seqLength, The length of the seqeunce 
 * @param mismatches, The number of mismatches allowed 
 * @return A vector of two bit combinations
 */
vector<uint64_t> computeMasksTwoBit(uint64_t seqLength, uint64_t mismatches) {
    vector<uint64_t> masks;

    // There will only be one valid combination when mismatches == seqLength.
    if (mismatches != seqLength) {

        // There are more mismatchs to assign.
        if (mismatches > 0) {
            // This loop reduces assigns a mismatch and reduces the remaining free positions and calls itself.
            for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches - 1)) {

                masks.push_back((1ULL << (seqLength - 1) * 2) + mask);
            }

            // In order to repeat the first loop process with a new sequence, shift the position of the first mismatch.
            for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches)) {
                masks.push_back(mask);
            }
        }
        else {
            // No more mismatches allowed, return a mask of 0 (Exit condition).
            masks.push_back(0ULL);
        }
    }
    /**
    * The only valid combination when mismatches == seqLength will be to assign remaining seq as mismatches.
    * E.g. mismatches=4 and seqLength=4, we want: 10 10 10 10
    * Do not recursivley call self as there are no more valid combinations from this point (Exit condition).
    */ 
    else {
        uint64_t tempMask = 0;
        for (int i = 0; i < seqLength; i++) {
            tempMask |= (1ULL << i * 2);
        }
        masks.push_back(tempMask);
    }
    return masks;
}

/**
 * calcMITLocalScore
 *
 * This function calculates the local MIT score based on the positions of mismatches.
 * https://dx.doi.org/10.1038/nbt.2647
 *
 * @param mismatch_array, An int array that contains the position of mismatches
 * @param length, The length of the param `mismatch_array`
 * @return The local MIT score
 */
double calcMITLocalScore(uint64_t* mismatch_array, uint64_t length) {
    uint64_t i;
    double T1 = 1.0, T2, T3, d = 0.0, score;
    /* Mismatch penalty array */
    double M[] = { 0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583 };

    /* 1st term */
    for (i = 0; i < length; ++i)
        T1 = T1 * (1.0 - M[mismatch_array[i]]);

    /* 2nd term */
    if (length == 1)
        d = 19.0;
    else {
        for (i = 0; i < length - 1; ++i)
            d += mismatch_array[i + 1] - mismatch_array[i];
        d = d / (length - 1);
    }
    T2 = 1.0 / ((19.0 - d) / 19.0 * 4.0 + 1);

    /* 3rd term */
    T3 = 1.0 / (length * length);

    /* Local score */
    score = T1 * T2 * T3 * 100;
    return score;
}

/**
 * predictMITLocalScore
 *
 * Calculates mistmatch postions before calling scoring equation.
 *
 * @param xoredSignatures, A mask that represents the mismatches between a potential target and off-target. (The result of the XOR operation in practice)
 * @return The local MIT score
 */
double predictMITLocalScore(uint64_t xoredSignatures)
{
    std::vector<uint64_t> mismatch_array(20);
    uint64_t m = 0;
    for (uint64_t j = 0; j < seqLength; j++) {
        if ((xoredSignatures >> (j * 2)) & 0x3) {
            mismatch_array[m++] = j;
        }
    }
    if (m == 0) return 0.0;
    return calcMITLocalScore(mismatch_array.data(), m);
}

int main(int argc, char** argv)
{
    // Check number of args
    if (argc < 5)
    {
        std::cerr << fmt::format("Usage: {} [offtargetSites.txt] [sliceconfig.txt] [sequence length] [sissltable]\n", argv[0]) << std::endl;
        exit(1);
    }


    // Check seq length 
    seqLength = atoi(argv[3]);
    if (seqLength > 21)
    {
        std::cerr << "Sequence length is greater than 21, which is the maximum supported currently\n" << std::endl;
        exit(1);
    }


    // Check offtarget sites exists
    path otFile(argv[1]);
    if (!exists(otFile))
    {
        std::cerr << fmt::format("Could not find the specified offtarget sites file: {}", otFile.string()) << std::endl;
        exit(1);
    }


    // Chek offtarget sites file size
    uintmax_t otFileSize = file_size(otFile);
    uint64_t seqLineLength = seqLength + getLineEnding(otFile).size();
    if (otFileSize % seqLineLength != 0)
    {
        std::cerr << fmt::format("fileSize: {}\n", otFileSize);
        std::cerr << fmt::format("Error: offtargetSites.txt file does is not a multiple of the expected line length ({})\n", seqLineLength);
        std::cerr << "The sequence length may be incorrect; alternatively, the line endings\n";
        std::cerr << fmt::format("may be something other than {}, or there may be junk at the end of the file.", getLineEnding(otFile)) << std::endl;
        exit(1);
    }

    // Slice config file exists
    path scFile(argv[2]);
    if (!exists(scFile))
    {
        std::cerr << fmt::format("Could not find the specified slice config file: {}", scFile.string()) << std::endl;
        exit(1);
    }


    // Check slice config file size
    uintmax_t scFileSize = file_size(scFile);
    uint64_t scLineLength = seqLength + getLineEnding(scFile).size();
    if (scFileSize % scLineLength != 0) {
        std::cerr << fmt::format("fileSize: {}\n", scFileSize);
        std::cerr << fmt::format("Error: sliceconfig.txt file does is not a multiple of the expected line length ({})\n", scLineLength);
        std::cerr << "The sequence length may be incorrect; alternatively, the line endings\n";
        std::cerr << fmt::format("may be something other than {}, or there may be junk at the end of the file.", getLineEnding(scFile)) << std::endl;
        exit(1);
    }

    // Read in and genereate slice masks
    ifstream scInFile;
    scInFile.open(argv[2], std::ios::in | std::ios::binary);
    vector<vector<uint64_t>> sliceMasks;
    vector<uint64_t> sliceMasksBinary;
    for (string line; std::getline(scInFile, line);)
    {
        vector<uint64_t> mask;
        uint64_t maskBinary = 0ULL;
        for (uint64_t j = 0; j < seqLength; j++)
        {
            if (line[j] == '1')
            {
                maskBinary |= 1ULL << j;
                mask.push_back(j);
            }   
        }
        sliceMasks.push_back(mask);
        sliceMasksBinary.push_back(maskBinary);
    }
    scInFile.close();
    size_t sliceCount = sliceMasks.size();


    // Begin counting off targets
    uint64_t seqCount = otFileSize / seqLineLength;
    std::cout << fmt::format("Number of sequences: {}", seqCount) << std::endl;

    uint64_t globalCount = 0;
    uint64_t offtargetsCount = 0;
    vector<uint64_t> seqSignatures;
    vector<uint32_t> seqSignaturesOccurrences;

    // Read off targets into memory
    ifstream otInFile;
    otInFile.open(argv[1], std::ios::in | std::ios::binary);
    vector<char> entireDataSet(otFileSize);
    otInFile.read(entireDataSet.data(), otFileSize);
    otInFile.close();

    uint64_t progressCount = 0;
    uint64_t offtargetId = 0;

    std::cout << "Counting occurrences..." << std::endl;
    while (progressCount < seqCount) {
        char* ptr = &entireDataSet[progressCount * seqLineLength];
        uint64_t signature = sequenceToSignature(ptr);
        // check how many times the off-target appears
        // (assumed the list is sorted)
        uint32_t occurrences = 1;
        while (memcmp(ptr, ptr + (seqLineLength * occurrences), seqLength) == 0) {
            occurrences++;
            if ((seqCount - progressCount - occurrences) < 100)
                std::cout << fmt::format("{}/{} : {}", (progressCount + occurrences), seqCount, offtargetsCount) << std::endl;
        }

        seqSignatures.push_back(signature);
        seqSignaturesOccurrences.push_back(occurrences);
        offtargetsCount++;
        if (progressCount % 10000 == 0)
            std::cout << fmt::format("{}/{} : {}", progressCount, seqCount, offtargetsCount) << std::endl;
        progressCount += occurrences;
    }
    std::cout << "Finished!" << std::endl;

    std::cout << "Writing index header to file..." << std::endl;
    std::ofstream isslIndex;
    isslIndex.open(argv[4], std::ios::out | std::ios::binary);
    isslIndex.write(reinterpret_cast<char*>(&offtargetsCount), sizeof(uint64_t));
    isslIndex.write(reinterpret_cast<char*>(&seqLength), sizeof(uint64_t));
    isslIndex.write(reinterpret_cast<char*>(&sliceCount), sizeof(size_t));
    std::cout << "Finished!" << std::endl;

    std::cout << "Writing offtargets to file..." << std::endl;
    isslIndex.write(reinterpret_cast<char*>(seqSignatures.data()), sizeof(uint64_t) * seqSignatures.size());
    std::cout << "Finished!" << std::endl;

    std::cout << "Writing slice masks to file..." << std::endl;
    for (uint64_t& maskBinary : sliceMasksBinary)
    {
        isslIndex.write(reinterpret_cast<char*>(&maskBinary), sizeof(uint64_t));
    }
    std::cout << "Finished!" << std::endl;

    isslIndex.close();

    std::cout << "Constructing index..." << std::endl;
    for (size_t i = 0; i < sliceMasks.size(); i++)
    {
        std::cout << fmt::format("\tBuilding slice list {}", i+1) << std::endl;
        size_t sliceListSize = 1ULL << (sliceMasks[i].size() * 2);
        vector<vector<uint64_t>> sliceList(sliceListSize);
        vector<std::mutex> mtxList(sliceListSize);
        uint32_t signatureId = 0;

        #pragma omp parallel for
        for (int64_t s = 0; s < seqSignatures.size(); s++) {
            uint32_t occurrences = seqSignaturesOccurrences[signatureId];
            uint32_t sliceVal = 0ULL;
            for (size_t j = 0; j < sliceMasks[i].size(); j++)
            {
                sliceVal |= ((seqSignatures[s] >> (sliceMasks[i][j] * 2)) & 3ULL) << (j * 2);
            }
            // seqSigIdVal represnets the sequence signature ID and number of occurrences of the associated sequence.
            // (((uint64_t)occurrences) << 32), the most significant 32 bits is the count of the occurrences.
            // (uint64_t)signatureId, the index of the sequence in `seqSignatures`
            uint64_t seqSigIdVal = (static_cast<uint64_t>(occurrences) << 32) | static_cast<uint64_t>(s);
            mtxList[sliceVal].lock();
            sliceList[sliceVal].push_back(seqSigIdVal);
            mtxList[sliceVal].unlock();
        }

        std::cout << "\tFinished!" << std::endl;

        std::cout << fmt::format("\tWriting slice list {} to file...", i + 1) << std::endl;

        vector<size_t> sliceListLengths(sliceListSize);
        #pragma omp parallel for
        for (int64_t j = 0; j < sliceListSize; j++) { // Slice limit given slice width
            sliceListLengths[j] = sliceList[j].size();
        }

        isslIndex.open(argv[4], std::ios::out | std::ios::binary | std::ios::app);
        // Write slice list lengths
        isslIndex.write(reinterpret_cast<char*>(sliceListLengths.data()), sizeof(size_t) * sliceListSize);
        // write slice list data
        for (size_t j = 0; j < sliceListSize; j++) { // Slice limit given slice width
            isslIndex.write(reinterpret_cast<char*>(sliceList[j].data()), sizeof(uint64_t) * sliceList[j].size());
        }
        isslIndex.close();
        std::cout << "\tFinished!" << std::endl;
    }
    std::cout << "Finished!" << std::endl;
    std::cout << "Done" << std::endl;
    return 0;
}
