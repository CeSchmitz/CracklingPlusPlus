#include "ISSLScoreOfftargets.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::unordered_map;

// Char to binary encoding
const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
// Binary to char encoding
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };

uint64_t sequenceToSignature(const std::string& seq, uint64_t seqLen)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLen; j++) {
        signature |= static_cast<uint64_t>(nucleotideIndex[seq[j]]) << (j * 2);
    }
    return signature;
}

string signatureToSequence(uint64_t sig, uint64_t seqLen)
{
    string sequence = string(seqLen, ' ');
    for (uint64_t j = 0; j < seqLen; j++) {
        sequence[j] = signatureIndex[(sig >> (j * 2)) & 0x3];
    }
    return sequence;
}

int main(int argc, char** argv)
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s [issltable] [query file] [max distance] [score-threshold] [score-method]\n", argv[0]);
        exit(1);
    }

    /** The maximum number of mismatches */
    int maxDist = atoi(argv[3]);

    /** The threshold used to exit scoring early */
    double threshold = atof(argv[4]);

    /** Scoring methods. To exit early:
     *      - only CFD must drop below `threshold`
     *      - only MIT must drop below `threshold`
     *      - both CFD and MIT must drop below `threshold`
     *      - CFD or MIT must drop below `threshold`
     *      - the average of CFD and MIT must below `threshold`
     */
    string argScoreMethod = argv[5];
    otScoreMethod scoreMethod;
    bool calcCfd = false;
    bool calcMit = false;
    if (!argScoreMethod.compare("and")) {
        scoreMethod = otScoreMethod::mitAndCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("or")) {
        scoreMethod = otScoreMethod::mitOrCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("avg")) {
        scoreMethod = otScoreMethod::avgMitCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("mit")) {
        scoreMethod = otScoreMethod::mit;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("cfd")) {
        scoreMethod = otScoreMethod::cfd;
        calcCfd = true;
    }
    else
    {
        fprintf(stderr, "Invalid scoring method. Acceptable options are: 'and', 'or', 'avg', 'mit', 'cfd'", argv[0]);
        exit(1);
    }

    /** Begin reading the binary encoded ISSL, structured as:
     *  - The header (3 items)
     *  - All binary-encoded off-target sites
     *  - Slice masks
     *  - Size of slice 1 lists
     *  - Contents of slice 1 lists
     *  ...
     *  - Size of slice N lists (N being the number of slices)
     *  - Contents of slice N lists
     */
    FILE* isslFp = fopen(argv[1], "rb");

    if (isslFp == NULL) {
        throw std::runtime_error("Error reading index: could not open file\n");
    }

    /** The index contains a fixed-sized header
     *      - the number of unique off-targets in the index
     *      - the length of an off-target
     *      - the number of slices
     */
    vector<size_t> slicelistHeader(3);

    if (fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: header invalid\n");
    }

    size_t offtargetsCount = slicelistHeader[0];
    size_t seqLength = slicelistHeader[1];
    size_t sliceCount = slicelistHeader[2];

    /** Load in all of the off-target sites */
    vector<uint64_t> offtargets(offtargetsCount);
    if (fread(offtargets.data(), sizeof(uint64_t), offtargetsCount, isslFp) == 0) {
        throw std::runtime_error("Error reading index: loading off-target sequences failed\n");
    }

    /** Read the slice masks and generate 2 bit masks */
    vector<vector<uint64_t>> sliceMasks;
    for (size_t i = 0; i < sliceCount; i++)
    {
        uint64_t maskBinary;
        fread(&maskBinary, sizeof(uint64_t), 1, isslFp);

        vector<uint64_t> mask;
        for (uint64_t j = 0; j < seqLength; j++)
        {
            if (maskBinary & (1ULL << j))
            {
                mask.push_back(j);
            }
        }
        sliceMasks.push_back(mask);
    }

    /** Calculate the total number of lists based on slice count and width */
    size_t sliceListCount = 0;
    for (size_t i = 0; i < sliceCount; i++)
    {
        sliceListCount += 1ULL << (sliceMasks[i].size() * 2);
    }

    /** The number of signatures embedded per slice. Store continguously */
    vector<size_t> allSlicelistSizes(sliceListCount);

    /** The contents of the slices. Stored contiguously
     *  Each signature (64-bit) is structured as:
     *      <occurrences 32-bit><off-target-id 32-bit>
     */
    vector<uint64_t> allSignatures(offtargetsCount * sliceCount);

    /** The number of signatures embedded per slice. Store continguously */
    sliceListCount = 0;
    for (size_t i = 0; i < sliceCount; i++)
    {
        size_t sliceListSize = 1ULL << (sliceMasks[i].size() * 2);
        if (fread(allSlicelistSizes.data() + sliceListCount, sizeof(size_t), sliceListSize, isslFp) == 0) {
            throw std::runtime_error("Error reading index: reading slice list sizes failed\n");
        }

        if (fread(allSignatures.data() + (offtargetsCount * i), sizeof(uint64_t), offtargetsCount, isslFp) == 0) {
            throw std::runtime_error("Error reading index: reading slice contents failed\n");
        }

        sliceListCount += 1ULL << (sliceMasks[i].size() * 2);
    }

    /** End reading the index */
    fclose(isslFp);

    cout << "ISSL Index Loaded." << endl;



    /** Prevent assessing an off-target site for multiple slices
     *
     *      Create enough 1-bit "seen" flags for the off-targets
     *      We only want to score a candidate guide against an off-target once.
     *      The least-significant bit represents the first off-target
     *      0 0 0 1   0 1 0 0   would indicate that the 3rd and 5th off-target have been seen.
     *      The CHAR_BIT macro tells us how many bits are in a byte (C++ >= 8 bits per byte)
     */
    uint64_t numOfftargetToggles = (offtargetsCount / ((size_t)sizeof(uint64_t) * (size_t)CHAR_BIT)) + 1;

    /** Start constructing index in memory
     *
     *      To begin, reverse the contiguous storage of the slices,
     *         into the following:
     *
     *         + Slice 0 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |----  ...
     *         |
     *         + Slice 1 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |---- ...
     *         | ...
     */

    vector<vector<uint64_t*>> sliceLists(sliceCount);
    // Assign sliceLists size based on each slice length
    for (size_t i = 0; i < sliceCount; i++)
    {
        sliceLists[i] = vector<uint64_t*>(1ULL << (sliceMasks[i].size() * 2));
    }

    uint64_t* offset = allSignatures.data();
    size_t sliceLimitOffset = 0;
    for (size_t i = 0; i < sliceCount; i++) {
        size_t sliceLimit = 1ULL << (sliceMasks[i].size() * 2);
        for (size_t j = 0; j < sliceLimit; j++) {
            size_t idx = sliceLimitOffset + j;
            sliceLists[i][j] = offset;
            offset += allSlicelistSizes[idx];
        }
        sliceLimitOffset += sliceLimit;
    }

    //TODO: rewrite
    /** Load query file (candidate guides)
     *      and prepare memory for calculated global scores
     */
    size_t seqLineLength = seqLength + 1;
    std::filesystem::path queryFile(argv[2]);
    size_t fileSize = std::filesystem::file_size(queryFile);
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "Error: query file is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t queryCount = fileSize / seqLineLength;
    FILE* fp = fopen(argv[2], "rb");
    vector<char> queryDataSet(fileSize);
    vector<uint64_t> querySignatures(queryCount);
    vector<double> querySignatureMitScores(queryCount);
    vector<double> querySignatureCfdScores(queryCount);

    if (fread(queryDataSet.data(), fileSize, 1, fp) < 1) {
        fprintf(stderr, "Failed to read in query file.\n");
        exit(1);
    }
    fclose(fp);

    /** Binary encode query sequences */
    #pragma omp parallel
    {
    #pragma omp for
        for (int i = 0; i < queryCount; i++) {
            char* ptr = &queryDataSet[i * seqLineLength];
            uint64_t signature = sequenceToSignature(ptr, 20);
            querySignatures[i] = signature;
        }
    }

    /** Begin scoring */
    omp_set_num_threads(32);
    #pragma omp parallel
    {
        vector<uint64_t> offtargetToggles(numOfftargetToggles);
        uint64_t* offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;
        /** For each candidate guide */
        // TODO: update to openMP > v2 (Use clang compiler)
        #pragma omp for
        for (int searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {

            auto searchSignature = querySignatures[searchIdx];

            /** Global scores */
            double totScoreMit = 0.0;
            double totScoreCfd = 0.0;

            double maximum_sum = (10000.0 - threshold * 100) / threshold;
            bool checkNextSlice = true;

            size_t sliceLimitOffset = 0;
            /** For each ISSL slice */
            for (size_t i = 0; i < sliceCount; i++) {
                vector<uint64_t>& sliceMask = sliceMasks[i];
                auto& sliceList = sliceLists[i];

                uint64_t searchSlice = 0ULL;
                for (int j = 0; j < sliceMask.size(); j++)
                {
                    searchSlice |= ((searchSignature >> (sliceMask[j] * 2)) & 3ULL) << (j * 2);
                }

                size_t idx = sliceLimitOffset + searchSlice;

                size_t signaturesInSlice = allSlicelistSizes[idx];
                uint64_t* sliceOffset = sliceList[searchSlice];

                /** For each off-target signature in slice */
                for (size_t j = 0; j < signaturesInSlice; j++) {
                    auto signatureWithOccurrencesAndId = sliceOffset[j];
                    auto signatureId = signatureWithOccurrencesAndId & 0xFFFFFFFFULL;
                    uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));

                    /** Prevent assessing the same off-target for multiple slices */
                    uint64_t seenOfftargetAlready = 0;
                    uint64_t* ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
                    seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;

                    if (!seenOfftargetAlready) {
                        *ptrOfftargetFlag |= (1ULL << (signatureId % 64));

                        /** Find the positions of mismatches
                            *
                            *  Search signature (SS):    A  A  T  T    G  C  A  T
                            *                           00 00 11 11   10 01 00 11
                            *
                            *        Off-target (OT):    A  T  A  T    C  G  A  T
                            *                           00 11 00 11   01 10 00 11
                            *
                            *                SS ^ OT:   00 00 11 11   10 01 00 11
                            *                         ^ 00 11 00 11   01 10 00 11
                            *                  (XORd) = 00 11 11 00   11 11 00 00
                            *
                            *        XORd & evenBits:   00 11 11 00   11 11 00 00
                            *                         & 10 10 10 10   10 10 10 10
                            *                   (eX)  = 00 10 10 00   10 10 00 00
                            *
                            *         XORd & oddBits:   00 11 11 00   11 11 00 00
                            *                         & 01 01 01 01   01 01 01 01
                            *                   (oX)  = 00 01 01 00   01 01 00 00
                            *
                            *         (eX >> 1) | oX:   00 01 01 00   01 01 00 00 (>>1)
                            *                         | 00 01 01 00   01 01 00 00
                            *            mismatches   = 00 01 01 00   01 01 00 00
                            *
                            *   popcount(mismatches):   4
                            */
                        uint64_t xoredSignatures = searchSignature ^ offtargets[signatureId];
                        uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAULL;
                        uint64_t oddBits = xoredSignatures & 0x5555555555555555ULL;
                        uint64_t mismatches = (evenBits >> 1) | oddBits;
                        uint64_t dist = popcount64(mismatches);

                        if (dist >= 0 && dist <= 4) {
                            // Begin calculating MIT score
                            if (calcMit) {
                                if (dist > 0) {
                                    totScoreMit += precalculatedMITScores.at(mismatches) * (double)occurrences;
                                }
                            }

                            // Begin calculating CFD score
                            if (calcCfd) {
                                /** "In other words, for the CFD score, a value of 0
                                    *      indicates no predicted off-target activity whereas
                                    *      a value of 1 indicates a perfect match"
                                    *      John Doench, 2016.
                                    *      https://www.nature.com/articles/nbt.3437
                                */
                                double cfdScore = 0;
                                if (dist == 0) {
                                    cfdScore = 1;
                                }
                                else {
                                    cfdScore = cfdPamPenalties[0b1010]; // PAM: NGG, TODO: do not hard-code the PAM

                                    for (size_t pos = 0; pos < 20; pos++) {
                                        size_t mask = pos << 4;

                                        /** Create the mask to look up the position-identity score
                                            *      In Python... c2b is char to bit
                                            *       mask = pos << 4
                                            *       mask |= c2b[sgRNA[pos]] << 2
                                            *       mask |= c2b[revcom(offTaret[pos])]
                                            *
                                            *      Find identity at `pos` for search signature
                                            *      example: find identity in pos=2
                                            *       Recall ISSL is inverted, hence:
                                            *                   3'-  T  G  C  C  G  A -5'
                                            *       start           11 10 01 01 10 00
                                            *       3UL << pos*2    00 00 00 11 00 00
                                            *       and             00 00 00 01 00 00
                                            *       shift           00 00 00 00 01 00
                                            */
                                        uint64_t searchSigIdentityPos = searchSignature;
                                        searchSigIdentityPos &= (3ULL << (pos * 2));
                                        searchSigIdentityPos = searchSigIdentityPos >> (pos * 2);
                                        searchSigIdentityPos = searchSigIdentityPos << 2;

                                        /** Find identity at `pos` for offtarget
                                            *      Example: find identity in pos=2
                                            *      Recall ISSL is inverted, hence:
                                            *                  3'-  T  G  C  C  G  A -5'
                                            *      start           11 10 01 01 10 00
                                            *      3UL<<pos*2      00 00 00 11 00 00
                                            *      and             00 00 00 01 00 00
                                            *      shift           00 00 00 00 00 01
                                            *      rev comp 3UL    00 00 00 00 00 10 (done below)
                                            */
                                        uint64_t offtargetIdentityPos = offtargets[signatureId];
                                        offtargetIdentityPos &= (3ULL << (pos * 2));
                                        offtargetIdentityPos = offtargetIdentityPos >> (pos * 2);

                                        /** Complete the mask
                                            *      reverse complement (^3UL) `offtargetIdentityPos` here
                                            */
                                        mask = (mask | searchSigIdentityPos | (offtargetIdentityPos ^ 3UL));

                                        if (searchSigIdentityPos >> 2 != offtargetIdentityPos) {
                                            cfdScore *= cfdPosPenalties[mask];
                                        }
                                    }
                                }
                                totScoreCfd += cfdScore * (double)occurrences;
                            }

                            /** Stop calculating global score early if possible */
                            if (scoreMethod == otScoreMethod::mitAndCfd) {
                                if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (scoreMethod == otScoreMethod::mitOrCfd) {
                                if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (scoreMethod == otScoreMethod::avgMitCfd) {
                                if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (scoreMethod == otScoreMethod::mit) {
                                if (totScoreMit > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (scoreMethod == otScoreMethod::cfd) {
                                if (totScoreCfd > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (!checkNextSlice)
                    break;
                sliceLimitOffset += 1ULL << (sliceMasks[i].size() * 2);
            }
            querySignatureMitScores[searchIdx] = 10000.0 / (100.0 + totScoreMit);
            querySignatureCfdScores[searchIdx] = 10000.0 / (100.0 + totScoreCfd);

            memset(offtargetToggles.data(), 0, sizeof(uint64_t) * offtargetToggles.size());
        }
    }

    /** Print global scores to stdout */
    for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
        auto querySequence = signatureToSequence(querySignatures[searchIdx], 20);
        printf("%s\t", querySequence.c_str());
        if (calcMit)
            printf("%f\t", querySignatureMitScores[searchIdx]);
        else
            printf("-1\t");

        if (calcCfd)
            printf("%f\n", querySignatureCfdScores[searchIdx]);
        else
            printf("-1\n");

    }

    return 0;
}
