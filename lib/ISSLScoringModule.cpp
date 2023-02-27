#include "../include/ISSLScoringModule.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::unordered_map;
namespace fs = std::filesystem;

const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };

ISSLScoringModule::ISSLScoringModule(cracklingConfig config) : specificityModule(config)
{
    this->config = config.offTarget;
    this->ISSLIndex = config.input.offtargetSites;
    /** Scoring methods. To exit early:
     *      - only CFD must drop below `threshold`
     *      - only MIT must drop below `threshold`
     *      - both CFD and MIT must drop below `threshold`
     *      - CFD or MIT must drop below `threshold`
     *      - the average of CFD and MIT must below `threshold`
     */
    switch (this->config.method)
    {
    case(otScoreMethod::mitAndCfd):
        calcCFD = true;
        calcMIT = true;
        break;
    case(otScoreMethod::mitOrCfd):
        calcCFD = true;
        calcMIT = true;
        break;
    case(otScoreMethod::avgMitCfd):
        calcCFD = true;
        calcMIT = true;
        break;
    case(otScoreMethod::mit):
        calcCFD = false;
        calcMIT = true;
        break;
    case(otScoreMethod::cfd):
        calcCFD = true;
        calcMIT = false;
        break;
    default:
        throw std::runtime_error("Scoring method is not valid.");
        break;
    }
}

void ISSLScoringModule::run(std::vector<guideResults>& candidateGuides)
{
	if (!toolIsSelected)
	{
		cout << "Off-target scoring has been configured not to run. Skipping Off-target scoring" << endl;
		return;
	}

	cout << "Loading ISSL Index." << endl;

    /** Begin reading the binary encoded ISSL, structured as:
     *      - a header (5 items)
     *      - precalcuated local MIT scores
     *      - length of all the slices
     *      - the positions within a slice
     *      - all binary-encoded off-target sites
     */
    FILE* isslFp = fopen(ISSLIndex.string().c_str(), "rb");

    if (isslFp == NULL) {
        throw std::runtime_error("Error reading index: could not open file\n");
    }

    /** The index contains a fixed-sized header
     *      - the number of unique off-targets in the index
     *      - the length of an off-target
     *      - the total number of off-targets
     *      - the number of slices
     *      - the number of precalculated MIT scores
     */
    vector<size_t> slicelistHeader(50);

    if (fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: header invalid\n");
    }

    size_t offtargetsCount = slicelistHeader[0];
    size_t seqLength = slicelistHeader[1];
    size_t seqCount = slicelistHeader[2];
    size_t sliceCount = slicelistHeader[3];
    size_t scoresCount = slicelistHeader[4];

    /** Read in the precalculated MIT scores
     *      - `mask` is a 2-bit encoding of mismatch positions
     *          For example,
     *              00 01 01 00 01  indicates mismatches in positions 1, 3 and 4
     *
     *      - `score` is the local MIT score for this mismatch combination
     */
    phmap::flat_hash_map<uint64_t, double> precalculatedScores;
    for (int i = 0; i < scoresCount; i++) {
        uint64_t mask = 0;
        double score = 0.0;
        fread(&mask, sizeof(uint64_t), 1, isslFp);
        fread(&score, sizeof(double), 1, isslFp);
        precalculatedScores.insert(pair<uint64_t, double>(mask, score));
    }


    /**
    * Read the slice lengths from header
    */
    vector<size_t> sliceLens;
    for (int i = 0; i < sliceCount; i++)
    {
        size_t sliceLen;
        fread(&sliceLen, sizeof(size_t), 1, isslFp);
        sliceLens.push_back(sliceLen);
    }

    vector<vector<int>> sliceMasks;
    for (int i = 0; i < sliceCount; i++) {
        vector<int> mask;
        for (int j = 0; j < sliceLens[i]; j++)
        {
            int pos;
            fread(&pos, sizeof(int), 1, isslFp);
            mask.push_back(pos);
        }
        sliceMasks.push_back(mask);
    }

    /** Load in all of the off-target sites */
    vector<uint64_t> offtargets(offtargetsCount);
    if (fread(offtargets.data(), sizeof(uint64_t), offtargetsCount, isslFp) == 0) {
        throw std::runtime_error("Error reading index: loading off-target sequences failed\n");
    }

    /** The number of signatures embedded per slice
     *
     *      These counts are stored contiguously
     *
     */
    size_t sliceListCount = 0;
    for (int i = 0; i < sliceCount; i++)
    {
        sliceListCount += 1ULL << (sliceLens[i] * 2);
    }
    vector<size_t> allSlicelistSizes(sliceListCount);

    if (fread(allSlicelistSizes.data(), sizeof(size_t), allSlicelistSizes.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: reading slice list sizes failed\n");
    }

    /** The contents of the slices
     *
     *      Stored contiguously
     *
     *      Each signature (64-bit) is structured as:
     *          <occurrences 32-bit><off-target-id 32-bit>
     */
    vector<uint64_t> allSignatures(offtargetsCount * sliceCount);

    if (fread(allSignatures.data(), sizeof(uint64_t), allSignatures.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: reading slice contents failed\n");
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

    // IDEA: Split each slice into it's own file. Use memory mapped files to allow parallel access without reading into memory. Should work due to the simple binary data stored in these files.
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
    for (int i = 0; i < sliceCount; i++)
    {
        sliceLists[i] = vector<uint64_t*>(1ULL << (sliceLens[i] * 2));
    }

    uint64_t* offset = allSignatures.data();
    size_t sliceLimitOffset = 0;
    for (size_t i = 0; i < sliceCount; i++) {
        size_t sliceLimit = 1ULL << (sliceLens[i] * 2);
        for (size_t j = 0; j < sliceLimit; j++) {
            size_t idx = sliceLimitOffset + j;
            sliceLists[i][j] = offset;
            offset += allSlicelistSizes[idx];
        }
        sliceLimitOffset += sliceLimit;
    }


    cout << "Beginning Off-target scoring." << endl;

    uint64_t testedCount = 0;
    uint64_t failedCount = 0;

    /** Begin scoring */
    omp_set_num_threads(config.threads);
    #pragma omp parallel
    {
        vector<uint64_t> offtargetToggles(numOfftargetToggles);
        uint64_t* offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;
        /** For each candidate guide */
        // TODO: update to openMP > v2 (Use clang compiler)
        #pragma omp for
        for (int guideIdx = 0; guideIdx < candidateGuides.size(); guideIdx++) {

            // Run time filtering
            if (!processGuide(candidateGuides[guideIdx])) { continue; }

            // Encode seqeunce to search signature
            auto searchSignature = sequenceToSignature(candidateGuides[guideIdx].seq, candidateGuides[guideIdx].seq.length());

            /** Global scores */
            double totScoreMit = 0.0;
            double totScoreCfd = 0.0;

            int numOffTargetSitesScored = 0;
            double maximum_sum = (10000.0 - config.scoreThreshold * 100) / config.scoreThreshold;

            bool checkNextSlice = true;

            size_t sliceLimitOffset = 0;
            /** For each ISSL slice */
            for (size_t i = 0; i < sliceCount; i++) {
                vector<int>& sliceMask = sliceMasks[i];
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
                        numOffTargetSitesScored += occurrences;

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
                        int dist = popcount64(mismatches);

                        if (dist >= 0 && dist <= config.maxDist) {
                            // Begin calculating MIT score
                            if (calcMIT) {
                                if (dist > 0) {
                                    totScoreMit += precalculatedScores[mismatches] * (double)occurrences;
                                }
                            }

                            // Begin calculating CFD score
                            if (calcCFD) {
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
                                else if (dist > 0 && dist <= config.maxDist) {
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
                                        searchSigIdentityPos &= (3UL << (pos * 2));
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
                                        offtargetIdentityPos &= (3UL << (pos * 2));
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
                        }
                    }
                }
                if (!checkNextSlice)
                    break;
                sliceLimitOffset += 1ULL << (sliceLens[i] * 2);
            }
            // Finished processing off targets, update results
            switch (this->config.method)
            {
            case(otScoreMethod::mitAndCfd):
                candidateGuides[guideIdx].mitOfftargetscore = 10000.0 / (100.0 + totScoreMit);
                candidateGuides[guideIdx].cfdOfftargetscore = 10000.0 / (100.0 + totScoreCfd);
                if (candidateGuides.at(guideIdx).mitOfftargetscore < config.scoreThreshold && candidateGuides.at(guideIdx).cfdOfftargetscore < config.scoreThreshold)
                {
                    candidateGuides[guideIdx].passedOffTargetScore = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[guideIdx].passedOffTargetScore = CODE_ACCEPTED; }
                break;
            case(otScoreMethod::mitOrCfd):
                candidateGuides[guideIdx].mitOfftargetscore = 10000.0 / (100.0 + totScoreMit);
                candidateGuides[guideIdx].cfdOfftargetscore = 10000.0 / (100.0 + totScoreCfd);
                if (candidateGuides.at(guideIdx).mitOfftargetscore < config.scoreThreshold || candidateGuides.at(guideIdx).cfdOfftargetscore < config.scoreThreshold)
                {
                    candidateGuides[guideIdx].passedOffTargetScore = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[guideIdx].passedOffTargetScore = CODE_ACCEPTED; }
                break;
            case(otScoreMethod::avgMitCfd):
                candidateGuides[guideIdx].mitOfftargetscore = 10000.0 / (100.0 + totScoreMit);
                candidateGuides[guideIdx].cfdOfftargetscore = 10000.0 / (100.0 + totScoreCfd);
                if ((candidateGuides.at(guideIdx).mitOfftargetscore + candidateGuides.at(guideIdx).cfdOfftargetscore) / 2 < config.scoreThreshold)
                {
                    candidateGuides[guideIdx].passedOffTargetScore = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[guideIdx].passedOffTargetScore = CODE_ACCEPTED; }
                break;
            case(otScoreMethod::mit):
                candidateGuides[guideIdx].mitOfftargetscore = 10000.0 / (100.0 + totScoreMit);
                candidateGuides[guideIdx].cfdOfftargetscore = -1;
                if (candidateGuides.at(guideIdx).mitOfftargetscore < config.scoreThreshold) 
                {
                    candidateGuides[guideIdx].passedOffTargetScore = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[guideIdx].passedOffTargetScore = CODE_ACCEPTED; }
                break;
            case(otScoreMethod::cfd):
                candidateGuides[guideIdx].mitOfftargetscore = -1;
                candidateGuides[guideIdx].cfdOfftargetscore = 10000.0 / (100.0 + totScoreCfd);
                if (candidateGuides.at(guideIdx).cfdOfftargetscore < config.scoreThreshold)
                {
                    candidateGuides[guideIdx].passedOffTargetScore = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[guideIdx].passedOffTargetScore = CODE_ACCEPTED; }
                break;
            default:
                throw std::runtime_error("Scoring method is not valid.");
                break;
            }
            testedCount++;
            memset(offtargetToggles.data(), 0, sizeof(uint64_t) * offtargetToggles.size());
        }
    }

    cout << fmt::format("\t{} of {} failed here.", failedCount, testedCount) << endl;
}


uint64_t ISSLScoringModule::sequenceToSignature(const std::string& seq, uint64_t seqLen)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLen; j++) {
        signature |= (uint64_t)(nucleotideIndex[seq[j]]) << (j * 2);
    }
    return signature;
}

string ISSLScoringModule::signatureToSequence(uint64_t sig, uint64_t seqLen)
{
    string sequence = string(seqLen, ' ');
    for (uint64_t j = 0; j < seqLen; j++) {
        sequence[j] = signatureIndex[(sig >> (j * 2)) & 0x3];
    }
    return sequence;
}


bool ISSLScoringModule::processGuide(const guideResults& guide)
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
		// Reject if the consensus was not passed
		if (guide.consensusCount < consensusN) { return false; }
		// Reject if bowtie2 was not passed
		if (guide.passedBowtie2 == CODE_REJECTED) { return false; }
	}

	// None of the failure conditions have been meet, return true
	return true;
}
