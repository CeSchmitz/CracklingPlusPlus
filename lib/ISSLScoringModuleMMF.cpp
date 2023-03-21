#include "../include/ISSLScoringModuleMMF.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::unordered_map;
using namespace boost::interprocess;

const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };

ISSLScoringModuleMMF::ISSLScoringModuleMMF(const cracklingConfig& config) : specificityModule(config)
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

void ISSLScoringModuleMMF::run(std::vector<guideResults>& candidateGuides)
{
    if (!toolIsSelected)
    {
        cout << "Off-target scoring has been configured not to run. Skipping Off-target scoring" << endl;
        return;
    }

    cout << "Mapping ISSL Index." << endl;

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
    file_mapping isslFp(ISSLIndex.string().c_str(), read_only);
    mapped_region region(isslFp, read_only);
    void* inFileFp = region.get_address();

    /** The index contains a fixed-sized header
     *      - the number of unique off-targets in the index
     *      - the length of an off-target
     *      - the number of slices
     */
    size_t* headerPtr = static_cast<size_t*>(inFileFp);
    size_t offtargetsCount = *headerPtr++;
    size_t seqLength = *headerPtr++;
    size_t sliceCount = *headerPtr++;

    /** Load in all of the off-target sites */
    uint64_t* offtargetsPtr = reinterpret_cast<uint64_t*>(headerPtr);


    /** Read the slice masks and generate 2 bit masks */
    uint64_t* sliceMasksPtr = static_cast<uint64_t*>(offtargetsPtr + offtargetsCount);
    vector<vector<uint64_t>> sliceMasks;
    for (size_t i = 0; i < sliceCount; i++)
    {
        vector<uint64_t> mask;
        for (uint64_t j = 0; j < seqLength; j++)
        {
            if (*sliceMasksPtr & (1ULL << j))
            {
                mask.push_back(j);
            }
        }
        sliceMasks.push_back(mask);
        sliceMasksPtr++;
    }

    /** The contents of the slices. Stored by slice
    * Contains:
    *   - Size of each list within the slice stored contiguously
    *   - The contents of all the lists stored contiguously
    */
    vector<size_t*> allSlicelistSizes(sliceCount);
    vector<uint64_t*> allSliceSignatures(sliceCount);
    size_t* listSizePtr = static_cast<size_t*>(sliceMasksPtr);
    uint64_t* signaturePtr = static_cast<uint64_t*>(sliceMasksPtr);
    for (size_t i = 0; i < sliceCount; i++)
    {
        allSlicelistSizes[i] = listSizePtr;
        signaturePtr = static_cast<uint64_t*>(listSizePtr + (1ULL << (sliceMasks[i].size() * 2)));
        allSliceSignatures[i] = signaturePtr;
        listSizePtr = static_cast<size_t*>(signaturePtr + offtargetsCount);
    }

    cout << "ISSL Index Mapped." << endl;

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

    for (size_t i = 0; i < sliceCount; i++) {
        uint64_t* sliceList = allSliceSignatures[i];
        size_t sliceLimit = 1ULL << (sliceMasks[i].size() * 2);
        for (size_t j = 0; j < sliceLimit; j++) {
            sliceLists[i][j] = sliceList;
            sliceList += allSlicelistSizes[i][j];
        }
    }

    cout << "Beginning Off-target scoring." << endl;
    std::atomic_ullong testedCount = 0;
    std::atomic_ullong failedCount = 0;

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
            uint64_t searchSignature = sequenceToSignature(candidateGuides[guideIdx].seq, 20);

            /** Global scores */
            double totScoreMit = 0.0;
            double totScoreCfd = 0.0;

            double maximum_sum = (10000.0 - config.scoreThreshold * 100) / config.scoreThreshold;

            bool checkNextSlice = true;

            size_t sliceLimitOffset = 0;
            /** For each ISSL slice */
            for (size_t i = 0; i < sliceCount; i++) {
                vector<uint64_t>& sliceMask = sliceMasks[i];
                auto& sliceList = sliceLists[i];

                uint64_t searchSlice = 0ULL;
                for (size_t j = 0; j < sliceMask.size(); j++)
                {
                    searchSlice |= ((searchSignature >> (sliceMask[j] * 2)) & 3ULL) << (j * 2);
                }

                size_t signaturesInSlice = allSlicelistSizes[i][searchSlice];
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
                        uint64_t xoredSignatures = searchSignature ^ offtargetsPtr[signatureId];
                        uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAULL;
                        uint64_t oddBits = xoredSignatures & 0x5555555555555555ULL;
                        uint64_t mismatches = (evenBits >> 1) | oddBits;
                        uint64_t dist = popcount64(mismatches);

                        if (dist >= 0 && dist <= config.maxDist) {
                            // Begin calculating MIT score
                            if (calcMIT) {
                                if (dist > 0) {
                                    totScoreMit += precalculatedMITScores.at(mismatches) * (double)occurrences;
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
                                        uint64_t offtargetIdentityPos = offtargetsPtr[signatureId];
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
                            if (this->config.method == otScoreMethod::mitAndCfd) {
                                if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (this->config.method == otScoreMethod::mitOrCfd) {
                                if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (this->config.method == otScoreMethod::avgMitCfd) {
                                if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (this->config.method == otScoreMethod::mit) {
                                if (totScoreMit > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (this->config.method == otScoreMethod::cfd) {
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

    cout << fmt::format(comma_locale, "\t{:L} of {:L} failed here.", failedCount, testedCount) << endl;
}


uint64_t ISSLScoringModuleMMF::sequenceToSignature(const std::string& seq, uint64_t seqLen)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLen; j++) {
        signature |= static_cast<uint64_t>(nucleotideIndex[seq[j]]) << (j * 2);
    }
    return signature;
}

string ISSLScoringModuleMMF::signatureToSequence(uint64_t sig, uint64_t seqLen)
{
    string sequence = string(seqLen, ' ');
    for (uint64_t j = 0; j < seqLen; j++) {
        sequence[j] = signatureIndex[(sig >> (j * 2)) & 0x3];
    }
    return sequence;
}


bool ISSLScoringModuleMMF::processGuide(const guideResults& guide)
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
        if (guide.passedBowtie2 == CODE_REJECTED) {
            return false;
        }
    }

    // None of the failure conditions have been meet, return true
    return true;
}
