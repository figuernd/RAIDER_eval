#ifndef RAIDER
#define RAIDER


#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <stdexcept>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <utility>
#include "Lmer.h"

using namespace std;


typedef unordered_map<size_t, Lmer*> LmerMap;


size_t baseToInt(char b) {
	switch (toupper(b)) {
	case 'A':
		return base_A;
	case 'C':
		return base_C;
	case 'G':
		return base_G;
	case 'T':
		return base_T;
	default:
		throw std::out_of_range("Unknown char base");
		return -1;
	}
}

size_t seedToInt(const char* seed, const uint L) {
	//TODO use bitset or vector<bool> to allow for potentially large bitsets
	size_t result = 0;
	for (uint i = 0; i < L; i++) {
		result <<= 2;
		result |= baseToInt(seed[i]);
	}
	return result;
}


char intToBase(size_t i) {
	switch (i) {
	case base_A:
		return 'A';
	case base_C:
		return 'C';
	case base_G:
		return 'G';
	case base_T:
		return 'T';
	default:
		throw std::out_of_range("Unknown int base");
		return -1;
	}
}

void intToSeed(size_t iSeed, char* seed, const uint L) {
	for (int i = L-1; i >= 0; i--) {
		seed[i] = intToBase(iSeed & 0b11);
		iSeed >>= 2;
	}
}

bool getNextSeedAndBase(seqan::Dna5String &sequence, uint &index, uint seqLength, uint L, char* seed, char &base) {
	for (uint i = 0; i < L && index + i < seqLength; i++) {
	  char b = toupper((char) sequence[index + i]);
	  if (b == 'N') {
	    index = index + i + 1;
	    i = -1;
	  } else {
	    seed[i] = sequence[index + i];
	  }
	}
	if (index+L < seqLength) {
	  base = sequence[index + L];
	  return true;
	}
	return false;
}

Lmer* getAndInsert(LmerMap &lmerMap, size_t seed) {
	Lmer*& p = lmerMap[seed];
	if (p == 0) {
		p = new Lmer(seed);
	}
	return p;
}



void countLmers(seqan::Dna5String &sequence, vector<Lmer*> &significantLmers, LmerMap &lmerMap, uint L, uint C) {
	const uint seqLength = seqan::length(sequence);
	uint index = -1;
	char unmasked[L];
	char base;

	while (getNextSeedAndBase(sequence, ++index, seqLength, L, unmasked, base)) {
		size_t seed = seedToInt(unmasked, L);
		Lmer* lmer = getAndInsert(lmerMap, seed);
		lmer->countSubsequent(base);
		if (lmer->getCount() == (int)C) {
			significantLmers.push_back(lmer);
		}
	}
}

uint determineNextBase(Lmer* lmer, float T) {
	size_t maxBase = lmer->getMaxSubsequentBase();
	unsigned short maxCount = lmer->getMaxSubsequentCount();
	if (((double)maxCount / (double)lmer->getCount()) >= T) {
		return maxBase;
	}
	return base_NULL;
}

size_t slide1Base(size_t instance, size_t nextBase, const size_t LMER_MASK) {
	instance <<= 2;
	instance |= nextBase;
	return instance &= LMER_MASK;
}

bool countCompatible(int count1, int count2, float I) {
	float min = (float) count1;
	float max = (float) count2;
	if (min > max)
		swap(min, max);
	return (min / max) >= I;
}

Lmer* determineNextLmer(Lmer* lmer, LmerMap &lmerMap, float I, float T, const size_t LMER_MASK) {
	size_t nextBase = determineNextBase(lmer, T);
	if (nextBase == base_NULL) {
		return nullptr;
	}
	size_t nextText = slide1Base(lmer->getFirstInstance(), nextBase, LMER_MASK);
	LmerMap::iterator it = lmerMap.find(nextText);
	if (it != lmerMap.end() && countCompatible(lmer->getCount(), it->second->getCount(), I)) {
		return it->second;
	}
	return nullptr;
}

void linkFamilies(vector<Lmer*> &lmers, LmerMap &lmerMap, float I, float T, const size_t LMER_MASK) {
	for (vector<Lmer*>::iterator it = lmers.begin(); it != lmers.end(); ++it) {
		Lmer* start = *it;
		if (start->state == Lmer::STATE_UNINITIALIZED) {
			start->state = Lmer::STATE_FAMILY_HEAD;
			Lmer* lmer = start;
			Lmer* next = determineNextLmer(lmer, lmerMap, I, T, LMER_MASK);
			while (next && next->state != Lmer::STATE_FAMILY_INTERIOR && next != start) {
				lmer->next = next;
				next->state = Lmer::STATE_FAMILY_INTERIOR;
				lmer = next;
				next = determineNextLmer(lmer, lmerMap, I, T, LMER_MASK);
			}
		}
	}
}


void assembleSeeds(vector<Lmer*> &significantLmers, uint L, vector<pair<string, int> > &results) {
	uint countFams = 0;
	uint countRepeats = 0;
	uint tableSize = significantLmers.size();

	for (uint i = 0; i < tableSize; i++) {
		Lmer* lmer = significantLmers[i];
		if (lmer->state == Lmer::STATE_FAMILY_HEAD) {
			countFams++;
			countRepeats += lmer->getCount();
			Lmer* nextLmer = lmer->next;
			std::stringstream ss;
			char seed[L + 1];
			seed[L] = 0;
			intToSeed(lmer->getFirstInstance(), seed, L);
			ss << seed;
			seed[1] = 0;
			while (nextLmer) {
				intToSeed(nextLmer->getFirstInstance(), seed, 1);
				ss << seed;
				nextLmer = nextLmer->next;
			}
			results.push_back(make_pair(ss.str(), lmer->getCount()));
		}
	}
}


/**
 * Creates an integer to logical AND to integer Lmers to ensure they
 * only contain an Lmer's worth of bases. Assumes 1 base occupies 2
 * bits of an integer. Thus if L=3, it will return 2*3=6 1's. For
 * example, if size_t were 8-bit and L=3, it would return an int
 * whose binary would be: 00111111. Logical AND'ing this to an
 * integer Lmer would ensure no data is kept in the first two bits.
 */
size_t getLmerMask(uint L) {
	return pow(2,2*L) - 1;
}

void raider(seqan::Dna5String &sequence, uint L, uint C, float I, float T, vector<pair<string, int> > &results) {
	const size_t LMER_MASK = getLmerMask(L);
	LmerMap lmerMap;
	vector<Lmer*> significantLmers;
	countLmers(sequence, significantLmers, lmerMap, L, C);
	linkFamilies(significantLmers, lmerMap, I, T, LMER_MASK);
	assembleSeeds(significantLmers, L, results);
}


#endif //RAIDER
