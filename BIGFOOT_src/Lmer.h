#ifndef RAIDER_LMER
#define RAIDER_LMER

#include <vector>
#include <iostream>
#include <assert.h>
#include <stdexcept>

#define base_A 0
#define base_C 1
#define base_G 2
#define base_T 3
#define base_NULL 4

using namespace std;

/**
 * Primary responsibility of this class is to represent an Lmer and keep count
 * of all bases that immediately follow it. Hence subsequentA, subsequentC, etc.
 * represent the counts of subsequent bases immediately following instances of
 * this Lmer, whose first instance is kept in firstInstance.
 */
class Lmer {
private:
	unsigned short subsequentA = 0;
	unsigned short subsequentC = 0;
	unsigned short subsequentG = 0;
	unsigned short subsequentT = 0;
	unsigned short *maxSubsequentBase = nullptr;

	size_t firstInstance;
	int countInstances = 0;

public:
	static const int STATE_UNINITIALIZED = 0;
	static const int STATE_FAMILY_HEAD = 1;
	static const int STATE_FAMILY_INTERIOR = 2;

	unsigned short state = STATE_UNINITIALIZED;
	Lmer* next = nullptr;

	Lmer(size_t firstInstance) {
		this->firstInstance = firstInstance;
	}


	uint countSubsequent(char base) {
		unsigned short *updated = nullptr;
		switch(base) {
		case 'A':
			subsequentA++;
			updated = &subsequentA;
			break;
		case 'C':
			subsequentC++;
			updated = &subsequentC;
			break;
		case 'G':
			subsequentG++;
			updated = &subsequentG;
			break;
		case 'T':
			subsequentT++;
			updated = &subsequentT;
			break;
		case 'N':
			break;
		default:
			throw std::invalid_argument("Invalid base (only A,C,G,T accepted)");
		}
		if (updated && (!maxSubsequentBase || *updated > *maxSubsequentBase)) {
			maxSubsequentBase = updated;
		}
		return ++countInstances;
	}


	size_t getFirstInstance() {
		return firstInstance;
	}


	int getCount() {
		return countInstances;
	}


	size_t getMaxSubsequentCount() {
		if (maxSubsequentBase == nullptr) {
			throw std::invalid_argument("Unable to get count - max subsequent pointer uninitialized");
		}
		return *maxSubsequentBase;
	}

	size_t getMaxSubsequentBase() {
		if (maxSubsequentBase == nullptr) {
			throw std::invalid_argument("Unable to get base - max subsequent pointer uninitialized");
		}
		if (maxSubsequentBase == &subsequentA)
			return base_A;
		if (maxSubsequentBase == &subsequentC)
			return base_C;
		if (maxSubsequentBase == &subsequentG)
			return base_G;
		if (maxSubsequentBase == &subsequentT)
			return base_T;
		throw std::out_of_range("Unknown base pointer");
	}


};

#endif //RAIDER_LMER
