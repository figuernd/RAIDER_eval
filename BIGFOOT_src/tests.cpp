#include "gtest/gtest.h"
#include "raider.h"
#include <vector>
#include <stdexcept>


using namespace std;


TEST(Lmer, countSubsequent) {
	Lmer lmer(0x0);

	EXPECT_EQ(lmer.getCount(), 0);
	lmer.countSubsequent('N');

	EXPECT_EQ(lmer.getCount(), 1);
	ASSERT_THROW(lmer.getMaxSubsequentCount(), std::invalid_argument);
	ASSERT_THROW(lmer.getMaxSubsequentBase(), std::invalid_argument);
	ASSERT_THROW(lmer.countSubsequent('Z'), std::invalid_argument);

	EXPECT_NO_THROW({
		EXPECT_EQ(lmer.countSubsequent('G'), 2);
		EXPECT_EQ(lmer.getMaxSubsequentCount(), 1);
		EXPECT_EQ(lmer.getMaxSubsequentBase(), base_G);

		EXPECT_EQ(lmer.countSubsequent('T'), 3);
		EXPECT_EQ(lmer.countSubsequent('T'), 4);
		EXPECT_EQ(lmer.getMaxSubsequentCount(), 2);
		EXPECT_EQ(lmer.getMaxSubsequentBase(), base_T);

		EXPECT_EQ(lmer.countSubsequent('A'), 5);
		EXPECT_EQ(lmer.countSubsequent('A'), 6);
		EXPECT_EQ(lmer.countSubsequent('A'), 7);
		EXPECT_EQ(lmer.getMaxSubsequentCount(), 3);
		EXPECT_EQ(lmer.getMaxSubsequentBase(), base_A);

		EXPECT_EQ(lmer.countSubsequent('C'), 8);
		EXPECT_EQ(lmer.countSubsequent('C'), 9);
		EXPECT_EQ(lmer.countSubsequent('C'), 10);
		EXPECT_EQ(lmer.countSubsequent('C'), 11);
		EXPECT_EQ(lmer.countSubsequent('T'), 12);
		EXPECT_EQ(lmer.getMaxSubsequentCount(), 4);
		EXPECT_EQ(lmer.getMaxSubsequentBase(), base_C);

		EXPECT_EQ(lmer.getCount(), 12);
	});
}


TEST(Lmer, getMaxSubsequentBase) {
	Lmer lmer(0x0);

	lmer.countSubsequent('A');
	lmer.countSubsequent('C');
	lmer.countSubsequent('T');
	lmer.countSubsequent('A');
	lmer.countSubsequent('G');
	lmer.countSubsequent('G');
	lmer.countSubsequent('G');

	EXPECT_EQ(lmer.getMaxSubsequentBase(), base_G);
}


TEST(raider, baseToInt) {
	EXPECT_NO_THROW({
		EXPECT_EQ(baseToInt('A'), base_A);
		EXPECT_EQ(baseToInt('C'), base_C);
		EXPECT_EQ(baseToInt('G'), base_G);
		EXPECT_EQ(baseToInt('g'), base_G);
		EXPECT_EQ(baseToInt('T'), base_T);
	});
	ASSERT_THROW(baseToInt('Z'), std::out_of_range);
}


TEST(raider, seedToInt) {
	const char* seed = "ACGTAC";
	uint L = 6;
	size_t iSeed = seedToInt(seed, L);

	EXPECT_EQ(iSeed, 433);
}

TEST(raider, intToBase) {
	EXPECT_NO_THROW({
		EXPECT_EQ(intToBase(base_A), 'A');
		EXPECT_EQ(intToBase(base_C), 'C');
		EXPECT_EQ(intToBase(base_G), 'G');
		EXPECT_EQ(intToBase(base_T), 'T');
	});
	ASSERT_THROW(baseToInt(9), std::out_of_range);
}


TEST(raider, intToSeed) {
	char seed[6] = {0};
	size_t iSeed = 433;
	uint L = 6;

	intToSeed(iSeed, seed, L);

	EXPECT_STREQ(seed, "ACGTAC");
}


TEST(raider, getLmerMask) {
	uint L = 3;

	size_t expected = 63; // ..00111111
	EXPECT_EQ(getLmerMask(L), expected);

	L = 7;
	expected = 16383; // ..0011111111111111
	EXPECT_EQ(getLmerMask(L), expected);
}


TEST(raider, getNextSeedAndBase) {
	seqan::Dna5String sequence = "CATNTTCA";
	uint seqLength = seqan::length(sequence);
	vector<Lmer*> significantLmers;
	LmerMap lmerMap;
	uint L = 3;
	uint index = 0;

	char seed[4] = {0};
	char base;

	ASSERT_TRUE(getNextSeedAndBase(sequence, index, seqLength, L, seed, base));
	ASSERT_STREQ(seed, "CAT");
	ASSERT_EQ(base, 'N');
	ASSERT_EQ(index, 0);

	ASSERT_TRUE(getNextSeedAndBase(sequence, ++index, seqLength, L, seed, base));
	ASSERT_STREQ(seed, "TTC");
	ASSERT_EQ(base, 'A');
	ASSERT_EQ(index, 4);

	ASSERT_FALSE(getNextSeedAndBase(sequence, ++index, seqLength, L, seed, base));
}

TEST(raider, slide1Base) {
	size_t LMER_MASK = getLmerMask(3);
	size_t instance = seedToInt("ACA", 3);
	size_t nextBase = base_T;
	size_t result = slide1Base(instance, nextBase, LMER_MASK);

	char str[4] = {0};
	intToSeed(result, str, 3);

	EXPECT_STREQ(str, "CAT");
}


TEST(raider, countLmers) {
	seqan::Dna5String sequence = "CATGTTCATGTT";
	vector<Lmer*> significantLmers;
	LmerMap lmerMap;
	uint L = 3;
	uint M = 2;

	countLmers(sequence, significantLmers, lmerMap, L, M);

	ASSERT_EQ(significantLmers.size(), 3);

	char seed[3] = {0};
	intToSeed(significantLmers[0]->getFirstInstance(), seed, L);
	EXPECT_STREQ(seed, "CAT");
	intToSeed(significantLmers[1]->getFirstInstance(), seed, L);
	EXPECT_STREQ(seed, "ATG");
	intToSeed(significantLmers[2]->getFirstInstance(), seed, L);
	EXPECT_STREQ(seed, "TGT");

	size_t iSeed = seedToInt("CAT", L);
	EXPECT_EQ(lmerMap[iSeed], significantLmers[0]);
	iSeed = seedToInt("ATG", L);
	EXPECT_EQ(lmerMap[iSeed], significantLmers[1]);
	iSeed = seedToInt("TGT", L);
	EXPECT_EQ(lmerMap[iSeed], significantLmers[2]);
}

TEST(raider, countCompatible) {
	float M = 0.9;

	EXPECT_TRUE(countCompatible(100, 110, M));
	EXPECT_FALSE(countCompatible(100, 115, M));
	EXPECT_TRUE(countCompatible(110, 100, M));
	EXPECT_FALSE(countCompatible(115, 100, M));
}

TEST(raider, determineNextBase) {
	float T = 0.51;
	Lmer* lmer = new Lmer(seedToInt("CAT", 3));

	lmer->countSubsequent('G');
	EXPECT_EQ(determineNextBase(lmer, T), base_G);

	lmer->countSubsequent('A');
	EXPECT_EQ(determineNextBase(lmer, T), base_NULL);

}

TEST(raider, determineNextLmer) {
	seqan::Dna5String sequence = "CATGTTCATGTT";
	vector<Lmer*> lmers;
	LmerMap lmerMap;
	float M = 1;
	uint L = 3;
	float T = 0.5;
	const size_t LMER_MASK = getLmerMask(L);

	countLmers(sequence, lmers, lmerMap, L, M);

	Lmer* lmer = lmers[0];

	Lmer* next = determineNextLmer(lmer, lmerMap, M, T, LMER_MASK);

	EXPECT_EQ(next, lmers[1]);
}

TEST(raider, linkFamilies) {
	seqan::Dna5String sequence = "CATGTTCATGTT";
	vector<Lmer*> lmers;
	LmerMap lmerMap;
	uint L = 3;
	float M = 1;
	float T = 0.5;
	const size_t LMER_MASK = getLmerMask(L);

	countLmers(sequence, lmers, lmerMap, L, M);
	linkFamilies(lmers, lmerMap, M, T, LMER_MASK);

	EXPECT_EQ(lmers[0]->next, lmers[1]);
}

