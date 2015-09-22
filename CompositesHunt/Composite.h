#include "Family.h"
#include <set>

#ifndef LIU_RAIDER_COMPOSITE
#define LIU_RAIDER_COMPOSITE

class Composite
{
public:
	Composite();
	void append(Family fam);
	void insert2Front(Family fam);
	bool hasFamily(Family fam);
	Family* getLastFamily();
	std::vector<Family> families;

private:
	std::set<int> famIdSet;
};

#endif
