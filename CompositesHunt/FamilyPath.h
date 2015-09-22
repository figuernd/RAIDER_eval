#include <vector>
#include "Family.h"

#ifndef LIU_RAIDER_FAMILY_PATH
#define LIU_RAIDER_FAMILY_PATH

class FamilyPath
{
public:
	FamilyPath(){
		path.clear();
		weight = -1;
	}
	float weight;
	std::vector<Family*> path;
};

#endif