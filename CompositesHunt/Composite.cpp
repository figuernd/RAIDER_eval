#include "Composite.h"

Composite::Composite() {

}

void Composite::append(Family fam){
	families.push_back(fam);
	famIdSet.insert(fam.getID());
} 

void Composite::insert2Front(Family fam){
	families.insert(families.begin(), fam);
	famIdSet.insert(fam.getID());
}

bool Composite::hasFamily(Family fam){
	return famIdSet.find(fam.getID()) != famIdSet.end() ? true:false;
}

Family* Composite::getLastFamily(){
	int famSize = families.size();
	return &families[famSize-1];
}
