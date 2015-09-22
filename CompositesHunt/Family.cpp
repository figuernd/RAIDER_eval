#include "Family.h"

bool isAdjacent(Family* fam1, Family* fam2, Family* fam3) {
	//TODO
	return true;
}

Family::Family() {
	id = -1;
	cnt = 0;
	minL = 0;
	maxWgt = 0;
	reWgt = 0;
	status = 0;
	hasParent = false;
	nextNode = NULL;
	nxtNoneInterLen = 0;
}

//TODO
Family::Family(Family* otherFam) {
	id = otherFam -> id;
	cnt = otherFam -> cnt;
	minL = otherFam -> minL;
	maxWgt =otherFam -> maxWgt;
	status = 0;
	std::copy(otherFam->startPositions.begin(), otherFam->startPositions.end(), std::back_inserter(startPositions));
	nextNode = NULL;
	nxtNoneInterLen = 0;
	hasParent = otherFam->hasParent;
}

void Family::increaseCount(){
	cnt++;
}

int Family::getCount() const{
	return cnt;
}

int Family::getTotalLen() const{
	return cnt * minL;
}

void Family::setMinL(int L) {
	minL = L;
}

void Family::addStartPosition(int startPos){
	startPositions.push_back(startPos);
}

void Family::addPosition(int startPos, int endPos){
	startPositions.push_back(startPos);
	endPositions.push_back(endPos);
}

void Family::setID(int ID){
	this->id = ID;
}

int Family::getID() const{
	return id;
}


void Family::dfsHunt(Family* preNode) {
	if (status == 1) {
		preNode->descendants[this] = -1;
		return;
	} else if (status == 2) return;

	if (descendants.size() == 0) {
		status = 2;
		return;
	}
	status = 1;
	for (std::unordered_map<Family*, float>::iterator iter = descendants.begin(); iter != descendants.end(); iter++) {
		Family *child = iter->first;
		if(isAdjacent(preNode, this, child)){
			child->dfsHunt(this);
			if (this->descendants[child] < 0) continue; // to determine if loop exists ! 
			float wgt = child -> maxWgt;
			if (wgt + descendants[child] > maxWgt) {
				maxWgt = wgt + descendants[child];
				nextNode = child;
				nxtNoneInterLen = deNoneInterLen[child];
			}
		}
	}
	
	if (nextNode != NULL){
		nextNode->hasParent = true;
		// if ( std::find(nextNode->parents.begin(), nextNode->parents.end(), this) == nextNode->parents.end() ) 
		nextNode->parents.push_back(this);
	}
	status = 2;
}

void Family::dfsHuntRePath(Family* nxtNode){
	if (status == 1) {
		nxtNode->descendants[this] = -1;
		return;
	} else if (status == 2) return;
	
	if (parents.empty()) {
		status = 2;
		return;
	}
	status = 1;
	Family *maxPar = NULL;
	for (size_t famIndex = 0; famIndex < parents.size(); famIndex++) {
		Family *par = parents[famIndex];
		if (par->descendants[this] < 0) continue; // detect loop
		par->dfsHuntRePath(this);
		int wgt = par->reWgt + 1;
		if (wgt > reWgt) {
			reWgt = wgt;
			maxPar = par;
		}
	}
	
	for (size_t famIndex = 0; famIndex < parents.size(); famIndex++) {
		Family *par = parents[famIndex];
		// if (par->descendants[this] < 0 || (1.0 * reWgt / tmpMap[par] > 2) || ( (par != maxPar) && (1.0 * (reWgt - tmpMap[par]) / reWgt < 0.01) ) ) 
		// TODO // par->descendants[this] < 0 ||
		if (par->descendants[this] < 0 || par != maxPar) 	par->descendants[this] = 0;
	}
	status = 2;
}

void Family::clear() {
	maxWgt = 0.0;
}

bool Family::operator<(const Family& otherFam) const {
	return this->id < otherFam.id;
}

//  bool Family::operator=(const Family& otherFam) const {
// 	return this->id == otherFam.id;
// }