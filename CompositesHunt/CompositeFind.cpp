#include "CompositeFind.h"
#include "sequence.hh"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>


bool compareFam(const Family& f1, const Family& f2)
{
	return f1.getCount() > f2.getCount();
}

CompositeFind::CompositeFind(){

}

void CompositeFind::splitEleLine(const std::string& src, std::vector<std::string>& dest){
	dest.clear();
	std::string buf;
	std::stringstream ss(src); 
	while (ss >> buf)
		dest.push_back(buf);
}

void CompositeFind::readInFamilies(std::string eleFilePath) {
	std::ifstream fin(eleFilePath.c_str());    
	std::string eleLine("");
	getline(fin,eleLine);

	std::vector<std::string> eleStrs;
	Family curFam;
	int preFamilyId = -1;
	int curFamilyId;
	int startPos;
	int endPos;

	while(getline(fin,eleLine)) {
		splitEleLine(eleLine, eleStrs);
		curFamilyId = atoi(eleStrs[0].c_str());
		startPos = atoi(eleStrs[4].c_str());
		if (curFamilyId != preFamilyId ) {
			if (preFamilyId != -1)
				families.push_back(curFam);
			curFam = Family();
			curFam.setID(curFamilyId);
			endPos = atoi(eleStrs[5].c_str());
			curFam.setMinL(endPos - startPos + 1);
			preFamilyId = curFamilyId;
		}
		curFam.increaseCount();
		curFam.addStartPosition(startPos);
	}
	families.push_back(curFam);
	fin.close();

	sort(families.begin(), families.end(), compareFam);
	std::cout<< "families sorted finished, families size:"<< families.size() <<std::endl;

	for (int j=0; j<(int)families.size(); j++) {
		int cnt = families[j].getCount();
		for (int i =0; i< cnt; i++)
			familyMap[families[j].startPositions[i]] = &families[j];
	}

}

void CompositeFind::buildGraph(std::string sequenceFilePath, std::string outFilePath) {
	std::ofstream outFile(outFilePath);
	//construct adjFamMap
	std::vector<int> sPos;
	//std::unordered_map<tpl_Fam, int, key_hash, key_equal> adjFamMap;
	std::for_each(familyMap.begin(), familyMap.end(), [&sPos](const std::pair<int, Family*> p) {sPos.push_back(p.first);}) ;
	
	Family* curFm;
	std::vector<int> startPs;
	std::map<Family*, int> endFamCache;
	for (size_t famIndex = 0; famIndex < families.size(); famIndex++) {
		curFm = &families[famIndex];
		startPs = curFm->startPositions;
		int bp;
		int rightPos;
		for (size_t i=0; i < startPs.size(); i++) {
			bp = startPs[i];
			rightPos = bp + curFm->minL - 1;
			for (int tmpFmStartPos=bp+1; tmpFmStartPos<rightPos; tmpFmStartPos++) {
				if (familyMap.find(tmpFmStartPos) != familyMap.end() && familyMap[tmpFmStartPos] != curFm) {
					Family* tmpFamily = familyMap[tmpFmStartPos];
					if (bp + curFm->minL >= tmpFmStartPos + tmpFamily->minL) continue;
					// 
					int interLen =  bp + curFm->minL- tmpFmStartPos ;
					if ( endFamCache.find(tmpFamily) != endFamCache.end() ) 
						endFamCache[tmpFamily] += interLen;
					else {
						endFamCache[tmpFamily] = interLen; 
						curFm->deNoneInterLen[tmpFamily] = tmpFamily->minL - interLen;
					}
				} 
			}
		}
		//TODO
		// keep families with overlapping ratio > 15%
		for (std::map<Family*, int>::iterator iter = endFamCache.begin(); iter != endFamCache.end(); iter++) {
			float ratio = iter->second * 1.0 / (curFm->getTotalLen());
			if (ratio >= 0.15) 
				curFm->descendants[iter->first] = ratio;
		}
		endFamCache.clear();
		startPs.clear();
	}
	
	std::cout << "rightAdjFamMap constructed"<< std::endl;

	for (size_t famIndex = 0; famIndex < families.size(); famIndex++) { 
		try {
			Family *famOrigin = &families[famIndex];
			if (famOrigin->status == 0)  // && !famOrigin->hasParent
				famOrigin->dfsHunt(NULL);
		} catch (std::exception e) {
			std::cout << "famIndex : " << famIndex << std::endl;
			std::cout << e.what() << std::endl;
		}
	}

	std::cout << "tree construction completed" << std::endl;

	for (size_t famIndex = 0; famIndex < families.size(); famIndex++) {
		Family *fam = &families[famIndex];
		fam->status = 0;
	}

	for (size_t famIndex = 0; famIndex < families.size(); famIndex++) {
		Family *fam = &families[famIndex];
		if (fam->parents.size() > 1) {
			fam->dfsHuntRePath(NULL);  
		}
	}

	// TODO
	Sequence sequence(sequenceFilePath);	
	sequence.read();

	for (size_t famIndex = 0; famIndex < families.size(); famIndex++) {
		Family* curNode = &families[famIndex];
		if ( !curNode->hasParent ) { // && (!curNode->containedInOtherFam)
			std::vector<Family*> compositeFams;
			bool hasPath = true;
			while (curNode != NULL) {
				compositeFams.push_back(curNode);
				if (curNode -> nextNode != NULL) {
					if (curNode -> descendants[curNode -> nextNode] < 0) {
						compositeFams.push_back(curNode -> nextNode);
						break;
					}  else if (curNode -> descendants[curNode -> nextNode] == 0) {
						hasPath = false;
						break;
					}
				} 
				curNode = curNode -> nextNode;
			}

			//TODO
			if (hasPath && (compositeFams.size() > 3)) {
				int preNoneInterLen = 0;
				outFile << ">" << families[famIndex].id <<  std::endl;
				for (size_t i = 0; i < compositeFams.size(); i++) {
					Family* fam = compositeFams[i];
					if (i == 0) {
						//outFile << "startPos : " << fam->startPositions[0] << "  len : "<< fam->minL; 
						outFile << sequence.subsequence(fam->startPositions[0], fam->minL); 
					} else {
						//outFile << "startPos : " << fam->startPositions[0] + fam->minL - preNoneInterLen << "  len : "<< preNoneInterLen; 
						outFile << sequence.subsequence(fam->startPositions[0] + fam->minL - preNoneInterLen, preNoneInterLen);
					}
					preNoneInterLen = fam->nxtNoneInterLen;
				}
				outFile << std::endl;
			} 
			compositeFams.clear();      // compositeFams.swap(std::vector<Family*>()); 
		} 
	}
}
